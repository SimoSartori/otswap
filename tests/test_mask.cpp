// SPDX-License-Identifier: GPL-2.0-or-later
/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * This program is free software; you can redistribute it and/or    *
 * modify it under the terms of the GNU General Public License as   *
 * published by the Free Software Foundation; either version 2 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU General Public License for more details.                     *
 *                                                                  *
 * You should have received a copy of the GNU General Public        *
 * License along with this program; if not, write to the Free       *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,      *
 * Boston, MA 02110-1301 USA.                                       *
 ********************************************************************/

/**
 *  @file tests/test_mask.cpp
 *
 *  @brief The Healpix veto mask, and the rejection of displacements whose
 *  great-circle path crosses masked sky.
 */

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <fitsio.h>
#include <arr.h>
#include <healpix_base.h>
#include <healpix_map.h>
#include <pointing.h>
#include <rangeset.h>
#include <vec3.h>

#include "arc.h"
#include "internal.h"
#include "check.h"

using namespace otswap;

// CosmoBolognaLib's arc_crosses_mask and what it calls, copied verbatim from
// the anonymous namespace at the top of OTreconstruction/OTreconstruction.cpp:
// the reference the port is checked against. Like otswap, it takes a pixel
// as observed when its value is > 0. It finds pixels with Healpix's own
// vec2pix, over the system atan2, where otswap uses its deterministic one;
// the two can differ only for a direction within rounding of a pixel
// boundary.
namespace cbl_reference {

  constexpr double kArcFloorInPixels = 1.e-6;


  bool edge_crossing_safe (const Healpix_Map<float>& mask, const int pix_p, const int pix_q,
                           const int pix_a, const int pix_b) {
    fix_arr<int, 8> nb_p, nb_q;
    mask.neighbors(pix_p, nb_p);
    if (nb_p[0] != pix_q && nb_p[2] != pix_q && nb_p[4] != pix_q && nb_p[6] != pix_q) return false;

    mask.neighbors(pix_q, nb_q);
    for (int i = 0; i < 8; ++i) {
      const int r = nb_p[i];
      if (r < 0 || r == pix_q || r == pix_a || r == pix_b || mask[r] > 0.) continue;
      for (int j = 0; j < 8; ++j)
        if (nb_q[j] == r) return false;
    }
    return true;
  }


  bool arc_segment_crosses_mask (const Healpix_Map<float>& mask,
                                 const vec3& p, const int pix_p, const vec3& q, const int pix_q,
                                 const int pix_a, const int pix_b, const double floor2, const double pixel2) {
    if (pix_p == pix_q) return false;
    const double d2 = (p - q).SquaredLength();
    if (d2 < floor2) return false;
    if (d2 < pixel2 && edge_crossing_safe(mask, pix_p, pix_q, pix_a, pix_b)) return false;

    const vec3 m = (p + q).Norm();
    const int pix_m = mask.vec2pix(m);
    if (pix_m != pix_a && pix_m != pix_b && !(mask[pix_m] > 0.)) return true;

    return arc_segment_crosses_mask(mask, p, pix_p, m, pix_m, pix_a, pix_b, floor2, pixel2)
        || arc_segment_crosses_mask(mask, m, pix_m, q, pix_q, pix_a, pix_b, floor2, pixel2);
  }


  bool arc_crosses_mask (const Healpix_Map<float>& mask,
                         const double x0, const double y0, const double z0,
                         const double x1, const double y1, const double z1,
                         const double pixel_size) {
    const double r0 = std::sqrt(x0*x0 + y0*y0 + z0*z0);
    const double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
    if (r0 == 0. || r1 == 0.) return false;

    const vec3 a(x0/r0, y0/r0, z0/r0);
    const vec3 b(x1/r1, y1/r1, z1/r1);
    const int pix_a = mask.vec2pix(a);
    const int pix_b = mask.vec2pix(b);
    const double arc_floor = kArcFloorInPixels * pixel_size;
    return arc_segment_crosses_mask(mask, a, pix_a, b, pix_b, pix_a, pix_b, arc_floor*arc_floor, pixel_size*pixel_size);
  }

}

namespace {

  constexpr double kPi = 3.14159265358979323846;

  // NSIDE 16 gives 3072 pixels, exactly three rows of the 1024-per-row
  // layout the reader expects.
  constexpr int kNside = 16;

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "otswap_test_" + name;
  }

  // Writes a Healpix map in the layout the reader assumes: a binary table
  // in HDU 2, one vector column of 1024 floats per row, with NSIDE and
  // ORDERING in that HDU's header.
  long write_map (const std::string& file, const std::vector<float>& pixels,
                  const char* ordering, int nsideValue = kNside)
  {
    std::remove(file.c_str());

    fitsfile* fptr = nullptr;
    int status = 0;
    const std::string create = "!" + file;

    fits_create_file(&fptr, create.c_str(), &status);

    const long nrows = (long)pixels.size() / 1024;
    char name[] = "SIGNAL";
    char form[] = "1024E";
    char* ttype[] = {name};
    char* tform[] = {form};

    fits_create_tbl(fptr, BINARY_TBL, nrows, 1, ttype, tform, nullptr, nullptr, &status);

    int nside = nsideValue;
    fits_write_key(fptr, TINT, "NSIDE", &nside, nullptr, &status);
    fits_write_key(fptr, TSTRING, "ORDERING", (void*)ordering, nullptr, &status);

    fits_write_col(fptr, TFLOAT, 1, 1, 1, (long)pixels.size(),
                   (void*)pixels.data(), &status);

    fits_close_file(fptr, &status);
    return status;
  }

  // A result holding one displacement per (Eulerian, Lagrangian) pair of
  // directions, in one realization, all valid.
  Result result_of (const std::vector<double>& eulerian, const std::vector<double>& lagrangian)
  {
    Result r;
    r.nObjects = eulerian.size() / 3;
    r.nRealizations = 1;
    r.matchedRandom = lagrangian;
    r.displacement.resize(lagrangian.size());
    for (std::size_t k = 0; k < lagrangian.size(); ++k)
      r.displacement[k] = lagrangian[k] - eulerian[k];
    r.valid.assign(r.nObjects, 1);
    internal::summarize(r);
    return r;
  }

  // A displacement from (ra, dec) at distance 1000 to (ra2, dec2) at 1000.
  void add_arc (std::vector<double>& eulerian, std::vector<double>& lagrangian,
                const double ra, const double dec, const double ra2, const double dec2)
  {
    double x, y, z;
    internal::to_cartesian(ra, dec, 1000., x, y, z);
    eulerian.insert(eulerian.end(), {x, y, z});
    internal::to_cartesian(ra2, dec2, 1000., x, y, z);
    lagrangian.insert(lagrangian.end(), {x, y, z});
  }

  vec3 unit (const double x, const double y, const double z)
  {
    const double r = std::sqrt(x*x + y*y + z*z);
    return vec3(x/r, y/r, z/r);
  }

  // The same vector as otswap's pixel code takes it.
  internal::Vec3 vec_of (const vec3& v)
  {
    return internal::Vec3(v.x, v.y, v.z);
  }

  // Every unobserved pixel met by points spaced at most step pixels apart
  // along the arc, the endpoint pixels excluded.
  std::set<std::int64_t> brute_force (const internal::PixelMask& map, const vec3& a,
                                      const vec3& b, const double step)
  {
    const double pixelSize = std::sqrt(4. * kPi / (double)map.base.npix());
    const double angle = std::acos(std::max(-1., std::min(1., dotprod(a, b))));
    const long n = std::max(1L, (long)std::ceil(angle / (step * pixelSize)));
    const std::int64_t pixA = internal::vec2pix(map.base, vec_of(a)),
                       pixB = internal::vec2pix(map.base, vec_of(b));
    const double s = std::sin(angle);

    std::set<std::int64_t> found;
    for (long i = 0; i <= n; ++i) {
      const double t = (double)i / (double)n;
      const vec3 w = (angle < 1.e-12) ? a
        : (a * (std::sin((1. - t) * angle) / s) + b * (std::sin(t * angle) / s)).Norm();
      const std::int64_t p = internal::vec2pix(map.base, vec_of(w));
      if (p != pixA && p != pixB && !internal::pixel_observed(map, p)) found.insert(p);
    }
    return found;
  }

  // True when points spaced step pixels apart along the arc, over the part
  // of it within three pixels of the pixel's centre, meet the pixel.
  bool arc_enters (const internal::PixelMask& map, const vec3& a, const vec3& b,
                   const std::int64_t pixel, const double step)
  {
    const double pixelSize = std::sqrt(4. * kPi / (double)map.base.npix());
    const double angle = std::acos(std::max(-1., std::min(1., dotprod(a, b))));
    const double s = std::sin(angle);
    const T_Healpix_Base<std::int64_t> healpix(map.base.nside(), map.base.nested() ? NEST : RING,
                                               SET_NSIDE);
    const vec3 centre = healpix.pix2vec(pixel);
    const double near = std::cos(3. * pixelSize);
    auto at = [&] (const double t) {
      return (a * (std::sin((1. - t) * angle) / s) + b * (std::sin(t * angle) / s)).Norm();
    };

    const long coarse = std::max(1L, (long)std::ceil(angle / (1.e-2 * pixelSize)));
    long first = -1, last = -1;
    for (long i = 0; i <= coarse; ++i)
      if (dotprod(at((double)i / (double)coarse), centre) >= near) {
        if (first < 0) first = i;
        last = i;
      }
    if (first < 0) return false;

    const double t0 = (double)std::max(0L, first - 1) / (double)coarse;
    const double t1 = (double)std::min(coarse, last + 1) / (double)coarse;
    const long fine = std::max(1L, (long)std::ceil((t1 - t0) * angle / (step * pixelSize)));
    for (long i = 0; i <= fine; ++i)
      if (internal::vec2pix(map.base, vec_of(at(t0 + (t1 - t0) * (double)i / (double)fine))) == pixel)
        return true;
    return false;
  }

  // A random unit vector tangent to the sphere at u.
  vec3 tangent (const vec3& u, std::mt19937& rng)
  {
    for (;;) {
      const vec3 g(internal::uniform_real(rng, -1., 1.), internal::uniform_real(rng, -1., 1.),
                   internal::uniform_real(rng, -1., 1.));
      const vec3 t = g - u * dotprod(g, u);
      if (t.SquaredLength() > 1.e-6) return t.Norm();
    }
  }

  vec3 random_direction (std::mt19937& rng)
  {
    for (;;) {
      const vec3 g(internal::uniform_real(rng, -1., 1.), internal::uniform_real(rng, -1., 1.),
                   internal::uniform_real(rng, -1., 1.));
      const double r2 = g.SquaredLength();
      if (r2 > 1.e-6 && r2 <= 1.) return g.Norm();
    }
  }

  // Pixels unobserved at random, never two that are neighbours.
  std::vector<float> isolated_pixels (const Healpix_Map<float>& geometry, std::mt19937& rng)
  {
    const int npix = geometry.Npix();
    std::vector<float> pixels((std::size_t)npix, 1.f);
    std::vector<int> order((std::size_t)npix);
    for (int p = 0; p < npix; ++p) order[(std::size_t)p] = p;
    internal::shuffle(order, rng);

    fix_arr<int, 8> nb;
    for (const int p : order) {
      if (internal::uniform_real(rng, 0., 1.) > 0.15) continue;
      geometry.neighbors(p, nb);
      bool alone = true;
      for (int k = 0; k < 8; ++k)
        if (nb[k] >= 0 && pixels[(std::size_t)nb[k]] == 0.f) alone = false;
      if (alone) pixels[(std::size_t)p] = 0.f;
    }
    return pixels;
  }

  // Unobserved discs of 0.5 to 5 pixels in radius around random centres,
  // returned in centres.
  std::vector<float> holes (const Healpix_Map<float>& geometry, std::mt19937& rng,
                            std::vector<vec3>& centres)
  {
    const double pixelSize = std::sqrt(4. * kPi / (double)geometry.Npix());
    std::vector<float> pixels((std::size_t)geometry.Npix(), 1.f);
    centres.clear();
    for (int h = 0; h < 300; ++h) {
      const vec3 c = random_direction(rng);
      centres.push_back(c);
      const double radius = internal::uniform_real(rng, 0.5, 5.) * pixelSize;
      rangeset<int> disc;
      geometry.query_disc(pointing(c), radius, disc);
      for (std::size_t r = 0; r < disc.nranges(); ++r)
        for (int p = disc.ivbegin(r); p < disc.ivend(r); ++p) pixels[(std::size_t)p] = 0.f;
    }
    return pixels;
  }

  // Arcs of 0.05 to 12 pixels, one in ten of 12 to 60, starting anywhere
  // or, when centres are given, within 8 pixels of one of them. Written as
  // Eulerian and Lagrangian positions at independent distances.
  void random_arcs (const std::size_t count, const double pixelSize,
                    const std::vector<vec3>& centres, std::mt19937& rng,
                    std::vector<double>& eulerian, std::vector<double>& lagrangian)
  {
    eulerian.clear();
    lagrangian.clear();
    for (std::size_t i = 0; i < count; ++i) {
      vec3 a = random_direction(rng);
      if (!centres.empty()) {
        const vec3& c = centres[internal::uniform_int(rng, (unsigned)centres.size())];
        const double off = internal::uniform_real(rng, 0., 8.) * pixelSize;
        a = (c * std::cos(off) + tangent(c, rng) * std::sin(off)).Norm();
      }
      const double length = (internal::uniform_int(rng, 10) == 0)
        ? internal::uniform_real(rng, 12., 60.) : internal::uniform_real(rng, 0.05, 12.);
      const double angle = length * pixelSize;
      const vec3 b = (a * std::cos(angle) + tangent(a, rng) * std::sin(angle)).Norm();

      const double r0 = internal::uniform_real(rng, 500., 1500.);
      const double r1 = internal::uniform_real(rng, 500., 1500.);
      eulerian.insert(eulerian.end(), {r0 * a.x, r0 * a.y, r0 * a.z});
      lagrangian.insert(lagrangian.end(), {r1 * b.x, r1 * b.y, r1 * b.z});
    }
  }

  // The adaptive search against CosmoBolognaLib's arc_crosses_mask and
  // against a brute-force sampler, on one mask.
  void compare_on (const std::string& label, const std::vector<float>& pixels, const int nside,
                   const Healpix_Ordering_Scheme scheme, const std::vector<vec3>& centres,
                   const unsigned seed)
  {
    Healpix_Map<float> map(nside, scheme, SET_NSIDE);
    for (int p = 0; p < map.Npix(); ++p) map[p] = pixels[(std::size_t)p];

    // The same mask as otswap holds it, for its internal search.
    internal::PixelMask bytes;
    internal::set_geometry(bytes, nside, scheme == NEST, label);
    const std::vector<double> asDouble(pixels.begin(), pixels.end());
    internal::mark_observed(bytes, 0, asDouble.data(), asDouble.size());

    const std::string file = temporary("arcs_" + std::to_string(seed) + ".fits");
    check(write_map(file, pixels, scheme == NEST ? "NESTED" : "RING", nside) == 0,
          label + ": the mask is written");
    const Mask mask(file);
    std::remove(file.c_str());

    const double pixelSize = std::sqrt(4. * kPi / (double)map.Npix());

    std::mt19937 rng(seed);
    std::vector<double> eulerian, lagrangian;
    random_arcs(600, pixelSize, centres, rng, eulerian, lagrangian);
    const Result arcs = result_of(eulerian, lagrangian);
    const std::size_t n = arcs.nObjects;

    // The directions exactly as rejectMaskCrossings forms them.
    std::vector<vec3> from(n), to(n);
    std::vector<double> e(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
      for (int k = 0; k < 3; ++k)
        e[3*i+k] = arcs.matchedRandom[3*i+k] - arcs.displacement[3*i+k];
      from[i] = unit(e[3*i], e[3*i+1], e[3*i+2]);
      to[i] = unit(arcs.matchedRandom[3*i], arcs.matchedRandom[3*i+1], arcs.matchedRandom[3*i+2]);
    }

    std::size_t crossing = 0, cblDisagree = 0, filterDisagree = 0, misses = 0;
    std::size_t confirmed = 0, unconfirmed = 0;
    std::size_t countDisagree = 0, multiple = 0;
    const unsigned limits[5] = {1, 2, 3, 5, 8};

    Result filtered = arcs;
    rejectMaskCrossings(filtered, mask, 0);

    // What the arc crosses: every pixel the sampler at 0.001 pixel meets,
    // and every other one the search reports that a scan at 1e-6 pixel
    // confirms the arc enters.
    std::vector<std::set<std::int64_t>> brute(n), crossed(n);
    std::vector<std::int64_t> found;

    for (std::size_t i = 0; i < n; ++i) {
      const bool cbl = cbl_reference::arc_crosses_mask(
        map, e[3*i], e[3*i+1], e[3*i+2],
        arcs.matchedRandom[3*i], arcs.matchedRandom[3*i+1], arcs.matchedRandom[3*i+2], pixelSize);
      const bool ours = internal::arc_unobserved_pixels(bytes, vec_of(from[i]), vec_of(to[i]), 0, found) > 0;
      if (cbl) ++crossing;
      if (cbl != ours) ++cblDisagree;
      if (filtered.valid[i] != (cbl ? 0 : 1)) ++filterDisagree;

      internal::arc_unobserved_pixels(bytes, vec_of(from[i]), vec_of(to[i]), UINT_MAX, found);
      const std::set<std::int64_t> all(found.begin(), found.end());
      brute[i] = brute_force(bytes, from[i], to[i], 1.e-3);
      crossed[i] = brute[i];
      for (const std::int64_t p : brute[i]) if (all.count(p) == 0) ++misses;
      for (const std::int64_t p : all)
        if (brute[i].count(p) == 0) {
          if (arc_enters(bytes, from[i], to[i], p, 1.e-6)) { ++confirmed; crossed[i].insert(p); }
          else ++unconfirmed;
        }
      if (crossed[i].size() > 1) ++multiple;

      for (const unsigned k : limits) {
        const bool exceeds = internal::arc_unobserved_pixels(bytes, vec_of(from[i]), vec_of(to[i]), k, found) > k;
        if (exceeds != (crossed[i].size() > k)) ++countDisagree;
      }
    }

    for (const unsigned k : limits) {
      Result r = arcs;
      rejectMaskCrossings(r, mask, k);
      for (std::size_t i = 0; i < n; ++i)
        if (r.valid[i] != (crossed[i].size() > k ? 0 : 1)) ++countDisagree;
    }

    std::cout << "    " << label << ": " << crossing << " of " << n << " arcs cross, "
              << multiple << " cross more than one pixel; " << confirmed
              << " pixels found that the sampler at 0.001 pixel steps over" << std::endl;

    check(crossing > n / 10 && crossing < n - n / 10,
          label + ": the arcs include many that cross and many that do not");
    check(cblDisagree == 0,
          label + ": with a limit of 0 the decision is CosmoBolognaLib's, arc by arc (" +
          std::to_string(cblDisagree) + " differ)");
    check(filterDisagree == 0,
          label + ": and rejectMaskCrossings with 0 clears exactly those (" +
          std::to_string(filterDisagree) + " differ)");
    check(misses == 0,
          label + ": no pixel the brute-force sampler meets is missed (" +
          std::to_string(misses) + " missed)");
    check(unconfirmed == 0,
          label + ": every other pixel reported is one a scan at 1e-6 pixel confirms (" +
          std::to_string(unconfirmed) + " not confirmed)");
    check(multiple > n / 20, label + ": many arcs cross more than one pixel");
    check(countDisagree == 0,
          label + ": for limits above 0 the decision follows the count of pixels crossed (" +
          std::to_string(countDisagree) + " differ)");
  }


  // Writes count values as a Healpix map in the layout the reader assumes,
  // in a column of type form ("1024D", "1024E", "1024J", "1024B") written
  // from datatype. With tnull, the column's TNULL is set to it.
  long write_typed_map (const std::string& file, const char* form, const int datatype,
                        const void* values, const long count, const char* ordering,
                        const long long nsideValue, const long long* tnull = nullptr)
  {
    std::remove(file.c_str());

    fitsfile* fptr = nullptr;
    int status = 0;
    const std::string create = "!" + file;
    fits_create_file(&fptr, create.c_str(), &status);

    char name[] = "SIGNAL";
    std::string formText = form;
    char* ttype[] = {name};
    char* tform[] = {&formText[0]};
    fits_create_tbl(fptr, BINARY_TBL, count / 1024, 1, ttype, tform, nullptr, nullptr, &status);

    long long nside = nsideValue;
    fits_write_key(fptr, TLONGLONG, "NSIDE", &nside, nullptr, &status);
    fits_write_key(fptr, TSTRING, "ORDERING", (void*)ordering, nullptr, &status);
    if (tnull != nullptr) {
      long long t = *tnull;
      fits_write_key(fptr, TLONGLONG, "TNULL1", &t, nullptr, &status);
    }
    if (count > 0)
      fits_write_col(fptr, datatype, 1, 1, 1, count, const_cast<void*>(values), &status);

    fits_close_file(fptr, &status);
    return status;
  }

  // Pixel values of every kind: 0, negative, NaN, 1e-300, fractional and
  // 1, about one in seven unobserved.
  std::vector<double> mixed_values (const std::size_t npix, const unsigned seed)
  {
    std::mt19937 rng(seed);
    std::vector<double> v(npix);
    for (double& x : v)
      switch (internal::uniform_int(rng, 20)) {
        case 0: x = 0.; break;
        case 1: x = -internal::uniform_real(rng, 0., 2.); break;
        case 2: x = std::numeric_limits<double>::quiet_NaN(); break;
        case 3: x = 1.e-300; break;
        case 4: case 5: case 6: case 7: case 8:
          x = internal::uniform_real(rng, 0.25, 1.); break;
        default: x = 1.;
      }
    return v;
  }

  // Ways two masks differ: NSIDE, sky area, and allows on a 721 x 361 grid
  // of directions, the poles included.
  std::size_t disagreements (const Mask& a, const Mask& b)
  {
    std::size_t differ = (a.nside() != b.nside() ? 1 : 0) +
                         (a.skyAreaDeg2() != b.skyAreaDeg2() ? 1 : 0);
    for (int i = 0; i < 721; ++i)
      for (int j = 0; j <= 360; ++j) {
        const double ra = 2. * kPi * (double)i / 720.;
        const double dec = (kPi/2.) * (double)(j - 180) / 180.;
        if (a.allows(ra, dec) != b.allows(ra, dec)) ++differ;
      }
    return differ;
  }

  // Valid entries that differ when the same random arcs are filtered
  // against each mask, at limits 0 and 2.
  std::size_t filter_disagreements (const Mask& a, const Mask& b, const unsigned seed)
  {
    const double pixelSize = std::sqrt(kPi / 3.) / (double)a.nside();
    std::mt19937 rng(seed);
    std::vector<double> eulerian, lagrangian;
    random_arcs(600, pixelSize, std::vector<vec3>(), rng, eulerian, lagrangian);
    std::size_t differ = 0;
    for (const unsigned limit : {0u, 2u}) {
      Result ra = result_of(eulerian, lagrangian), rb = ra;
      rejectMaskCrossings(ra, a, limit);
      rejectMaskCrossings(rb, b, limit);
      for (std::size_t k = 0; k < ra.valid.size(); ++k)
        if (ra.valid[k] != rb.valid[k]) ++differ;
    }
    return differ;
  }

}

int main ()
{
  const T_Healpix_Base<int> base(kNside, RING, SET_NSIDE);
  const long npix = base.Npix();

  // A band in declination is forbidden; everything else is allowed.
  std::vector<float> pixels((std::size_t)npix, 1.f);
  long forbidden = 0;
  for (long p = 0; p < npix; ++p) {
    const pointing dir = base.pix2ang(p);
    const double dec = kPi/2. - dir.theta;
    if (std::fabs(dec) < 0.2) { pixels[(std::size_t)p] = 0.f; ++forbidden; }
  }

  const std::string file = temporary("mask.fits");
  check(write_map(file, pixels, "RING") == 0, "the synthetic mask is written");

  group("the mask is read, and reports its geometry");
  {
    const Mask mask(file);
    check(mask.nside() == kNside, "NSIDE comes from the header");

    const double fullSky = 4. * kPi * (180./kPi) * (180./kPi);
    const double expected = fullSky * (double)(npix - forbidden) / (double)npix;
    check_close(mask.skyAreaDeg2(), expected, 1.e-6,
                "the sky area counts the allowed pixels");

    check(mask.allows(1.0, 1.2), "a direction outside the band is allowed");
    check(!mask.allows(1.0, 0.0), "a direction inside the band is not");
    check(mask.allows(1.0 + 2.*kPi, 1.2),
          "right ascension is folded, so an unwrapped angle gives the same answer");
  }

  group("a file that is not a Healpix map in the expected layout is refused");
  {
    check_throws([] { const Mask m("/nonexistent/otswap/missing.fits"); (void)m; }, "a missing file raises");

    // Two rows of 1024 is 2048 pixels, which is not 12 * NSIDE^2 for any
    // NSIDE, so the layout check must fire.
    const std::string bad = temporary("badmask.fits");
    std::vector<float> wrong(2048, 1.f);
    {
      std::remove(bad.c_str());
      fitsfile* fptr = nullptr;
      int status = 0;
      const std::string create = "!" + bad;
      fits_create_file(&fptr, create.c_str(), &status);
      char name[] = "SIGNAL";
      char form[] = "1024E";
      char* ttype[] = {name};
      char* tform[] = {form};
      fits_create_tbl(fptr, BINARY_TBL, 2, 1, ttype, tform, nullptr, nullptr, &status);
      fits_write_col(fptr, TFLOAT, 1, 1, 1, 2048, (void*)wrong.data(), &status);
      fits_close_file(fptr, &status);
    }
    check_throws([&] { const Mask m(bad); (void)m; }, "a pixel count that is not 12 * NSIDE^2 raises");
    std::remove(bad.c_str());
  }

  group("any pixel value is accepted, and a pixel is observed when its value exceeds 0");
  {
    // Pixel 17 onwards hold the values under test; every other pixel is 1.
    const float values[] = {0.5f, 1.e-30f, 2.f, 1.f, 0.f, -0.f, -1.f, -1.6375e30f,
                            std::numeric_limits<float>::quiet_NaN(),
                            -std::numeric_limits<float>::infinity(),
                            std::numeric_limits<float>::infinity()};
    const bool observed[] = {true, true, true, true, false, false, false, false,
                             false, false, true};
    const std::size_t nValues = sizeof(values) / sizeof(values[0]);

    std::vector<float> mixed((std::size_t)npix, 1.f);
    for (std::size_t k = 0; k < nValues; ++k) mixed[17 + k] = values[k];

    const std::string file2 = temporary("fractional.fits");
    check(write_map(file2, mixed, "RING") == 0, "the mask of mixed values is written");

    long nObserved = 0;
    for (long p = 0; p < npix; ++p)
      if (mixed[(std::size_t)p] > 0.f) ++nObserved;

    try {
      const Mask m(file2);
      for (std::size_t k = 0; k < nValues; ++k) {
        const pointing centre = base.pix2ang(17 + (long)k);
        check(m.allows(centre.phi, kPi/2. - centre.theta) == observed[k],
              "the pixel holding " + std::to_string(values[k]) + " is " +
              (observed[k] ? "observed" : "unobserved"));
      }
      const double fullSky = 4. * kPi * (180./kPi) * (180./kPi);
      check(m.skyAreaDeg2() == fullSky * (double)nObserved / (double)npix,
            "the sky area counts the pixels above 0, whatever their values");
    }
    catch (const Error& e) {
      check(false, std::string("a mask of mixed values is read without error (") + e.what() + ")");
    }

    std::remove(file2.c_str());
  }

  group("the > 0 rule holds for any value, the same in memory and in FITS columns of type "
        "D, E, J with TNULL, and B");
  {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double fullSky = 4. * kPi * (180./kPi) * (180./kPi);

    // Pixel 17 onwards hold the values under test; every other pixel is 1.
    const double values[] = {0.5, 1.e-30, 1.e-50, 3.4e38, 1.e300, inf,
                             0., -0., -1., -1.e300, -inf, -1.6375e30, nan};
    const bool observed[] = {true, true, true, true, true, true,
                             false, false, false, false, false, false, false};
    const std::size_t nValues = sizeof(values) / sizeof(values[0]);

    std::vector<double> all(12 * kNside * kNside, 1.);
    for (std::size_t k = 0; k < nValues; ++k) all[17 + k] = values[k];
    long nObserved = 0;
    for (const double v : all) if (v > 0.) ++nObserved;

    const Mask memory(all, PixelOrdering::Ring);
    for (std::size_t k = 0; k < nValues; ++k) {
      const pointing centre = base.pix2ang(17 + (long)k);
      check(memory.allows(centre.phi, kPi/2. - centre.theta) == observed[k],
            "in memory, the pixel holding " + std::to_string(values[k]) + " is " +
            (observed[k] ? "observed" : "unobserved"));
    }
    check(memory.skyAreaDeg2() == fullSky * (double)nObserved / (double)npix,
          "in memory, the sky area counts the pixels above 0");

    // D: every value, +-inf and +-1e300 included, as double.
    const std::string fileD = temporary("rule_D.fits");
    check(write_typed_map(fileD, "1024D", TDOUBLE, all.data(), npix, "RING", kNside) == 0,
          "the D map is written");
    try {
      const Mask fromD(fileD);
      check(disagreements(fromD, memory) == 0,
            "a D column, +-inf and +-1e300 included, reads and gives the same mask as in memory");
    }
    catch (const Error& e) {
      check(false, std::string("a D column holding +-inf and +-1e300 reads (") + e.what() + ")");
    }
    std::remove(fileD.c_str());

    // E: the values a float holds.
    const float finf = std::numeric_limits<float>::infinity();
    const float fvalues[] = {0.5f, 1.e-30f, 3.4e38f, finf, 0.f, -0.f, -1.f, -finf, -1.6375e30f,
                             std::numeric_limits<float>::quiet_NaN()};
    std::vector<float> asFloat((std::size_t)npix, 1.f);
    std::vector<double> floatTwin((std::size_t)npix, 1.);
    for (std::size_t k = 0; k < sizeof(fvalues) / sizeof(fvalues[0]); ++k) {
      asFloat[17 + k] = fvalues[k];
      floatTwin[17 + k] = (double)fvalues[k];
    }
    const std::string fileE = temporary("rule_E.fits");
    check(write_typed_map(fileE, "1024E", TFLOAT, asFloat.data(), npix, "RING", kNside) == 0,
          "the E map is written");
    check(disagreements(Mask(fileE), Mask(floatTwin, PixelOrdering::Ring)) == 0,
          "an E column gives the same mask as its values in memory");
    std::remove(fileE.c_str());

    // J, with -999 as TNULL: a null pixel is unobserved.
    const long long tnull = -999;
    const int ivalues[] = {1, 0, -1, 7, -999, 2147483647, -2147483647 - 1};
    std::vector<int> asInt((std::size_t)npix, 1);
    std::vector<double> intTwin((std::size_t)npix, 1.);
    for (std::size_t k = 0; k < sizeof(ivalues) / sizeof(ivalues[0]); ++k) {
      asInt[17 + k] = ivalues[k];
      intTwin[17 + k] = ivalues[k] == -999 ? nan : (double)ivalues[k];
    }
    const std::string fileJ = temporary("rule_J.fits");
    check(write_typed_map(fileJ, "1024J", TINT, asInt.data(), npix, "RING", kNside, &tnull) == 0,
          "the J map is written");
    const Mask fromJ(fileJ);
    check(disagreements(fromJ, Mask(intTwin, PixelOrdering::Ring)) == 0,
          "a J column gives the same mask as its values in memory");
    const pointing nullPixel = base.pix2ang(17 + 4);
    check(!fromJ.allows(nullPixel.phi, kPi/2. - nullPixel.theta),
          "a J null pixel is unobserved");
    std::remove(fileJ.c_str());

    // B: bytes 0 and 1.
    std::vector<unsigned char> asByte((std::size_t)npix, 1);
    std::vector<double> byteTwin((std::size_t)npix, 1.);
    for (long p = 0; p < npix; p += 5) { asByte[(std::size_t)p] = 0; byteTwin[(std::size_t)p] = 0.; }
    const std::string fileB = temporary("rule_B.fits");
    check(write_typed_map(fileB, "1024B", TBYTE, asByte.data(), npix, "RING", kNside) == 0,
          "the B map is written");
    check(disagreements(Mask(fileB), Mask(byteTwin, PixelOrdering::Ring)) == 0,
          "a B column gives the same mask as its values in memory");
    std::remove(fileB.c_str());
  }

  group("the in-memory and FITS forms give the same mask, RING and NESTED, at NSIDE 16, 64 "
        "and 256, rejectMaskCrossings included; the ordering is honoured");
  {
    unsigned seed = 50;
    for (const int nside : {16, 64, 256}) {
      const std::size_t n = 12 * (std::size_t)nside * (std::size_t)nside;
      const std::vector<double> values = mixed_values(n, seed);
      const std::string label = "NSIDE " + std::to_string(nside);

      const Mask ring(values, PixelOrdering::Ring), nested(values, PixelOrdering::Nested);
      for (const bool nest : {false, true}) {
        const std::string file = temporary("twin.fits");
        check(write_typed_map(file, "1024D", TDOUBLE, values.data(), (long)n,
                              nest ? "NESTED" : "RING", nside) == 0, label + ": the map is written");
        const Mask fits(file);
        std::remove(file.c_str());
        const Mask& memory = nest ? nested : ring;
        const std::string what = label + (nest ? " NESTED" : " RING");
        check(memory.nside() == nside, what + ": NSIDE follows from the length");
        check(disagreements(fits, memory) == 0,
              what + ": NSIDE, the sky area and allows are the same in both forms");
        check(filter_disagreements(fits, memory, seed + (nest ? 1 : 0)) == 0,
              what + ": rejectMaskCrossings decides the same for both forms");
      }
      check(disagreements(ring, nested) > 1000,
            label + ": the same values read as RING and as NESTED give different masks");
      ++seed;
    }

    const std::vector<double> twelve = {1., 0., 1., 1., 0.5, 1., 1., -1., 1., 1., 1., 1.};
    const Mask one(twelve, PixelOrdering::Nested);
    check(one.nside() == 1, "12 values make an NSIDE 1 map");
    check_close(one.skyAreaDeg2(), 4. * kPi * (180./kPi) * (180./kPi) * 10. / 12., 1.e-9,
                "and its sky area counts its 10 observed pixels");
  }

  group("NSIDE: 12 NSIDE^2 values, NSIDE at most 2^29 and a power of 2 for NESTED, "
        "checked before anything is read or allocated; the same in both forms");
  {
    for (const std::size_t n : {12u, 48u, 192u, 108u})
      for (const PixelOrdering o : {PixelOrdering::Ring, PixelOrdering::Nested}) {
        const std::string what = std::to_string(n) + " values in " +
                                 (o == PixelOrdering::Ring ? "RING" : "NESTED");
        if (n == 108 && o == PixelOrdering::Nested) {
          check_throws([&] { const Mask m(std::vector<double>(n, 1.), o); (void)m; },
                       what + " (NSIDE 3) is refused: NESTED needs a power of 2");
          continue;
        }
        const Mask m(std::vector<double>(n, 1.), o);
        check(m.nside() * m.nside() * 12 == (int)n, what + " are accepted, with their NSIDE");
      }

    for (const std::size_t n : {0u, 11u, 13u, 24u})
      check_throws([&] { const Mask m(std::vector<double>(n, 1.), PixelOrdering::Ring); (void)m; },
                   std::to_string(n) + " values are not 12 NSIDE^2 and are refused");

    check_throws([] { const Mask m(nullptr, 12, PixelOrdering::Ring); (void)m; },
                 "a null pointer is refused");
    check_throws([] {
                   const Mask m(std::vector<double>(12, 1.), static_cast<PixelOrdering>(7));
                   (void)m;
                 }, "an ordering that is neither enumerator is refused");

    // Lengths far beyond memory: the check fires before a value is read
    // (the pointer holds one) or a byte allocated.
    const double one = 1.;
    for (const std::uint64_t nside : {(std::uint64_t(1) << 29) + 1, std::uint64_t(1) << 30})
      for (const PixelOrdering o : {PixelOrdering::Ring, PixelOrdering::Nested})
        check_throws([&] {
                       const Mask m(&one, (std::size_t)(12 * nside * nside), o);
                       (void)m;
                     }, "NSIDE " + std::to_string(nside) + " is refused by length, above 2^29");

    // From FITS: NSIDE 48 fills 27 rows of 1024, so it can be written in
    // the layout; RING is read, NESTED refused with otswap's Error.
    const std::size_t n48 = 12 * 48 * 48;
    const std::vector<double> ones(n48, 1.);
    const std::string file = temporary("nside48.fits");
    for (const bool nest : {false, true}) {
      check(write_typed_map(file, "1024D", TDOUBLE, ones.data(), (long)n48,
                            nest ? "NESTED" : "RING", 48) == 0, "the NSIDE 48 map is written");
      if (nest)
        check_throws([&] { const Mask m(file); (void)m; },
                     "a NESTED FITS map of NSIDE 48 is refused, with otswap's Error");
      else
        check(Mask(file).nside() == 48, "a RING FITS map of NSIDE 48 is read");
    }

    // A header NSIDE above 2^29 is refused before the table is read.
    for (const long long nside : {(1LL << 29) + 1, 1LL << 30, -16LL}) {
      check(write_typed_map(file, "1024D", TDOUBLE, ones.data(), 3072, "RING", nside) == 0,
            "the map with a header NSIDE of " + std::to_string(nside) + " is written");
      check_throws([&] { const Mask m(file); (void)m; },
                   "a header NSIDE of " + std::to_string(nside) + " is refused");
    }
    std::remove(file.c_str());
  }

  group("a FITS map longer than one block of 64 rows, not a multiple of it, reads as "
        "the same values in memory");
  {
    // NSIDE 96: 110,592 pixels in 108 rows, one block of 64 and one of 44.
    const std::size_t n = 12 * 96 * 96;
    const std::vector<double> values = mixed_values(n, 77);
    const std::string file = temporary("blocks.fits");
    check(write_typed_map(file, "1024D", TDOUBLE, values.data(), (long)n, "RING", 96) == 0,
          "the NSIDE 96 map is written");
    const Mask fits(file), memory(values, PixelOrdering::Ring);
    std::remove(file.c_str());
    check(disagreements(fits, memory) == 0, "both forms give the same mask");

    long nObserved = 0;
    for (const double v : values) if (v > 0.) ++nObserved;
    check(fits.skyAreaDeg2() == 4. * kPi * (180./kPi) * (180./kPi) * (double)nObserved / (double)n,
          "and the sky area counts every observed pixel of both blocks");
  }

  group("the values are read once and not kept");
  {
    std::vector<double> values = mixed_values(12 * 16 * 16, 91);
    const std::vector<double> original = values;
    const Mask mask(values, PixelOrdering::Ring);
    std::fill(values.begin(), values.end(), 0.);
    check(disagreements(mask, Mask(original, PixelOrdering::Ring)) == 0,
          "changing the source values afterwards changes nothing");
  }

  group("allows refuses a direction that is not one, with otswap's Error");
  {
    const Mask mask(file);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    check_throws([&] { (void)mask.allows(nan, 0.); }, "a NaN right ascension raises");
    check_throws([&] { (void)mask.allows(inf, 0.); }, "an infinite right ascension raises");
    check_throws([&] { (void)mask.allows(1., nan); }, "a NaN declination raises");
    check_throws([&] { (void)mask.allows(1., std::nextafter(kPi/2., 2.)); },
                 "a declination above pi/2 raises");
    check_throws([&] { (void)mask.allows(1., -std::nextafter(kPi/2., 2.)); },
                 "a declination below -pi/2 raises");
    check(mask.allows(1., kPi/2.) && mask.allows(1., -kPi/2.),
          "the poles themselves are accepted, and outside the band");
  }

  group("a displacement crossing the masked band is rejected, one clear of it kept");
  {
    const Mask mask(file);

    // The first stays well above the band; the second runs from below it
    // to above it, so its great-circle arc must pass through forbidden
    // pixels.
    std::vector<double> e, l;
    add_arc(e, l, 1.0, 0.8, 1.0, 0.9);
    add_arc(e, l, 1.0, -0.6, 1.0, 0.6);
    Result r = result_of(e, l);

    rejectMaskCrossings(r, mask);

    check(r.valid.size() == 2, "valid keeps one entry per displacement");
    check(r.valid[0] == 1 && r.validRealizations[0] == 1,
          "the displacement clear of the band survives");
    check(r.valid[1] == 0 && r.validRealizations[1] == 0,
          "the one crossing the band is rejected by the default, a limit of 0");

    for (int k = 0; k < 3; ++k)
      check_close(r.meanDisplacement[k], r.displacement[k], 1.e-12,
                  "the survivor keeps its mean displacement");

    for (int k = 3; k < 6; ++k)
      check(std::isnan(r.meanDisplacement[k]),
            "an object left with no valid realization has a NaN mean");
  }

  group("the mean is over the valid realizations only");
  {
    const Mask mask(file);

    // One object, two realizations: the first clear of the band, the
    // second across it.
    std::vector<double> e, l;
    add_arc(e, l, 1.0, 0.8, 1.0, 0.9);
    add_arc(e, l, 1.0, -0.6, 1.0, 0.6);
    Result two = result_of(e, l);
    two.nObjects = 1;
    two.nRealizations = 2;
    two.valid.assign(2, 1);
    internal::summarize(two);

    for (int k = 0; k < 3; ++k)
      check(two.meanDisplacement[k] == (two.displacement[k] + two.displacement[3+k]) / 2.,
            "unfiltered, the mean is over both realizations");

    rejectMaskCrossings(two, mask);
    check(two.valid[0] == 1 && two.valid[1] == 0, "the crossing realization alone is rejected");
    check(two.validRealizations[0] == 1, "validRealizations counts the survivor");
    for (int k = 0; k < 3; ++k)
      check(two.meanDisplacement[k] == two.displacement[k], "and the mean is the survivor's");
  }

  group("a tolerance wide enough admits the crossing");
  {
    const Mask mask(file);

    std::vector<double> e, l;
    add_arc(e, l, 1.0, -0.6, 1.0, 0.6);
    Result r = result_of(e, l);

    rejectMaskCrossings(r, mask, 10000);
    check(r.valid[0] == 1 && r.validRealizations[0] == 1,
          "a tolerance above the pixels crossed keeps the displacement");
  }

  group("the count is of distinct pixels, and depends on NSIDE");
  {
    // The band of |dec| < 0.2, crossed along a meridian, at NSIDE 16 and
    // at NSIDE 32.
    const T_Healpix_Base<int> fine(2 * kNside, RING, SET_NSIDE);
    std::vector<float> finePixels((std::size_t)fine.Npix(), 1.f);
    for (int p = 0; p < fine.Npix(); ++p)
      if (std::fabs(kPi/2. - fine.pix2ang(p).theta) < 0.2) finePixels[(std::size_t)p] = 0.f;
    const std::string fineFile = temporary("mask32band.fits");
    check(write_map(fineFile, finePixels, "RING", 2 * kNside) == 0, "the NSIDE 32 band is written");

    const Mask coarse(file), fineMask(fineFile);
    std::remove(fineFile.c_str());

    // The smallest limit that keeps the displacement is its pixel count.
    auto count = [] (const Mask& mask) {
      std::vector<double> e, l;
      add_arc(e, l, 1.0, -0.6, 1.0, 0.6);
      for (unsigned k = 0; k < 1000; ++k) {
        Result r = result_of(e, l);
        rejectMaskCrossings(r, mask, k);
        if (r.valid[0]) return k;
      }
      return 1000u;
    };

    const unsigned atCoarse = count(coarse), atFine = count(fineMask);
    std::cout << "    the band crossing counts " << atCoarse << " pixels at NSIDE " << kNside
              << " and " << atFine << " at NSIDE " << 2 * kNside << std::endl;
    check(atCoarse >= 4, "the arc crosses several unobserved pixels");
    check(atFine >= (3 * atCoarse) / 2,
          "doubling NSIDE raises the count of the same arc by about a factor of two");

    std::vector<double> e, l;
    add_arc(e, l, 1.0, -0.6, 1.0, 0.6);
    Result r = result_of(e, l);
    rejectMaskCrossings(r, fineMask, atCoarse);
    check(r.valid[0] == 0,
          "so a limit that keeps the arc at one NSIDE rejects it at the other");
  }

  group("a malformed or already filtered result is refused");
  {
    const Mask mask(file);

    Result empty;
    check_throws([&] { rejectMaskCrossings(empty, mask); }, "an empty result raises");

    Result ragged;
    ragged.nObjects = 2;
    ragged.nRealizations = 1;
    ragged.displacement.resize(3);
    ragged.matchedRandom.resize(6);
    ragged.valid.assign(2, 1);
    check_throws([&] { rejectMaskCrossings(ragged, mask); },
                 "arrays of the wrong size raise");

    std::vector<double> e, l;
    add_arc(e, l, 1.0, 0.8, 1.0, 0.9);
    Result noValid = result_of(e, l);
    noValid.valid.clear();
    check_throws([&] { rejectMaskCrossings(noValid, mask); }, "a valid of the wrong size raises");

    Result notBinary = result_of(e, l);
    notBinary.valid[0] = 2;
    check_throws([&] { rejectMaskCrossings(notBinary, mask); },
                 "a valid entry other than 0 or 1 raises");

    Result filtered = result_of(e, l);
    filtered.filteredNside = kNside * 2;
    check_throws([&] { rejectMaskCrossings(filtered, mask); },
                 "a result filtered against a different NSIDE raises");
  }

  group("the filter only clears entries, records its mask, and is idempotent");
  {
    const Mask mask(file);

    std::vector<double> e, l;
    add_arc(e, l, 1.0, 0.8, 1.0, 0.9);
    add_arc(e, l, 1.0, -0.6, 1.0, 0.6);
    add_arc(e, l, 2.0, 0.8, 2.0, 0.9);
    Result r = result_of(e, l);

    // An entry cleared beforehand, by the caller or another filter, stays
    // cleared though its arc is clear of the band.
    r.valid[2] = 0;

    check(r.filteredNside == 0, "an unfiltered result records no mask");
    rejectMaskCrossings(r, mask);
    check(r.filteredNside == kNside, "the filter records the NSIDE it used");
    check(r.valid[0] == 1 && r.valid[1] == 0 && r.valid[2] == 0,
          "the crossing is cleared and the entry cleared beforehand stays cleared");
    check(r.validRealizations[2] == 0 && std::isnan(r.meanDisplacement[6]),
          "and the summary follows valid");

    const Result once = r;
    rejectMaskCrossings(r, mask);
    check(r.valid == once.valid, "filtering again with the same mask changes nothing");
    check(r.validRealizations == once.validRealizations, "nor the valid counts");
    check(r.meanDisplacement[0] == once.meanDisplacement[0] &&
          r.meanDisplacement[1] == once.meanDisplacement[1] &&
          r.meanDisplacement[2] == once.meanDisplacement[2],
          "and leaves the mean displacement where it was");

    rejectMaskCrossings(r, mask, 10000);
    check(r.valid == once.valid, "a looser limit afterwards restores nothing");

    // A mask at a different NSIDE must now be refused.
    const std::string other = temporary("mask32.fits");
    if (write_map(other, std::vector<float>(12288, 1.f), "RING", 32) == 0) {
      const Mask wrongSide(other);
      check_throws([&] { rejectMaskCrossings(r, wrongSide); },
                   "a mask of a different NSIDE is refused");
      std::remove(other.c_str());
    }
  }

  std::remove(file.c_str());

  group("with a limit of 0 the search decides as CosmoBolognaLib's arc_crosses_mask, "
        "misses nothing a sampler at 0.001 pixel finds, and counts as it does");
  {
    unsigned seed = 1;
    for (const int nside : {64, 256})
      for (const Healpix_Ordering_Scheme scheme : {RING, NEST}) {
        const Healpix_Map<float> geometry(nside, scheme, SET_NSIDE);
        const std::string name = "NSIDE " + std::to_string(nside) +
                                 (scheme == NEST ? " NESTED" : " RING");

        std::mt19937 rng(1000 + seed);
        const std::vector<vec3> none;
        compare_on(name + ", isolated pixels", isolated_pixels(geometry, rng), nside, scheme,
                   none, seed++);

        std::vector<vec3> centres;
        const std::vector<float> holed = holes(geometry, rng, centres);
        compare_on(name + ", holes", holed, nside, scheme, centres, seed++);
      }
  }

  return report("test_mask");
}
