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
 *  @file tests/test_rsd.cpp
 *
 *  @brief The redshift-space correction of otswap/RSD.h: each step against
 *  a direct computation, the whole chain in both geometries, the bias
 *  table reader, and the agreement of a Result with its tracers.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fitsio.h>
#include <meshsearch/MeshGrid.h>

#include "otswap/OT.h"
#include "internal.h"
#include "check.h"

using namespace otswap;

namespace {

  const double kNaN = std::numeric_limits<double>::quiet_NaN();
  const double kInf = std::numeric_limits<double>::infinity();

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "otswap_test_rsd_" + name;
  }

  bool same_bits (const std::vector<double>& a, const std::vector<double>& b)
  {
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  }

  bool close (const double a, const double b, const double relative)
  {
    return std::fabs(a - b) <= relative * std::max(1., std::max(std::fabs(a), std::fabs(b)));
  }

  std::vector<double> uniform_points (const std::size_t n, const double size, std::mt19937& rng)
  {
    std::vector<double> p(3 * n);
    for (double& v : p) v = internal::uniform_real(rng, 0., size);
    return p;
  }

  // The average written out over every pair, with the system exp.
  std::vector<double> direct_average (const std::vector<double>& p, const std::vector<double>& v,
                                      const std::vector<unsigned>& r, const double sigma,
                                      const bool weighted, std::vector<unsigned>& nn,
                                      std::vector<unsigned>& nr)
  {
    const std::size_t n = v.size();
    std::vector<double> out(n, kNaN);
    nn.assign(n, 0);
    nr.assign(n, 0);
    const double cutoff2 = (3. * sigma) * (3. * sigma);
    for (std::size_t i = 0; i < n; ++i) {
      double sum = 0., norm = 0.;
      for (std::size_t j = 0; j < n; ++j) {
        if (r[j] == 0) continue;
        const double dx = p[3*j] - p[3*i], dy = p[3*j+1] - p[3*i+1], dz = p[3*j+2] - p[3*i+2];
        const double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 > cutoff2) continue;
        const double w = (weighted ? (double)r[j] : 1.) * std::exp(-d2 / (2. * sigma * sigma));
        sum += w * v[j];
        norm += w;
        ++nn[i];
        nr[i] += r[j];
      }
      if (nn[i] > 0) out[i] = sum / norm;
    }
    return out;
  }

  // A Result for N tracers holding only what the correction reads: the
  // mean displacements and the valid realizations.
  Result result_with (const std::vector<double>& mean, const std::vector<unsigned>& valid)
  {
    Result r;
    r.nObjects = valid.size();
    r.nRealizations = 1;
    r.meanDisplacement = mean;
    r.validRealizations = valid;
    return r;
  }

  // Comoving distance 3000 z Mpc/h on [0, zMax], constant growth rate f.
  DistanceTable linear_table (const double zMax, const double f)
  {
    std::vector<double> z, d, g;
    for (int k = 0; k <= 1000; ++k) {
      z.push_back(zMax * k / 1000.);
      d.push_back(3000. * zMax * k / 1000.);
      g.push_back(f);
    }
    return DistanceTable(z, d, g);
  }

  std::string captured_clog (const std::function<void()>& f)
  {
    std::ostringstream text;
    std::streambuf* old = std::clog.rdbuf(text.rdbuf());
    try {
      f();
    }
    catch (...) {
      std::clog.rdbuf(old);
      throw;
    }
    std::clog.rdbuf(old);
    return text.str();
  }

  std::size_t lines (const std::string& s)
  {
    std::size_t n = 0;
    for (const char c : s) if (c == '\n') ++n;
    return n;
  }

  template <typename F>
  bool throws_naming (F&& f, const std::string& text)
  {
    try {
      f();
    }
    catch (const Error& e) {
      return std::string(e.what()).find(text) != std::string::npos;
    }
    return false;
  }

  void write_text (const std::string& file, const std::string& text)
  {
    std::ofstream out(file);
    out << text;
  }

  // A full-sky Healpix mask of NSIDE 16, every pixel observed, in the
  // layout Mask reads: one vector column of 1024 floats per row in HDU 2.
  void write_full_mask (const std::string& file)
  {
    std::remove(file.c_str());
    std::vector<float> pixels(12 * 16 * 16, 1.f);
    fitsfile* fptr = nullptr;
    int status = 0;
    const std::string create = "!" + file;
    fits_create_file(&fptr, create.c_str(), &status);
    char name[] = "SIGNAL";
    char form[] = "1024E";
    char* ttype[] = {name};
    char* tform[] = {form};
    fits_create_tbl(fptr, BINARY_TBL, 3, 1, ttype, tform, nullptr, nullptr, &status);
    int nside = 16;
    char ordering[] = "RING";
    fits_write_key(fptr, TINT, "NSIDE", &nside, nullptr, &status);
    fits_write_key(fptr, TSTRING, "ORDERING", ordering, nullptr, &status);
    fits_write_col(fptr, TFLOAT, 1, 1, 1, (long)pixels.size(), pixels.data(), &status);
    fits_close_file(fptr, &status);
  }

  // Sky coordinates of n objects in a patch, redshifts uniform in [z0, z1].
  std::vector<double> sky_points (const std::size_t n, const double z0, const double z1,
                                  std::mt19937& rng)
  {
    std::vector<double> sky(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
      sky[3*i]   = internal::uniform_real(rng, 0.2, 0.7);
      sky[3*i+1] = internal::uniform_real(rng, -0.25, 0.25);
      sky[3*i+2] = internal::uniform_real(rng, z0, z1);
    }
    return sky;
  }

  std::vector<std::size_t> inside (const std::vector<double>& sky, const double lo, const double hi)
  {
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < sky.size() / 3; ++i)
      if (sky[3*i+2] >= lo && sky[3*i+2] <= hi) keep.push_back(i);
    return keep;
  }

  std::vector<double> rows (const std::vector<double>& a, const std::vector<std::size_t>& keep)
  {
    std::vector<double> out;
    for (const std::size_t i : keep) out.insert(out.end(), a.begin() + 3*i, a.begin() + 3*i + 3);
    return out;
  }

  bool same_double (const double a, const double b)
  {
    return std::memcmp(&a, &b, sizeof a) == 0;
  }

}

int main ()
{
  std::mt19937 rng(20261006);

  group("projection: the radial and the axis component, NaN carried through, bad input refused");
  {
    const std::vector<double> p {3., 4., 0., 0., 0., -2.};
    const std::vector<double> d {1., 2., 3., kNaN, 0., 0.};
    const std::vector<double> radial = lineOfSightProjection(p, d);
    check(radial.size() == 2 && radial[0] == (3.*1. + 4.*2. + 0.*3.) / 5., "r.d/|r| for a position");
    check(std::isnan(radial[1]), "NaN for a NaN displacement");

    const std::vector<double> axis = lineOfSightProjection(d, 2);
    check(axis[0] == 3. && axis[1] == 0., "the axis component, copied");

    check_throws([&] { lineOfSightProjection({0., 0., 0.}, {1., 1., 1.}); }, "a position at the origin raises");
    check_throws([&] { lineOfSightProjection(p, {1., 2., 3.}); }, "a size mismatch raises");
    check_throws([&] { lineOfSightProjection(p, {1., 2., kInf, 0., 0., 0.}); }, "an infinite displacement raises");
    check_throws([&] { lineOfSightProjection(d, 3); }, "an axis beyond 2 raises");
    check(throws_naming([&] { lineOfSightProjection({1., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0.}); }, "object 1"),
          "the message names the object at the origin");
  }

  const std::size_t n = 3000;
  const std::vector<double> points = uniform_points(n, 100., rng);
  std::vector<double> values(n);
  std::vector<unsigned> valid(n);
  for (std::size_t i = 0; i < n; ++i) {
    valid[i] = internal::uniform_int(rng, 4);   // a quarter of the objects invalid
    values[i] = valid[i] > 0 ? internal::uniform_real(rng, -5., 5.) : kNaN;
  }

  group("average: equal to the sum over every pair, weighted or not, with the same diagnostics");
  {
    for (const bool weighted : {false, true}) {
      std::vector<unsigned> nn, nr, dn, dr;
      const std::vector<double> got = neighbourAverage(points, values, valid, 6., weighted, nn, nr);
      const std::vector<double> expected = direct_average(points, values, valid, 6., weighted, dn, dr);
      std::size_t bad = 0, nans = 0;
      for (std::size_t i = 0; i < n; ++i) {
        if (std::isnan(expected[i])) { if (!std::isnan(got[i])) ++bad; ++nans; continue; }
        if (!close(got[i], expected[i], 1.e-12)) ++bad;
      }
      const std::string mode = weighted ? " (weighted by realizations)" : "";
      check(bad == 0, std::to_string(bad) + " averages differ from the direct sum" + mode);
      check(nn == dn && nr == dr, "the diagnostics are the direct counts" + mode);
      check(nans < n / 10, "few objects are left without a neighbour" + mode);
    }
  }

  group("sigma = 0 returns the input unchanged, with diagnostics 1 and the own count for a valid object");
  {
    std::vector<double> anything = values;
    anything[0] = 123.;          // an invalid object's value is returned as given
    std::vector<unsigned> nn, nr;
    const std::vector<double> got = neighbourAverage(points, anything, valid, 0., false, nn, nr);
    check(same_bits(got, anything), "the values come back bit for bit");
    bool diagnostics = true;
    for (std::size_t i = 0; i < n; ++i)
      diagnostics = diagnostics && nn[i] == (valid[i] > 0 ? 1u : 0u) && nr[i] == valid[i];
    check(diagnostics, "and the diagnostics count the object itself");
  }

  group("an invalid object receives the average of its valid neighbours; NaN only without any");
  {
    // Three valid objects around the origin, an invalid one among them and
    // an invalid one 100 away; they span a volume, as the grid needs.
    const std::vector<double> p {0., 0., 0.,  1., 0., 0.3,  0., 1., -0.2,  0.5, 0.5, 0.1,  100., 0., 5.};
    const std::vector<double> v {1., 2., 3., kNaN, kNaN};
    const std::vector<unsigned> r {1, 1, 1, 0, 0};
    std::vector<unsigned> nn, nr, dn, dr;
    const std::vector<double> got = neighbourAverage(p, v, r, 2., false, nn, nr);
    const std::vector<double> expected = direct_average(p, v, r, 2., false, dn, dr);
    check(close(got[3], expected[3], 1.e-14) && nn[3] == 3, "the invalid object amid valid ones is averaged");
    check(std::isnan(got[4]) && nn[4] == 0 && nr[4] == 0, "the isolated one receives NaN");
    check(got[0] != 1. && nn[0] == 3, "a valid object is averaged with its neighbours too");
  }

  group("normalisation: a constant comes back unchanged everywhere, at the edges of the catalogue too");
  {
    const std::vector<double> constant(n, 0.7);
    std::vector<unsigned> all(n, 2);
    const std::vector<double> got = neighbourAverage(points, constant, all, 8., true);
    double worst = 0.;
    for (const double g : got) worst = std::max(worst, std::fabs(g - 0.7));
    check(worst < 1.e-14, "the largest deviation from the constant is " + std::to_string(worst));
  }

  group("the average does not depend on the internal cell within rounding, and its counts not at all");
  {
    const double cell = internal::neighbour_cell(points);
    check(cell > 0. && std::isfinite(cell), "the cell is positive");
    std::vector<unsigned> nn0, nr0;
    const std::vector<double> reference =
      internal::neighbour_average(points, values, valid, 6., false, cell, nn0, nr0);
    for (const double factor : {0.3, 0.5, 2., 5.}) {
      std::vector<unsigned> nn, nr;
      const std::vector<double> got =
        internal::neighbour_average(points, values, valid, 6., false, cell * factor, nn, nr);
      std::size_t bad = 0;
      for (std::size_t i = 0; i < n; ++i)
        if (!(std::isnan(got[i]) && std::isnan(reference[i])) && !close(got[i], reference[i], 1.e-13)) ++bad;
      check(bad == 0 && nn == nn0 && nr == nr0,
            "a cell " + std::to_string(factor) + " times as large changes nothing beyond rounding");
    }
  }

#ifdef _OPENMP
  group("the average is the same bit for bit with 1 and with 4 threads");
  {
    const int before = omp_get_max_threads();
    omp_set_num_threads(1);
    const std::vector<double> one = neighbourAverage(points, values, valid, 6., true);
    omp_set_num_threads(4);
    const std::vector<double> four = neighbourAverage(points, values, valid, 6., true);
    omp_set_num_threads(before);
    check(same_bits(one, four), "1 and 4 threads agree");
  }
#endif

  group("the internal cell is 4 (V/N)^(1/3); positions spanning no volume are refused by the grid");
  {
    const double cell = internal::neighbour_cell(points);
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) { lo[k] = hi[k] = points[(std::size_t)k]; }
    for (std::size_t i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], points[3*i+k]); hi[k] = std::max(hi[k], points[3*i+k]); }
    const double volume = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
    check(close(cell, 4. * std::cbrt(volume / (double)n), 1.e-15), "the cell of a catalogue filling its box");

    std::vector<double> plane = uniform_points(400, 50., rng), line = uniform_points(400, 50., rng);
    for (std::size_t i = 0; i < 400; ++i) { plane[3*i+2] = 7.; line[3*i+1] = 0.; line[3*i+2] = 0.; }
    const std::vector<double> single {1., 2., 3.};
    const std::vector<double> same {2., 2., 2.,  2., 2., 2.,  2., 2., 2.};
    for (const auto& c : {std::make_pair(std::string("a single object"), single),
                          std::make_pair(std::string("a plane"), plane),
                          std::make_pair(std::string("a line"), line),
                          std::make_pair(std::string("coincident objects"), same)}) {
      const std::size_t m = c.second.size() / 3;
      const std::vector<double> v(m, 1.);
      const std::vector<unsigned> r(m, 1);
      check(internal::neighbour_cell(c.second) == 0., c.first + ": the cell is 0");
      bool refused = false;
      try {
        neighbourAverage(c.second, v, r, 3.);
      }
      catch (const meshsearch::Error&) {
        refused = true;
      }
      check(refused, c.first + ": the grid refuses it with meshsearch's error, before any parallel region");
      check(neighbourAverage(c.second, v, r, 0.) == v, c.first + ": with sigma = 0 there is no grid, and no error");
    }
  }

  group("factor: b(z) interpolated between nodes, extrapolated from the end segments, f/(b + 3f/5)");
  {
    const DistanceTable table = linear_table(3., 0.8);
    const std::vector<double> zb {0.5, 1., 2.}, b {1., 1.5, 3.5};
    const std::vector<double> z {0.75, 1., 1.5, 2.5, 0.25};
    const double bz[] = {1.25, 1.5, 2.5, 4.5, 0.75};
    std::size_t extrapolated = 99;
    const std::vector<double> got = rsdFactor(z, table, zb, b, extrapolated);
    bool ok = true;
    for (std::size_t i = 0; i < z.size(); ++i)
      ok = ok && close(got[i], 0.8 / (bz[i] + 3. * 0.8 / 5.), 1.e-14);
    check(ok, "the factor at nodes, between them and beyond both ends");
    check(extrapolated == 2, "the two redshifts beyond the table are counted");

    const std::string once = captured_clog([&] { rsdFactor(z, table, zb, b); });
    check(lines(once) == 1 && once.find("2 of 5") != std::string::npos,
          "one line on std::clog per call that extrapolates, with the count");
    const std::string none = captured_clog([&] { rsdFactor({0.6, 1.9}, table, zb, b); });
    check(none.empty(), "and nothing when every redshift is inside the table");
    check(captured_clog([&] { rsdFactor(z, table, zb, b, extrapolated); }).empty(),
          "the counting overload writes nothing");

    check(throws_naming([&] { rsdFactor({1.5, 0.1}, table, {1., 2.}, {0.2, 1.2}); }, "object 1"),
          "an extrapolated bias that is not positive raises, naming the object");
    check_throws([&] { rsdFactor(z, table, {1.}, {1.}); }, "a single node raises");
    check_throws([&] { rsdFactor(z, table, {1., 2.}, {1.}); }, "unequal lengths raise");
    check_throws([&] { rsdFactor(z, table, {1., 1.}, {1., 2.}); }, "redshifts that do not increase raise");
    check_throws([&] { rsdFactor(z, table, {1., 2.}, {1., 0.}); }, "a bias of 0 raises");
    check_throws([&] { rsdFactor(z, table, {1., kNaN}, {1., 2.}); }, "a NaN node raises");
    check(throws_naming([&] { rsdFactor({1., 3.5}, table, zb, b); }, "object 1"),
          "a redshift outside the distance table raises, naming the object");
    check_throws([&] { rsdFactor({kNaN}, table, zb, b); }, "a NaN redshift raises");
    const DistanceTable noGrowth({0., 1.}, {0., 3000.});
    check_throws([&] { rsdFactor({0.5}, noGrowth, zb, b); }, "a table without growth rate raises");

    check(close(rsdFactorBox(1., table, 2.), 0.8 / (2. + 0.48), 1.e-15), "the box factor");
    check_throws([&] { rsdFactorBox(1., table, 0.); }, "a box bias of 0 raises");
    check_throws([&] { rsdFactorBox(kNaN, table, 1.); }, "a NaN box redshift raises");
    check_throws([&] { rsdFactorBox(1., noGrowth, 1.); }, "a box table without growth rate raises");
  }

  group("shift: radial moves the distance by the shift, the axis moves one coordinate; NaN carried");
  {
    const std::vector<double> p {3., 4., 12.,  1., 1., 1.};
    const std::vector<double> radial = shiftAlongLineOfSight(p, {2., kNaN});
    const double r = std::sqrt(radial[0]*radial[0] + radial[1]*radial[1] + radial[2]*radial[2]);
    check(close(r, 15., 1.e-15) && close(radial[0] / radial[2], 3. / 12., 1.e-15),
          "a positive shift moves outwards, along the line of sight");
    check(std::isnan(radial[3]) && std::isnan(radial[4]) && std::isnan(radial[5]), "a NaN shift gives a NaN row");

    const std::vector<double> axis = shiftAlongLineOfSight(p, {2., -1.}, 1);
    check(axis[0] == 3. && axis[1] == 6. && axis[2] == 12. && axis[4] == 0., "the axis coordinate alone moves");
    check_throws([&] { shiftAlongLineOfSight(p, {kInf, 0.}); }, "an infinite shift raises");
    check_throws([&] { shiftAlongLineOfSight(p, {1.}); }, "a size mismatch raises");
    check_throws([&] { shiftAlongLineOfSight({0., 0., 0.}, {1.}); }, "a position at the origin raises");
  }

  group("sign: a tracer whose displacement points away from the observer moves away");
  {
    const DistanceTable table = linear_table(1.5, 0.8);
    const std::vector<double> sky {0.3, 0.2, 0.5};
    const Result outward = result_with({0., 0., 0.}, {1});
    Result r = outward;
    double x, y, z;
    internal::to_cartesian(0.3, 0.2, 1., x, y, z);
    r.meanDisplacement = {5. * x, 5. * y, 5. * z};
    const RealSpaceCatalog c = realSpaceLightcone(r, sky, table, {0., 1.5}, {1.5, 1.5}, 0.);
    check(c.positions[2] > 0.5, "the lightcone redshift grows");
    check(close(c.positions[2], 0.5 + 5. * 0.8 / (1.5 + 0.48) / 3000., 1.e-12),
          "by the factor times the displacement, in distance");

    const Result box = result_with({0., 0., -2.}, {1});
    const RealSpaceCatalog cb = realSpaceBox(box, {10., 10., 10.}, 2, 0.5, table, 1.5, 0.);
    check(cb.positions[2] < 10. && cb.positions[0] == 10., "the box coordinate follows the displacement's sign");
  }

  group("box and lightcone agree on an equivalent setup: a small patch far along the x axis");
  {
    const DistanceTable table = linear_table(1.5, 0.8);
    const std::size_t m = 600;
    std::vector<double> cart(3 * m), sky(3 * m), mean(3 * m);
    for (std::size_t i = 0; i < m; ++i) {
      const double x = 2000. + internal::uniform_real(rng, 0., 40.);
      const double y = internal::uniform_real(rng, -20., 20.), z = internal::uniform_real(rng, -20., 20.);
      double ra, dec, d;
      internal::to_sky(x, y, z, ra, dec, d);
      sky[3*i] = ra; sky[3*i+1] = dec; sky[3*i+2] = d / 3000.;
      internal::to_cartesian(ra, dec, table.distanceAt(d / 3000.), cart[3*i], cart[3*i+1], cart[3*i+2]);
      for (int k = 0; k < 3; ++k) mean[3*i+k] = internal::uniform_real(rng, -2., 2.);
    }
    const Result r = result_with(mean, std::vector<unsigned>(m, 1));
    const RealSpaceCatalog light = realSpaceLightcone(r, sky, table, {0., 1.5}, {1.5, 1.5}, 8.);
    const RealSpaceCatalog box = realSpaceBox(r, cart, 0, 0.7, table, 1.5, 8.);
    double worst = 0.;
    for (std::size_t i = 0; i < m; ++i) {
      const double radial = table.distanceAt(light.positions[3*i+2]) - table.distanceAt(sky[3*i+2]);
      const double along = box.positions[3*i] - cart[3*i];
      worst = std::max(worst, std::fabs(radial - along));
    }
    check(worst < 0.03, "the radial and the x shifts differ by at most " + std::to_string(worst) +
                        " Mpc/h, for shifts of order 1");
    check(light.nNeighbours == box.nNeighbours, "with the same neighbours");
  }

  group("the whole chain: right ascension and declination copied, uncorrected only without valid neighbours");
  {
    const DistanceTable table = linear_table(1.5, 0.8);
    std::vector<double> sky {0.1, 0.1, 0.5,  0.1001, 0.1, 0.5,  0.1, 0.1001, 0.5,  1.0, -0.3, 0.6};
    const Result r = result_with({1., 1., 1.,  kNaN, kNaN, kNaN,  2., 0., 0.,  kNaN, kNaN, kNaN}, {2, 0, 1, 0});
    const RealSpaceCatalog c = realSpaceLightcone(r, sky, table, {0., 1.5}, {1.5, 1.5}, 10.);
    check(c.nObjects == 4 && c.positions[0] == 0.1 && c.positions[1] == 0.1 && c.positions[4] == 0.1,
          "the angles are the input ones");
    check(!std::isnan(c.positions[5]) && c.nNeighbours[1] == 2 && c.nRealizationsAveraged[1] == 3,
          "the invalid tracer amid valid ones is corrected with their average");
    check(c.uncorrected == std::vector<std::size_t>{3} && std::isnan(c.positions[9]) &&
          std::isnan(c.positions[11]), "the isolated invalid tracer alone is uncorrected, a NaN row");

    const RealSpaceCatalog zero = realSpaceLightcone(r, sky, table, {0., 1.5}, {1.5, 1.5}, 0.);
    check((zero.uncorrected == std::vector<std::size_t>{1, 3}), "with sigma = 0 every invalid tracer is uncorrected");

    check_throws([&] { realSpaceLightcone(r, {0.1, 0.1, 0.5}, table, {0., 1.5}, {1.5, 1.5}, 1.); },
                 "a result for another number of tracers raises");
    Result bad = r;
    bad.meanDisplacement.pop_back();
    check_throws([&] { realSpaceLightcone(bad, sky, table, {0., 1.5}, {1.5, 1.5}, 1.); },
                 "a malformed result raises");
    check_throws([&] { realSpaceBox(r, {1., 1., 1., 2., 2., 2., 3., 3., 3., 4., 4., 4.}, 3, 0.5, table, 1.5, 1.); },
                 "a box axis beyond 2 raises");
  }

  group("a corrected distance outside the table raises, naming the object and suggesting a wider table");
  {
    const DistanceTable table = linear_table(0.5, 0.8);
    double x, y, z;
    internal::to_cartesian(0.2, 0.1, 1., x, y, z);
    const Result r = result_with({0., 0., 0.,  400. * x, 400. * y, 400. * z}, {1, 1});
    const std::vector<double> sky {0.5, 0.1, 0.3,  0.2, 0.1, 0.49};
    bool named = false, suggests = false;
    try {
      realSpaceLightcone(r, sky, table, {0., 0.5}, {1.5, 1.5}, 0.);
    }
    catch (const Error& e) {
      named = std::string(e.what()).find("object 1") != std::string::npos;
      suggests = std::string(e.what()).find("wider redshift range") != std::string::npos;
    }
    check(named && suggests, "beyond the table's end");

    const Result inward = result_with({0., 0., 0.,  -5000. * x, -5000. * y, -5000. * z}, {1, 1});
    check(throws_naming([&] { realSpaceLightcone(inward, sky, table, {0., 0.5}, {1.5, 1.5}, 0.); }, "not positive"),
          "and below zero distance");
  }

  group("a Result agrees with its tracers: matchedRandom - displacement is the tracer, within rounding");
  {
    const std::vector<double> tracers = uniform_points(400, 80., rng);
    Config config;
    config.nRealizations = 3;
    config.seed = 7;
    const Result box = reconstructBox(tracers, {}, 80. / std::cbrt(400.), config);

    const DistanceTable table(0.3, 0.7, -1., 0., 0., 1.5, 4000);
    std::vector<double> sky(3 * 600), randomsSky(3 * 2400);
    for (std::size_t i = 0; i < sky.size(); i += 3) {
      sky[i] = internal::uniform_real(rng, 0.2, 0.7);
      sky[i+1] = internal::uniform_real(rng, -0.25, 0.25);
      sky[i+2] = internal::uniform_real(rng, 0.3, 0.6);
    }
    for (std::size_t i = 0; i < randomsSky.size(); i += 3) {
      randomsSky[i] = internal::uniform_real(rng, 0.2, 0.7);
      randomsSky[i+1] = internal::uniform_real(rng, -0.25, 0.25);
      randomsSky[i+2] = internal::uniform_real(rng, 0.3, 0.6);
    }
    const Result light = reconstructLightcone(sky, randomsSky, 800., 1, table, config);
    const std::vector<double> cart = toCartesian(sky, table);

    const std::pair<const Result*, const std::vector<double>*> cases[] = {{&box, &tracers}, {&light, &cart}};
    for (const auto& c : cases) {
      const Result& r = *c.first;
      const std::vector<double>& t = *c.second;
      double worst = 0.;
      for (unsigned rec = 0; rec < r.nRealizations; ++rec)
        for (std::size_t k = 0; k < t.size(); ++k) {
          const std::size_t at = (std::size_t)rec * t.size() + k;
          worst = std::max(worst, std::fabs(r.matchedRandom[at] - r.displacement[at] - t[k]));
        }
      check(worst <= 1.e-12 * 3000., "every realization, largest difference " + std::to_string(worst));
    }
  }

  group("readBiasTable: ASCII and FITS, and every error");
  {
    const std::string ascii = temporary("bias.txt");
    write_text(ascii, "# z b\n0.1 1.2\n\n0.5 1.6   # trailing comment ignored by strtod\n1.0 2.0\n");
    const io::BiasTable a = io::readBiasTable(ascii);
    check(a.redshift == std::vector<double>({0.1, 0.5, 1.0}) && a.bias == std::vector<double>({1.2, 1.6, 2.0}),
          "ASCII: columns 0 and 1, comments and blank lines skipped");

    const std::string fits = temporary("bias.fits");
    io::write(fits, {{"REDSHIFT", 'D', "", {0.2, 0.4}}, {"BIAS", 'D', "", {1.1, 1.3}}});
    const io::BiasTable f = io::readBiasTable(fits);
    check(f.redshift == std::vector<double>({0.2, 0.4}) && f.bias == std::vector<double>({1.1, 1.3}),
          "FITS: the REDSHIFT and BIAS columns");

    check_throws([] { io::readBiasTable("/nonexistent/otswap/bias.txt"); }, "a missing file raises");

    const std::string wrong = temporary("wrong.fits");
    io::write(wrong, {{"REDSHIFT", 'D', "", {0.2, 0.4}}, {"B", 'D', "", {1.1, 1.3}}});
    check_throws([&] { io::readBiasTable(wrong); }, "FITS without a BIAS column raises");

    const std::string nonFinite = temporary("nan.fits");
    io::write(nonFinite, {{"REDSHIFT", 'D', "", {0.2, kNaN}}, {"BIAS", 'D', "", {1.1, 1.3}}});
    check_throws([&] { io::readBiasTable(nonFinite); }, "FITS with a NaN raises");

    const std::pair<const char*, const char*> cases[] = {
      {"0.1\n0.5 1.6\n", "row too short"},
      {"0.1 1.2\n", "at least 2"},
      {"0.1 1.2\nnan 1.6\n", "data row 2"},
      {"0.5 1.2\n0.5 1.6\n", "data row 2"},
      {"0.5 1.2\n0.4 1.6\n", "not greater"},
      {"0.1 1.2\n0.5 0\n", "data row 2"},
      {"0.1 -1\n0.5 1\n", "data row 1"},
      {"0.1 1.2\n0.5 nan\n", "data row 2"},
      {"0.1 1.2\n0.5 inf\n", "data row 2"},
    };
    for (const auto& c : cases) {
      write_text(ascii, c.first);
      check(throws_naming([&] { io::readBiasTable(ascii); }, c.second),
            std::string("ASCII \"") + c.first + "\" raises with \"" + c.second + "\"");
    }

    for (const std::string& file : {ascii, fits, wrong, nonFinite}) std::remove(file.c_str());
  }

  group("redshift cut: the kept tracers are reconstructed exactly as the pre-cut catalogue, "
        "the cut ones keep a flagged NaN row");
  {
    const DistanceTable table(0.3, 0.7, -1., 0., 0., 1.5, 4000);
    const std::vector<double> sky = sky_points(900, 0.2, 0.7, rng);
    const std::vector<double> randomsSky = sky_points(3600, 0.2, 0.7, rng);
    Config config;
    config.nRealizations = 2;
    config.seed = 11;
    config.verbose = false;
    const RedshiftCut cut {0.3, 0.6};

    const Result withCut = reconstructLightcone(sky, randomsSky, 800., 1, table, config, cut);
    const std::vector<std::size_t> keep = inside(sky, 0.3, 0.6), keepRandoms = inside(randomsSky, 0.3, 0.6);
    const Result kept = reconstructLightcone(rows(sky, keep), rows(randomsSky, keepRandoms), 800., 1,
                                             table, config);

    const std::size_t n = 900, m = keep.size();
    check(withCut.nObjects == n && withCut.outsideRedshiftCut.size() == n, "one row per input tracer");
    bool flags = true, keptSame = true, cutNaN = true;
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const bool isKept = k < m && keep[k] == i;
      flags = flags && withCut.outsideRedshiftCut[i] == (isKept ? 0 : 1);
      for (unsigned rec = 0; rec < 2; ++rec) {
        const std::size_t to = (std::size_t)rec * n + i;
        if (isKept) {
          const std::size_t from = (std::size_t)rec * m + k;
          keptSame = keptSame && withCut.valid[to] == kept.valid[from];
          for (int c = 0; c < 3; ++c)
            keptSame = keptSame && same_double(withCut.displacement[3*to+c], kept.displacement[3*from+c]) &&
                       same_double(withCut.matchedRandom[3*to+c], kept.matchedRandom[3*from+c]);
        }
        else {
          cutNaN = cutNaN && withCut.valid[to] == 0;
          for (int c = 0; c < 3; ++c)
            cutNaN = cutNaN && std::isnan(withCut.displacement[3*to+c]) && std::isnan(withCut.matchedRandom[3*to+c]);
        }
      }
      if (isKept) {
        keptSame = keptSame && withCut.validRealizations[i] == kept.validRealizations[k];
        for (int c = 0; c < 3; ++c)
          keptSame = keptSame && same_double(withCut.meanDisplacement[3*i+c], kept.meanDisplacement[3*k+c]);
        ++k;
      }
      else {
        cutNaN = cutNaN && withCut.validRealizations[i] == 0 && std::isnan(withCut.meanDisplacement[3*i]);
      }
    }
    check(m > 300 && m < 700, "the cut keeps part of the tracers");
    check(flags, "the flag marks exactly the tracers outside [0.3, 0.6]");
    check(keptSame, "the kept tracers' rows are, bit for bit, those of the catalogue cut beforehand");
    check(cutNaN, "a cut tracer has NaN rows, no valid realization and a NaN mean");

    const std::vector<double> cart = toCartesian(sky, table), randoms = toCartesian(randomsSky, table);
    const Result cartesian = reconstructLightcone(cart, randoms, sky, randomsSky, 800., 1, table, config, cut);
    check(same_bits(cartesian.displacement, withCut.displacement) &&
          cartesian.outsideRedshiftCut == withCut.outsideRedshiftCut,
          "the Cartesian overload cuts the same rows and gives the same result");

    const Result none = reconstructLightcone(sky, randomsSky, 800., 1, table, config);
    const Result open = reconstructLightcone(sky, randomsSky, 800., 1, table, config, RedshiftCut{});
    check(same_bits(none.displacement, open.displacement) &&
          none.outsideRedshiftCut == std::vector<std::uint8_t>(n, 0) &&
          open.outsideRedshiftCut == none.outsideRedshiftCut,
          "the default cuts nothing, and flags nothing");

    std::vector<double> far = sky;
    far[3*5+2] = 5.;
    const Result farCut = reconstructLightcone(far, randomsSky, 800., 1, table, config, cut);
    check(farCut.outsideRedshiftCut[5] == 1, "a cut tracer beyond the distance table is accepted, never converted");

    check_throws([&] { reconstructLightcone(sky, randomsSky, 800., 1, table, config, RedshiftCut{kNaN, 1.}); },
                 "a NaN bound raises");
    check_throws([&] { reconstructLightcone(sky, randomsSky, 800., 1, table, config, RedshiftCut{0.6, 0.3}); },
                 "min > max raises");
    check(throws_naming([&] { reconstructLightcone(sky, randomsSky, 800., 1, table, config, RedshiftCut{0.3, 0.301}); },
                        "tracers"), "a cut keeping too few tracers raises, with the counts");
    std::vector<double> skewed = sky_points(2400, 0.6000001, 0.7, rng);
    const std::vector<double> some = sky_points(1200, 0.3, 0.6, rng);
    skewed.insert(skewed.end(), some.begin(), some.end());
    Config three = config;
    three.nRealizations = 3;
    check(throws_naming([&] { reconstructLightcone(sky, skewed, 800., 1, table, three, cut); }, "randoms"),
          "a cut keeping too few randoms raises, with the counts");

    group("rejectMaskCrossings leaves cut tracers alone, and refuses a NaN anywhere else");
    const std::string maskFile = temporary("full.fits");
    write_full_mask(maskFile);
    const Mask mask(maskFile);
    std::remove(maskFile.c_str());
    Result filtered = withCut;
    rejectMaskCrossings(filtered, mask);
    check(filtered.valid == withCut.valid && filtered.outsideRedshiftCut == withCut.outsideRedshiftCut,
          "a full-sky mask rejects nothing, and the flags are kept");
    const std::size_t aCut = (std::size_t)(std::find(withCut.outsideRedshiftCut.begin(),
                                                     withCut.outsideRedshiftCut.end(), 1) -
                                           withCut.outsideRedshiftCut.begin());
    Result validCut = withCut;
    validCut.valid[aCut] = 1;
    check_throws([&] { rejectMaskCrossings(validCut, mask); }, "a cut tracer with a valid entry is malformed");
    Result nanKept = withCut;
    nanKept.displacement[3 * keep[0]] = kNaN;
    check_throws([&] { rejectMaskCrossings(nanKept, mask); }, "a NaN in a kept tracer's row is malformed");
    Result badFlags = withCut;
    badFlags.outsideRedshiftCut.pop_back();
    check_throws([&] { rejectMaskCrossings(badFlags, mask); }, "a flag array of the wrong size is malformed");

    group("without a cut, invalid realizations and tracers with no valid realization are handled as before");
    for (const bool emptyFlags : {false, true}) {
      // Object 0 loses one realization, object 1 both: its mean is NaN, its
      // displacement rows stay finite, as after any filter.
      Result invalid = none;
      if (emptyFlags) invalid.outsideRedshiftCut.clear();
      invalid.valid[0] = 0;
      invalid.valid[1] = 0;
      invalid.valid[n + 1] = 0;
      internal::summarize(invalid);
      const std::string label = emptyFlags ? " (flags left empty)" : " (flags all 0)";
      check(invalid.validRealizations[1] == 0 && std::isnan(invalid.meanDisplacement[3]),
            "the tracer without a valid realization has a NaN mean" + label);

      Result filtered = invalid;
      bool threw = false;
      try {
        rejectMaskCrossings(filtered, mask);
      }
      catch (const Error&) {
        threw = true;
      }
      check(!threw && filtered.valid == invalid.valid && std::isnan(filtered.meanDisplacement[3]) &&
            filtered.validRealizations == invalid.validRealizations,
            "it is accepted and left as it was" + label);

      Result nanInvalid = invalid;
      nanInvalid.displacement[3 * 1] = kNaN;
      check_throws([&] { rejectMaskCrossings(nanInvalid, mask); },
                   "a NaN in an uncut row is refused even where the realization is invalid, as before" + label);
    }

    group("the correction leaves cut tracers out: uncorrected, diagnostics 0, never a neighbour");
    const std::vector<double> zb {0.2, 0.8}, b {1.2, 1.8};
    const RealSpaceCatalog all = realSpaceLightcone(withCut, sky, table, zb, b, 10.);
    const RealSpaceCatalog sub = realSpaceLightcone(kept, rows(sky, keep), table, zb, b, 10.);
    bool same = true;
    k = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (k < m && keep[k] == i) {
        for (int c = 0; c < 3; ++c) same = same && same_double(all.positions[3*i+c], sub.positions[3*k+c]);
        same = same && all.nNeighbours[i] == sub.nNeighbours[k] &&
               all.nRealizationsAveraged[i] == sub.nRealizationsAveraged[k];
        ++k;
      }
      else {
        same = same && std::isnan(all.positions[3*i+2]) && all.nNeighbours[i] == 0 &&
               std::binary_search(all.uncorrected.begin(), all.uncorrected.end(), i);
      }
    }
    check(same, "the kept tracers are corrected exactly as the catalogue cut beforehand; the cut ones are listed");
    check(all.uncorrected.size() == (n - m) + sub.uncorrected.size(), "and nothing else is uncorrected");
    const RealSpaceCatalog beyond = realSpaceLightcone(farCut, far, table, zb, b, 10.);
    check(std::isnan(beyond.positions[3*5+2]), "a cut tracer beyond the distance table is not converted");
  }

  group("the selection report: one line per array listing the selections applied, one for the "
        "crossings; empty when nothing was applied");
  {
    SelectionCounts c;
    check(c.message().empty(), "no selection, no report");
    c.tracers = 100;
    c.randoms = 400;
    c.redshiftCutApplied = true;
    c.redshiftCut = RedshiftCut{0.9, 1.08};
    c.tracersOutsideRedshiftCut = 7;
    c.randomsOutsideRedshiftCut = 30;
    check(c.message() == "otswap: kept 93 of 100 tracers: 7 outside the redshift cut [0.9, 1.08]\n"
                         "otswap: kept 370 of 400 randoms: 30 outside the redshift cut [0.9, 1.08]\n",
          "the cut alone: " + c.message());
    c.redshiftCutApplied = false;
    c.tracersOutsideRedshiftCut = c.randomsOutsideRedshiftCut = 0;
    c.maskApplied = true;
    c.tracersOutsideMask = 5;
    c.randomsOutsideMask = 20;
    check(c.message() == "otswap: kept 95 of 100 tracers: 5 outside the mask\n"
                         "otswap: kept 380 of 400 randoms: 20 outside the mask\n",
          "the mask alone: " + c.message());
    c.tracersOutsideRedshiftCut = 7;
    c.randomsOutsideRedshiftCut = 30;
    c.tracersOutsideBoth = 2;
    c.redshiftCutApplied = true;
    c.redshiftCut = RedshiftCut{-kInf, 1e-5};
    c.crossingsRejected = true;
    c.maxUnobservedPixelsCrossed = 2;
    c.displacements = 180;
    c.displacementsCrossingMask = 11;
    check(c.message() == "otswap: kept 90 of 100 tracers: 7 outside the redshift cut [-inf, 1e-05], "
                         "5 outside the mask (2 outside both)\n"
                         "otswap: kept 350 of 400 randoms: 30 outside the redshift cut [-inf, 1e-05], "
                         "20 outside the mask (0 outside both)\n"
                         "otswap: rejected 11 of 180 displacements crossing more than 2 unobserved pixels\n",
          "both, an open bound, and the crossings: " + c.message());
  }

  {
    const DistanceTable table(0.3, 0.7, -1., 0., 0., 1.5, 4000);
    const std::vector<double> sky = sky_points(900, 0.2, 0.7, rng);
    const std::vector<double> randomsSky = sky_points(3600, 0.2, 0.7, rng);
    Config config;
    config.nRealizations = 2;
    config.seed = 11;
    config.verbose = false;
    const RedshiftCut cut {0.3, 0.6};

    // NSIDE 64, RING, one pixel in nine unobserved, with values of every
    // kind.
    std::vector<double> values(12 * 64 * 64, 1.);
    for (std::size_t p = 0; p < values.size(); p += 9) values[p] = (p % 2) ? kNaN : -1.;
    const Mask mask(values, PixelOrdering::Ring);

    const std::size_t n = 900;
    std::vector<std::size_t> keep, keepRandoms;
    std::vector<std::uint8_t> expectCut(n), expectMask(n);
    std::size_t both = 0;
    for (std::size_t i = 0; i < n; ++i) {
      expectCut[i] = !(sky[3*i+2] >= cut.min && sky[3*i+2] <= cut.max);
      expectMask[i] = !mask.allows(sky[3*i], sky[3*i+1]);
      both += expectCut[i] && expectMask[i];
      if (!expectCut[i] && !expectMask[i]) keep.push_back(i);
    }
    std::size_t randomsOutsideMask = 0;
    for (std::size_t i = 0; i < randomsSky.size() / 3; ++i) {
      const bool outsideCut = !(randomsSky[3*i+2] >= cut.min && randomsSky[3*i+2] <= cut.max);
      const bool outsideMask = !mask.allows(randomsSky[3*i], randomsSky[3*i+1]);
      randomsOutsideMask += outsideMask;
      if (!outsideCut && !outsideMask) keepRandoms.push_back(i);
    }
    const std::size_t m = keep.size();

    Config unfiltered = config;
    unfiltered.rejectCrossings = false;
    const Result masked = reconstructLightcone(sky, randomsSky, mask, 1, table, config, cut);
    const Result plain = reconstructLightcone(sky, randomsSky, mask, 1, table, unfiltered, cut);
    const Result kept = reconstructLightcone(rows(sky, keep), rows(randomsSky, keepRandoms),
                                             mask.skyAreaDeg2(), 1, table, config);

    group("mask in reconstructLightcone: both selections on every object, the flags independent");
    check(masked.outsideMask == expectMask && masked.outsideRedshiftCut == expectCut,
          "outsideMask is !Mask::allows and outsideRedshiftCut the cut, each on every tracer");
    check(both > 0 && std::count(expectMask.begin(), expectMask.end(), 1) > (std::ptrdiff_t)both &&
          m > 300 && m < 700, "some tracers carry both flags, some one, and part is kept");

    group("the kept tracers are reconstructed exactly as the catalogue selected beforehand, with "
          "the mask's area; the rest keep a NaN row");
    bool keptSame = true, outNaN = true;
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const bool isKept = k < m && keep[k] == i;
      for (unsigned rec = 0; rec < 2; ++rec) {
        const std::size_t to = (std::size_t)rec * n + i;
        if (isKept) {
          const std::size_t from = (std::size_t)rec * m + k;
          keptSame = keptSame && plain.valid[to] == kept.valid[from];
          for (int c = 0; c < 3; ++c)
            keptSame = keptSame && same_double(plain.displacement[3*to+c], kept.displacement[3*from+c]) &&
                       same_double(plain.matchedRandom[3*to+c], kept.matchedRandom[3*from+c]);
        }
        else {
          outNaN = outNaN && plain.valid[to] == 0 && masked.valid[to] == 0;
          for (int c = 0; c < 3; ++c)
            outNaN = outNaN && std::isnan(plain.displacement[3*to+c]) && std::isnan(plain.matchedRandom[3*to+c]);
        }
      }
      if (isKept) ++k;
      else outNaN = outNaN && plain.validRealizations[i] == 0 && std::isnan(plain.meanDisplacement[3*i]);
    }
    check(keptSame, "the kept rows are, bit for bit, those of the selected catalogue");
    check(outNaN, "a tracer left out has NaN rows, no valid realization and a NaN mean");

    group("with rejectCrossings, the result is filtered as rejectMaskCrossings filters it, with "
          "the threshold of Config; without, it is not");
    Result byHand = plain;
    rejectMaskCrossings(byHand, mask, 0);
    check(masked.valid == byHand.valid && masked.filteredNside == 64 && plain.filteredNside == 0,
          "the filter is applied once, with the mask, and recorded");
    check(std::count(masked.valid.begin(), masked.valid.end(), 1) <
          std::count(plain.valid.begin(), plain.valid.end(), 1),
          "it rejects some displacements here");
    Config two = config;
    two.maxUnobservedPixelsCrossed = 2;
    Result byHandTwo = plain;
    rejectMaskCrossings(byHandTwo, mask, 2);
    check(reconstructLightcone(sky, randomsSky, mask, 1, table, two, cut).valid == byHandTwo.valid,
          "maxUnobservedPixelsCrossed is the threshold");

    group("the counts in Result::selection");
    const SelectionCounts& s = masked.selection;
    const std::size_t validBefore = (std::size_t)std::count(plain.valid.begin(), plain.valid.end(), 1);
    const std::size_t validAfter = (std::size_t)std::count(masked.valid.begin(), masked.valid.end(), 1);
    check(s.redshiftCutApplied && s.redshiftCut.min == 0.3 && s.redshiftCut.max == 0.6 && s.maskApplied,
          "the selections applied");
    check(s.tracers == n &&
          s.tracersOutsideRedshiftCut == (std::size_t)std::count(expectCut.begin(), expectCut.end(), 1) &&
          s.tracersOutsideMask == (std::size_t)std::count(expectMask.begin(), expectMask.end(), 1) &&
          s.tracersOutsideBoth == both &&
          s.tracers - (s.tracersOutsideRedshiftCut + s.tracersOutsideMask - s.tracersOutsideBoth) == m,
          "the tracer counts");
    check(s.randoms == 3600 && s.randomsOutsideMask == randomsOutsideMask &&
          s.randoms - (s.randomsOutsideRedshiftCut + s.randomsOutsideMask - s.randomsOutsideBoth) ==
            keepRandoms.size(),
          "the random counts");
    check(s.crossingsRejected && s.maxUnobservedPixelsCrossed == 0 && s.displacements == 2 * m &&
          validBefore == 2 * m && s.displacementsCrossingMask == validBefore - validAfter,
          "the crossing counts");
    check(!plain.selection.crossingsRejected && plain.selection.displacements == 0,
          "no crossing count without the filter");

    group("the Cartesian mask overload drops the same rows and gives the same result");
    const Result cartesian = reconstructLightcone(toCartesian(sky, table), toCartesian(randomsSky, table),
                                                  sky, randomsSky, mask, 1, table, config, cut);
    check(same_bits(cartesian.displacement, masked.displacement) && cartesian.valid == masked.valid &&
          cartesian.outsideMask == masked.outsideMask &&
          cartesian.outsideRedshiftCut == masked.outsideRedshiftCut &&
          cartesian.selection.message() == masked.selection.message(),
          "bit for bit, flags and counts included");

    group("Config::verbose writes Result::selection.message() to std::clog, and nothing without a "
          "selection");
    Config loud = config;
    loud.verbose = true;
    Result reported;
    const std::string text = captured_clog([&] {
      reported = reconstructLightcone(sky, randomsSky, mask, 1, table, loud, cut);
    });
    check(text == reported.selection.message() && lines(text) == 3, "three lines: " + text);
    check(captured_clog([&] { reconstructLightcone(sky, randomsSky, mask, 1, table, config, cut); }).empty(),
          "verbose false writes nothing");
    check(captured_clog([&] { reconstructLightcone(sky, randomsSky, 800., 1, table, loud); }).empty(),
          "the area overload without a cut writes nothing");
    const std::string cutText = captured_clog([&] {
      reconstructLightcone(sky, randomsSky, 800., 1, table, loud, cut);
    });
    check(lines(cutText) == 2 && cutText.find("outside the redshift cut [0.3, 0.6]") != std::string::npos &&
          cutText.find("mask") == std::string::npos,
          "the cut alone reports in two lines: " + cutText);
    const Result unmasked = reconstructLightcone(sky, randomsSky, 800., 1, table, config);
    check(unmasked.outsideMask == std::vector<std::uint8_t>(n, 0) && !unmasked.selection.maskApplied &&
          !unmasked.selection.redshiftCutApplied && unmasked.selection.message().empty(),
          "the area overload flags nothing outside the mask, and applies nothing");
    const Result box = reconstructBox(toCartesian(sky, table), {}, 20., config);
    check(box.outsideMask == std::vector<std::uint8_t>(n, 0) && box.selection.tracers == 0 &&
          box.selection.message().empty(), "a box result: no flag, all counts 0");

    group("too few objects left by the mask raise, naming the selections and the counts");
    std::vector<double> mostlyUnobserved(12 * 64 * 64, 0.);
    for (std::size_t p = 0; p < mostlyUnobserved.size(); p += 20) mostlyUnobserved[p] = 1.;
    const Mask sparse(mostlyUnobserved, PixelOrdering::Ring);
    check(throws_naming([&] { reconstructLightcone(sky, randomsSky, sparse, 1, table, config); },
                        "the mask keeps"), "the mask alone");
    check(throws_naming([&] { reconstructLightcone(sky, randomsSky, sparse, 1, table, config, cut); },
                        "the redshift cut [0.300000, 0.600000] and the mask keeps"), "both");

    group("rejectMaskCrossings leaves masked tracers alone, and refuses one with a valid entry");
    Result again = masked;
    rejectMaskCrossings(again, mask, 0);
    check(again.valid == masked.valid && again.outsideMask == masked.outsideMask,
          "filtering again with the same mask and threshold changes nothing");
    std::size_t aMasked = n;
    for (std::size_t i = 0; i < n && aMasked == n; ++i)
      if (expectMask[i] && !expectCut[i]) aMasked = i;
    Result validMasked = masked;
    validMasked.valid[aMasked] = 1;
    check(throws_naming([&] { rejectMaskCrossings(validMasked, mask); }, "outside the redshift cut or the mask"),
          "a masked tracer with a valid entry is malformed");
    Result badMask = masked;
    badMask.outsideMask.pop_back();
    check(throws_naming([&] { rejectMaskCrossings(badMask, mask); }, "outsideMask holds"),
          "an outsideMask of the wrong size is malformed");

    group("the correction leaves masked tracers out: uncorrected, diagnostics 0, never a neighbour");
    const std::vector<double> zb {0.2, 0.8}, b {1.2, 1.8};
    const RealSpaceCatalog all = realSpaceLightcone(plain, sky, table, zb, b, 10.);
    const RealSpaceCatalog sub = realSpaceLightcone(kept, rows(sky, keep), table, zb, b, 10.);
    bool same = true;
    k = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (k < m && keep[k] == i) {
        for (int c = 0; c < 3; ++c) same = same && same_double(all.positions[3*i+c], sub.positions[3*k+c]);
        same = same && all.nNeighbours[i] == sub.nNeighbours[k];
        ++k;
      }
      else {
        same = same && std::isnan(all.positions[3*i+2]) && all.nNeighbours[i] == 0 &&
               std::binary_search(all.uncorrected.begin(), all.uncorrected.end(), i);
      }
    }
    check(same, "the kept tracers are corrected exactly as the selected catalogue; the others are listed");
    Result maskedHasValid = plain;
    maskedHasValid.validRealizations[aMasked] = 1;
    check(throws_naming([&] { realSpaceLightcone(maskedHasValid, sky, table, zb, b, 10.); },
                        "outside the redshift cut or the mask"),
          "a masked tracer with valid realizations is malformed");

    group("lagrangianSky: toSky of the tracer plus meanDisplacement, NaN rows for the flagged and those "
          "without a valid realization; emptied by a filter that rejects");
    const double kPi = 3.14159265358979323846;
    const std::vector<double> tracerCart = toCartesian(sky, table);
    const auto lagrangian = [&] (const Result& r) {
      std::vector<double> at(3 * n, kNaN);
      for (std::size_t i = 0; i < n; ++i)
        if (!internal::excluded(r, i) && r.validRealizations[i] > 0)
          for (int c = 0; c < 3; ++c) at[3*i+c] = tracerCart[3*i+c] + r.meanDisplacement[3*i+c];
      return toSky(at, table);
    };
    check(same_bits(masked.lagrangianSky, lagrangian(masked)) &&
          same_bits(plain.lagrangianSky, lagrangian(plain)) &&
          same_bits(unmasked.lagrangianSky, lagrangian(unmasked)),
          "bit for bit, with the mask and the filter, without the filter, and with the area");
    bool rowsRight = true, someNaN = false;
    for (std::size_t i = 0; i < n; ++i) {
      const double* s = &masked.lagrangianSky[3*i];
      const bool none = internal::excluded(masked, i) || masked.validRealizations[i] == 0;
      someNaN = someNaN || (none && !internal::excluded(masked, i));
      rowsRight = rowsRight && (none ? std::isnan(s[0]) && std::isnan(s[1]) && std::isnan(s[2])
                                     : s[0] >= 0. && s[0] < 2. * kPi && std::isfinite(s[1]) &&
                                       std::isfinite(s[2]));
    }
    check(rowsRight && someNaN, "NaN rows exactly for the flagged tracers and those the filter left "
          "without a valid realization; finite rows elsewhere, right ascension in [0, 2 pi)");
    check(same_bits(cartesian.lagrangianSky, masked.lagrangianSky),
          "the Cartesian mask overload gives the same coordinates");
    check(box.lagrangianSky.empty(), "a box result has none");
    Result refiltered = plain;
    rejectMaskCrossings(refiltered, mask, 0);
    check(refiltered.lagrangianSky.empty(), "a filter that rejects a displacement empties it");
    Result unchanged = masked;
    rejectMaskCrossings(unchanged, mask, 0);
    check(same_bits(unchanged.lagrangianSky, masked.lagrangianSky), "one that rejects nothing keeps it");
    Result badSky = masked;
    badSky.lagrangianSky.pop_back();
    check(throws_naming([&] { rejectMaskCrossings(badSky, mask); }, "lagrangianSky holds"),
          "a lagrangianSky of the wrong size is malformed for the filter");
    check(throws_naming([&] { realSpaceLightcone(badSky, sky, table, zb, b, 10.); }, "lagrangianSky holds"),
          "and for the correction");
    Result noSky = masked;
    noSky.lagrangianSky.clear();
    check(same_bits(realSpaceLightcone(noSky, sky, table, zb, b, 10.).positions,
                    realSpaceLightcone(masked, sky, table, zb, b, 10.).positions),
          "an empty one is accepted, and the correction does not read it");

    group("a declination outside [-pi/2, pi/2] is refused, with or without a mask; pi/2 is accepted");
    const double above = std::nextafter(kPi / 2., 2.);
    std::vector<double> bad = sky;
    bad[3*7+1] = above;
    std::vector<double> badRandoms = randomsSky;
    badRandoms[3*9+1] = -above;
    const std::vector<double> cart = toCartesian(sky, table), cartRandoms = toCartesian(randomsSky, table);
    check(throws_naming([&] { reconstructLightcone(bad, randomsSky, 800., 1, table, config); },
                        "the tracer sky array holds a declination of 1.570796 at object 7, outside [-pi/2, pi/2]"),
          "the sky overload, tracers");
    check(throws_naming([&] { reconstructLightcone(sky, badRandoms, 800., 1, table, config); },
                        "the random sky array holds a declination of -1.570796 at object 9"),
          "the sky overload, randoms");
    check(throws_naming([&] { reconstructLightcone(cart, cartRandoms, bad, randomsSky, 800., 1, table, config); },
                        "the tracer sky array holds a declination"), "the Cartesian overload");
    check(throws_naming([&] { reconstructLightcone(sky, badRandoms, mask, 1, table, config); },
                        "the random sky array holds a declination"), "the sky mask overload");
    check(throws_naming([&] { reconstructLightcone(cart, cartRandoms, sky, badRandoms, mask, 1, table, config); },
                        "the random sky array holds a declination"), "the Cartesian mask overload");
    check(throws_naming([&] { toCartesian(bad, table); }, "the sky array holds a declination"), "toCartesian");
    check(throws_naming([&] { realSpaceLightcone(plain, bad, table, zb, b, 10.); },
                        "the tracer sky array holds a declination"), "realSpaceLightcone");
    const std::vector<double> poles {1., kPi / 2., 0.5, 2., -kPi / 2., 0.5};
    bool accepted = true;
    try {
      toCartesian(poles, table);
    }
    catch (const Error&) {
      accepted = false;
    }
    check(accepted, "the poles themselves are accepted");
  }

  {
    group("lagrangianSky: a mean position below the table's smallest distance has a NaN redshift and "
          "keeps its right ascension and declination");
    // Every object in a shell 0.2 Mpc/h thick at the table's lower end: the
    // mean of matched randoms lies on a chord, closer to the observer.
    const DistanceTable shell(0.3, 0.7, -1., 0., 0.5, 0.6, 2000);
    const DistanceTable fromZero(0.3, 0.7, -1., 0., 0., 0.6, 2000);
    const std::vector<double> sky = sky_points(600, 0.5, 0.5001, rng);
    const std::vector<double> randomsSky = sky_points(2400, 0.5, 0.5001, rng);
    Config config;
    config.nRealizations = 4;
    config.seed = 5;
    config.verbose = false;
    const Result r = reconstructLightcone(sky, randomsSky, 800., 1, shell, config);
    const std::vector<double> cart = toCartesian(sky, shell);
    std::size_t below = 0;
    bool kept = true, raises = true;
    for (std::size_t i = 0; i < 600; ++i) {
      const std::vector<double> at {cart[3*i] + r.meanDisplacement[3*i],
                                    cart[3*i+1] + r.meanDisplacement[3*i+1],
                                    cart[3*i+2] + r.meanDisplacement[3*i+2]};
      const std::vector<double> s = toSky(at, fromZero);
      const double* l = &r.lagrangianSky[3*i];
      if (s[2] < 0.5) {
        ++below;
        kept = kept && std::isnan(l[2]) && same_double(l[0], s[0]) && same_double(l[1], s[1]);
        raises = raises && throws_naming([&] { toSky(at, shell); }, "covering a wider redshift range");
      }
      else {
        kept = kept && same_double(l[0], s[0]) && same_double(l[1], s[1]) && std::isfinite(l[2]);
      }
    }
    check(below > 0 && below < 600, std::to_string(below) + " of 600 below the table");
    check(kept, "their redshift NaN, right ascension and declination as toSky's; the rest finite");
    check(raises, "toSky itself raises for them, suggesting a wider table");
  }

  return report("test_rsd");
}
