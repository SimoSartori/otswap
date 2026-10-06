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
 *  @file src/Mask.cpp
 *
 *  @brief The Healpix veto mask, and the displacement filter built on it.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include <fitsio.h>

#include <healpix_base.h>
#include <healpix_map.h>

#include "arc.h"
#include "detmath.h"
#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;

  const double kFullSkyDeg2 = 4. * kPi * (180. / kPi) * (180. / kPi);

  // Pixel values are read 1024 to a FITS row; the layout is checked
  // against the pixel count before anything is read.
  constexpr long kPixelsPerRow = 1024;

  // A sub-arc shorter than this fraction of the pixel size is not split
  // further.
  constexpr double kArcFloorInPixels = 1.e-6;

  // Healpix's loc2pix is a protected member of T_Healpix_Base. Named
  // through a derived class it yields a pointer to member, which applies to
  // any Healpix_Base, so the vendored files need no change.
  struct Locator : Healpix_Base {
    static int pixel (const Healpix_Base& base, const double z, const double phi,
                      const double sth, const bool haveSth)
    {
      return (base.*(&Locator::loc2pix))(z, phi, sth, haveSth);
    }
  };

  // True when the sub-arc from a point in pixP to a point in pixQ, shorter
  // than a pixel, needs no further splitting: pixQ shares an edge with
  // pixP, and no unobserved pixel other than pixQ and the two endpoint
  // pixels neighbours both.
  bool edge_crossing_safe (const Healpix_Map<float>& map, const int pixP, const int pixQ,
                           const int pixA, const int pixB)
  {
    fix_arr<int, 8> nbP, nbQ;
    map.neighbors(pixP, nbP);
    if (nbP[0] != pixQ && nbP[2] != pixQ && nbP[4] != pixQ && nbP[6] != pixQ) return false;

    map.neighbors(pixQ, nbQ);
    for (int i = 0; i < 8; ++i) {
      const int r = nbP[i];
      if (r < 0 || r == pixQ || r == pixA || r == pixB ||
          otswap::internal::pixel_observed(map, r)) continue;
      for (int j = 0; j < 8; ++j)
        if (nbQ[j] == r) return false;
    }
    return true;
  }

  struct ArcSearch {
    const Healpix_Map<float>& map;
    int pixA, pixB;
    double floor2, pixel2;
    unsigned limit;
    std::vector<int>& found;
  };

  // Examines the sub-arc from p to q, whose end pixels have been examined
  // already. Returns true once more than search.limit pixels are found.
  bool segment_exceeds (const ArcSearch& search,
                        const vec3& p, const int pixP, const vec3& q, const int pixQ)
  {
    if (pixP == pixQ) return false;
    const double d2 = (p - q).SquaredLength();
    if (d2 < search.floor2) return false;
    if (d2 < search.pixel2 &&
        edge_crossing_safe(search.map, pixP, pixQ, search.pixA, search.pixB)) return false;

    const vec3 m = (p + q).Norm();
    const int pixM = otswap::internal::vec2pix(search.map, m);
    if (pixM != search.pixA && pixM != search.pixB &&
        !otswap::internal::pixel_observed(search.map, pixM) &&
        std::find(search.found.begin(), search.found.end(), pixM) == search.found.end()) {
      search.found.push_back(pixM);
      if (search.found.size() > search.limit) return true;
    }

    return segment_exceeds(search, p, pixP, m, pixM) || segment_exceeds(search, m, pixM, q, pixQ);
  }

}


// ============================================================================


class otswap::Mask::Impl {

public:

  Healpix_Map<float> map;
  long allowed = 0;

  bool allows (const double ra, const double dec) const
  {
    if (!std::isfinite(ra))
      throw Error("the right ascension is " + std::to_string(ra) + "; it must be finite");
    if (!(std::fabs(dec) <= kPi/2.))
      throw Error("the declination is " + std::to_string(dec) +
                  " radians; it must lie in [-pi/2, pi/2]");

    const double theta = kPi/2. - dec;
    return internal::pixel_observed(map, internal::ang2pix(map, theta, internal::normalize_ra(ra)));
  }

};


// ============================================================================


otswap::Mask::Mask (const std::string& fitsFile)
{
  fitsfile* fptr = nullptr;
  int status = 0;
  int hdutype = 0;
  int nside = 0;

  auto fail = [&] (const std::string& action) {
    char text[FLEN_STATUS];
    fits_get_errstatus(status, text);
    int closeStatus = 0;
    if (fptr != nullptr) fits_close_file(fptr, &closeStatus);
    throw Error(action + " in the mask file " + fitsFile + ": " + text +
                " (cfitsio status " + std::to_string(status) + ")");
  };

  if (fits_open_file(&fptr, fitsFile.c_str(), READONLY, &status))
    fail("cannot open");

  if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) fail("cannot move to HDU 2");

  char ordering[FLEN_VALUE] = "RING";

  if (fits_read_key(fptr, TINT, "NSIDE", &nside, NULL, &status)) {
    if (status != KEY_NO_EXIST) fail("cannot read NSIDE from HDU 2");
    status = 0;
    if (fits_movabs_hdu(fptr, 1, &hdutype, &status)) fail("cannot move to HDU 1");
    if (fits_read_key(fptr, TINT, "NSIDE", &nside, NULL, &status)) {
      if (status != KEY_NO_EXIST) fail("cannot read NSIDE from HDU 1");
      status = 0;
    }
  }

  if (fits_read_key(fptr, TSTRING, "ORDERING", ordering, NULL, &status)) {
    if (status != KEY_NO_EXIST) fail("cannot read ORDERING");
    status = 0;
    std::strcpy(ordering, "RING");
  }

  if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) fail("cannot move to HDU 2");

  long nrows = 0;
  if (fits_get_num_rows(fptr, &nrows, &status)) fail("cannot read the number of rows");

  const long totalPixels = nrows * kPixelsPerRow;

  if (nside == 0) nside = (int)std::sqrt((double)totalPixels / 12.0);

  if (nside <= 0 || 12L * (long)nside * (long)nside != totalPixels) {
    status = 0;
    int closeStatus = 0;
    if (fptr != nullptr) fits_close_file(fptr, &closeStatus);
    throw Error("the mask file " + fitsFile + " holds " + std::to_string(totalPixels) +
                " pixels over " + std::to_string(nrows) + " rows of " +
                std::to_string(kPixelsPerRow) + ", which is not 12 * NSIDE^2 for NSIDE = " +
                std::to_string(nside) + "; it is not a Healpix map in the expected layout");
  }

  std::vector<float> values((std::size_t)totalPixels);
  for (long row = 1; row <= nrows; ++row)
    if (fits_read_col(fptr, TFLOAT, 1, row, 1, kPixelsPerRow, NULL,
                      &values[(std::size_t)((row-1) * kPixelsPerRow)], NULL, &status))
      fail("cannot read pixel row " + std::to_string(row));

  if (fits_close_file(fptr, &status)) {
    fptr = nullptr;
    fail("cannot close");
  }

  const Healpix_Ordering_Scheme scheme =
    std::strncmp(ordering, "NEST", 4) == 0 ? NEST : RING;

  auto impl = std::make_shared<Impl>();
  impl->map.SetNside(nside, scheme);

  // Any value is accepted. A pixel is observed when its value exceeds 0,
  // fractional values included; 0, negative values, NaN and Healpix's
  // UNSEEN are unobserved. The values are kept as read, and the sky area
  // counts the observed pixels: no value is used as a weight.
  for (long i = 0; i < impl->map.Npix(); ++i) {
    const float v = values[(std::size_t)i];
    impl->map[i] = v;
    if (v > internal::kMaskAllowedAbove) ++impl->allowed;
  }

  m_impl = impl;
}


// ============================================================================


bool otswap::Mask::allows (const double rightAscension, const double declination) const
{
  return m_impl->allows(rightAscension, declination);
}


// ============================================================================


double otswap::Mask::skyAreaDeg2 () const
{
  return kFullSkyDeg2 * (double)m_impl->allowed / (double)m_impl->map.Npix();
}


// ============================================================================


int otswap::Mask::nside () const
{
  return m_impl->map.Nside();
}


// ============================================================================


int otswap::internal::vec2pix (const Healpix_Base& base, const vec3& v)
{
  const double xl = 1./v.Length();
  const double phi = (v.x == 0. && v.y == 0.) ? 0.0 : det_atan2(v.y, v.x);
  const double nz = v.z*xl;
  if (std::abs(nz) > 0.99)
    return Locator::pixel(base, nz, phi, std::sqrt(v.x*v.x + v.y*v.y)*xl, true);
  else
    return Locator::pixel(base, nz, phi, 0, false);
}


// ============================================================================


int otswap::internal::ang2pix (const Healpix_Base& base, const double theta, const double phi)
{
  return ((theta < 0.01) || (theta > 3.14159-0.01)) ?
    Locator::pixel(base, det_cos(theta), phi, det_sin(theta), true) :
    Locator::pixel(base, det_cos(theta), phi, 0., false);
}


// ============================================================================


bool otswap::internal::pixel_observed (const Healpix_Map<float>& map, const int pixel)
{
  return map[pixel] > kMaskAllowedAbove;
}


// ============================================================================


std::size_t otswap::internal::arc_unobserved_pixels (const Healpix_Map<float>& map,
                                                     const vec3& a, const vec3& b,
                                                     const unsigned limit,
                                                     std::vector<int>& found)
{
  found.clear();

  const int pixA = vec2pix(map, a);
  const int pixB = vec2pix(map, b);
  const double pixelSize = std::sqrt(4. * kPi / (double)map.Npix());
  const double arcFloor = kArcFloorInPixels * pixelSize;

  const ArcSearch search {map, pixA, pixB, arcFloor * arcFloor, pixelSize * pixelSize,
                          limit, found};
  segment_exceeds(search, a, pixA, b, pixB);

  return found.size();
}


// ============================================================================


void otswap::rejectMaskCrossings (Result& result, const Mask& mask,
                                  const unsigned maxForbiddenPixels)
{
  const std::size_t nObjects = result.nObjects;
  const unsigned nRealizations = result.nRealizations;
  const std::size_t nDisplacements = (std::size_t)nRealizations * nObjects;
  const std::size_t expected = 3 * nDisplacements;

  if (nObjects == 0 || nRealizations == 0)
    throw Error("the result is empty");

  if (result.displacement.size() != expected || result.matchedRandom.size() != expected)
    throw Error("the result is malformed: displacement and matchedRandom should each hold " +
                std::to_string(expected) + " entries, but hold " +
                std::to_string(result.displacement.size()) + " and " +
                std::to_string(result.matchedRandom.size()));

  if (result.valid.size() != nDisplacements)
    throw Error("the result is malformed: valid holds " + std::to_string(result.valid.size()) +
                " entries against " + std::to_string(nDisplacements) + " displacements");

  for (std::size_t k = 0; k < nDisplacements; ++k)
    if (result.valid[k] > 1)
      throw Error("the result is malformed: valid holds " + std::to_string(result.valid[k]) +
                  " at entry " + std::to_string(k) + "; every entry must be 0 or 1");

  // Tracers outside a redshift cut took no part in the reconstruction:
  // their rows are NaN and never valid, and the filter leaves them alone.
  const std::vector<std::uint8_t>& cut = result.outsideRedshiftCut;
  if (!cut.empty() && cut.size() != nObjects)
    throw Error("the result is malformed: outsideRedshiftCut holds " + std::to_string(cut.size()) +
                " entries for " + std::to_string(nObjects) + " objects");

  for (std::size_t k = 0; k < nDisplacements; ++k) {
    const std::size_t i = k % nObjects;
    if (!cut.empty() && cut[i] != 0) {
      if (result.valid[k] != 0)
        throw Error("the result is malformed: object " + std::to_string(i) + " lies outside the "
                    "redshift cut but has a valid displacement at entry " + std::to_string(k));
      continue;
    }
    for (std::size_t c = 0; c < 3; ++c)
      if (!std::isfinite(result.displacement[3*k+c]) || !std::isfinite(result.matchedRandom[3*k+c]))
        throw Error("the result is malformed: a non-finite coordinate at entry " +
                    std::to_string(3*k+c));
  }

  if (result.filteredNside != 0 && result.filteredNside != mask.nside())
    throw Error("this result was already filtered against a mask of NSIDE " +
                std::to_string(result.filteredNside) + ", and this mask has NSIDE " +
                std::to_string(mask.nside()) + "; a pixel count is not comparable "
                "between the two, so filter an unfiltered result instead");

  const Healpix_Map<float>& map = mask.m_impl->map;

#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t i = 0; i < nObjects; ++i) {

    std::vector<int> found;

    for (unsigned rec = 0; rec < nRealizations; ++rec) {

      const std::size_t d = (std::size_t)rec * nObjects + i;
      if (!result.valid[d]) continue;

      const std::size_t at = 3 * d;

      const double lx = result.matchedRandom[at];
      const double ly = result.matchedRandom[at+1];
      const double lz = result.matchedRandom[at+2];

      const double ex = lx - result.displacement[at];
      const double ey = ly - result.displacement[at+1];
      const double ez = lz - result.displacement[at+2];

      const double r0 = std::sqrt(ex*ex + ey*ey + ez*ez);
      const double r1 = std::sqrt(lx*lx + ly*ly + lz*lz);
      if (r0 == 0. || r1 == 0.) continue;

      const vec3 a(ex/r0, ey/r0, ez/r0);
      const vec3 b(lx/r1, ly/r1, lz/r1);

      if ((a + b).SquaredLength() == 0. ||
          internal::arc_unobserved_pixels(map, a, b, maxForbiddenPixels, found) > maxForbiddenPixels)
        result.valid[d] = 0;
    }
  }

  internal::summarize(result);
  result.filteredNside = mask.nside();
}
