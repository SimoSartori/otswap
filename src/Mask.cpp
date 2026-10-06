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
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <fitsio.h>

#include <healpix_base.h>

#include "arc.h"
#include "detmath.h"
#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;

  const double kFullSkyDeg2 = 4. * kPi * (180. / kPi) * (180. / kPi);

  // Pixel values are read 1024 to a FITS row; the layout is checked
  // against the pixel count before anything is read.
  constexpr LONGLONG kPixelsPerRow = 1024;

  // Rows read at a time: the values of one block, 512 KiB of doubles, are
  // turned into bytes before the next block is read.
  constexpr LONGLONG kRowsPerBlock = 64;

  // A sub-arc shorter than this fraction of the pixel size is not split
  // further.
  constexpr double kArcFloorInPixels = 1.e-6;

  // Healpix's loc2pix is a protected member of T_Healpix_Base. Named
  // through a derived class it yields a pointer to member, which applies to
  // any base of the same index type, so the vendored files need no change.
  template <typename I>
  struct Locator : T_Healpix_Base<I> {
    static I pixel (const T_Healpix_Base<I>& base, const double z, const double phi,
                    const double sth, const bool haveSth)
    {
      return (base.*(&Locator<I>::loc2pix))(z, phi, sth, haveSth);
    }
  };

  using Pixel = std::int64_t;

  // True when the sub-arc from a point in pixP to a point in pixQ, shorter
  // than a pixel, needs no further splitting: pixQ shares an edge with
  // pixP, and no unobserved pixel other than pixQ and the two endpoint
  // pixels neighbours both.
  bool edge_crossing_safe (const otswap::internal::PixelMask& mask, const Pixel pixP,
                           const Pixel pixQ, const Pixel pixA, const Pixel pixB)
  {
    fix_arr<Pixel, 8> nbP, nbQ;
    mask.base.neighbors(pixP, nbP);
    if (nbP[0] != pixQ && nbP[2] != pixQ && nbP[4] != pixQ && nbP[6] != pixQ) return false;

    mask.base.neighbors(pixQ, nbQ);
    for (int i = 0; i < 8; ++i) {
      const Pixel r = nbP[i];
      if (r < 0 || r == pixQ || r == pixA || r == pixB ||
          otswap::internal::pixel_observed(mask, r)) continue;
      for (int j = 0; j < 8; ++j)
        if (nbQ[j] == r) return false;
    }
    return true;
  }

  struct ArcSearch {
    const otswap::internal::PixelMask& mask;
    Pixel pixA, pixB;
    double floor2, pixel2;
    unsigned limit;
    std::vector<Pixel>& found;
  };

  // Examines the sub-arc from p to q, whose end pixels have been examined
  // already. Returns true once more than search.limit pixels are found.
  bool segment_exceeds (const ArcSearch& search,
                        const vec3& p, const Pixel pixP, const vec3& q, const Pixel pixQ)
  {
    if (pixP == pixQ) return false;
    const double d2 = (p - q).SquaredLength();
    if (d2 < search.floor2) return false;
    if (d2 < search.pixel2 &&
        edge_crossing_safe(search.mask, pixP, pixQ, search.pixA, search.pixB)) return false;

    const vec3 m = (p + q).Norm();
    const Pixel pixM = otswap::internal::vec2pix(search.mask.base, m);
    if (pixM != search.pixA && pixM != search.pixB &&
        !otswap::internal::pixel_observed(search.mask, pixM) &&
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

  internal::PixelMask mask;

  bool allows (const double ra, const double dec) const
  {
    if (!std::isfinite(ra))
      throw Error("the right ascension is " + std::to_string(ra) + "; it must be finite");
    if (!(std::fabs(dec) <= kPi/2.))
      throw Error("the declination is " + std::to_string(dec) +
                  " radians; it must lie in [-pi/2, pi/2]");

    const double theta = kPi/2. - dec;
    return internal::pixel_observed(mask, internal::ang2pix(mask.base, theta,
                                                            internal::normalize_ra(ra)));
  }

};


// ============================================================================


otswap::Mask::Mask (const std::string& fitsFile)
{
  const std::string source = "the mask file " + fitsFile;

  fitsfile* fptr = nullptr;
  int status = 0;
  int hdutype = 0;
  LONGLONG nside = 0;

  auto fail = [&] (const std::string& action) {
    char text[FLEN_STATUS];
    fits_get_errstatus(status, text);
    throw Error(action + " in " + source + ": " + text +
                " (cfitsio status " + std::to_string(status) + ")");
  };

  if (fits_open_file(&fptr, fitsFile.c_str(), READONLY, &status))
    fail("cannot open");

  auto impl = std::make_shared<Impl>();

  try {
    if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) fail("cannot move to HDU 2");

    char ordering[FLEN_VALUE] = "RING";

    if (fits_read_key(fptr, TLONGLONG, "NSIDE", &nside, NULL, &status)) {
      if (status != KEY_NO_EXIST) fail("cannot read NSIDE from HDU 2");
      status = 0;
      if (fits_movabs_hdu(fptr, 1, &hdutype, &status)) fail("cannot move to HDU 1");
      if (fits_read_key(fptr, TLONGLONG, "NSIDE", &nside, NULL, &status)) {
        if (status != KEY_NO_EXIST) fail("cannot read NSIDE from HDU 1");
        status = 0;
        nside = 0;
      }
    }

    if (fits_read_key(fptr, TSTRING, "ORDERING", ordering, NULL, &status)) {
      if (status != KEY_NO_EXIST) fail("cannot read ORDERING");
      status = 0;
      std::strcpy(ordering, "RING");
    }

    const bool nested = std::strncmp(ordering, "NEST", 4) == 0;

    if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) fail("cannot move to HDU 2");

    LONGLONG nrows = 0;
    if (fits_get_num_rowsll(fptr, &nrows, &status)) fail("cannot read the number of rows");

    // The most rows a map of NSIDE 2^29 fills; more cannot be a valid map,
    // and the bound keeps the pixel count below overflow.
    const LONGLONG maxRows = 12 * internal::kMaxNside * internal::kMaxNside / kPixelsPerRow;
    if (nrows > maxRows)
      throw Error(source + " holds " + std::to_string(nrows) + " rows of " +
                  std::to_string(kPixelsPerRow) + " pixels, more than a map of NSIDE 2^29 "
                  "fills; it is not a Healpix map in the expected layout");

    const LONGLONG totalPixels = nrows * kPixelsPerRow;

    if (nside == 0) nside = (LONGLONG)std::llround(std::sqrt((double)totalPixels / 12.0));

    if (nside >= 1 && nside <= internal::kMaxNside && 12 * nside * nside != totalPixels)
      throw Error(source + " holds " + std::to_string(totalPixels) +
                  " pixels over " + std::to_string(nrows) + " rows of " +
                  std::to_string(kPixelsPerRow) + ", which is not 12 * NSIDE^2 for NSIDE = " +
                  std::to_string(nside) + "; it is not a Healpix map in the expected layout");

    int typecode = 0;
    LONGLONG repeat = 0, width = 0;
    if (fits_get_coltypell(fptr, 1, &typecode, &repeat, &width, &status))
      fail("cannot read the type of the pixel column");

    internal::set_geometry(impl->mask, nside, nested, source);

    // A floating column is read as stored, NaN and +-inf included: cfitsio's
    // null check would turn +-inf into the null value too. An integer
    // column is read with the check, so that its null value (TNULL) reads
    // as NaN, unobserved.
    const bool floating = typecode == TFLOAT || typecode == TDOUBLE;
    double nullValue = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> block((std::size_t)(kRowsPerBlock * kPixelsPerRow));
    for (LONGLONG row = 1; row <= nrows; row += kRowsPerBlock) {
      const LONGLONG count = std::min(kRowsPerBlock, nrows - row + 1) * kPixelsPerRow;
      int anyNull = 0;
      if (fits_read_col(fptr, TDOUBLE, 1, row, 1, count, floating ? nullptr : &nullValue,
                        block.data(), &anyNull, &status))
        fail("cannot read pixel rows from " + std::to_string(row));
      internal::mark_observed(impl->mask, (row - 1) * kPixelsPerRow, block.data(),
                              (std::size_t)count);
    }

    if (fits_close_file(fptr, &status)) {
      fptr = nullptr;
      fail("cannot close");
    }
  }
  catch (...) {
    if (fptr != nullptr) {
      int closeStatus = 0;
      fits_close_file(fptr, &closeStatus);
    }
    throw;
  }

  m_impl = impl;
}


// ============================================================================


otswap::Mask::Mask (const double* values, const std::size_t count, const PixelOrdering ordering)
{
  if (ordering != PixelOrdering::Ring && ordering != PixelOrdering::Nested)
    throw Error("the pixel ordering is neither PixelOrdering::Ring nor PixelOrdering::Nested");

  const std::string source = "the map of " + std::to_string(count) + " values";

  if (count == 0)
    throw Error("the map is empty; a full-sky Healpix map holds 12 * NSIDE^2 values");
  if (values == nullptr)
    throw Error(source + " is given as a null pointer");

  // NSIDE from the length, in integers: the root of count / 12, if exact.
  const std::uint64_t perFace = (std::uint64_t)count / 12;
  std::uint64_t nside = (std::uint64_t)std::llround(std::sqrt((double)perFace));
  while (nside > 0 && nside * nside > perFace) --nside;
  while ((nside + 1) * (nside + 1) <= perFace) ++nside;

  if (count % 12 != 0 || nside * nside != perFace)
    throw Error(source + " is not 12 * NSIDE^2 for any integer NSIDE: it implies NSIDE = " +
                std::to_string(std::sqrt((double)count / 12.)) +
                "; a full-sky Healpix map holds 12 * NSIDE^2 values");

  auto impl = std::make_shared<Impl>();
  // count / 12 < 2^63, so its root fits an int64 whatever its size.
  internal::set_geometry(impl->mask, (std::int64_t)nside, ordering == PixelOrdering::Nested,
                         source);
  internal::mark_observed(impl->mask, 0, values, count);

  m_impl = impl;
}


// ============================================================================


otswap::Mask::Mask (const std::vector<double>& values, const PixelOrdering ordering)
  : Mask(values.data(), values.size(), ordering)
{}


// ============================================================================


bool otswap::Mask::allows (const double rightAscension, const double declination) const
{
  return m_impl->allows(rightAscension, declination);
}


// ============================================================================


double otswap::Mask::skyAreaDeg2 () const
{
  return kFullSkyDeg2 * (double)m_impl->mask.allowed / (double)m_impl->mask.base.Npix();
}


// ============================================================================


int otswap::Mask::nside () const
{
  return (int)m_impl->mask.base.Nside();
}


// ============================================================================


void otswap::internal::check_nside (const std::int64_t nside, const bool nested,
                                    const std::string& source)
{
  if (nside < 1 || nside > kMaxNside)
    throw Error(source + " has NSIDE = " + std::to_string(nside) + "; NSIDE must lie in [1, " +
                std::to_string(kMaxNside) + "], 2^29 being Healpix's limit");
  if (nested && (nside & (nside - 1)) != 0)
    throw Error(source + " is in NESTED ordering with NSIDE = " + std::to_string(nside) +
                ", which is not a power of 2; the NESTED scheme needs one");
}


// ============================================================================


void otswap::internal::set_geometry (PixelMask& mask, const std::int64_t nside, const bool nested,
                                     const std::string& source)
{
  check_nside(nside, nested, source);
  mask.base.SetNside(nside, nested ? NEST : RING);
  mask.observed.assign((std::size_t)mask.base.Npix(), 0);
  mask.allowed = 0;
}


// ============================================================================


void otswap::internal::mark_observed (PixelMask& mask, const std::int64_t first,
                                      const double* values, const std::size_t count)
{
  std::uint8_t* out = mask.observed.data() + first;
  std::int64_t allowed = 0;
  for (std::size_t i = 0; i < count; ++i) {
    const std::uint8_t observed = values[i] > kMaskAllowedAbove ? 1 : 0;
    out[i] = observed;
    allowed += observed;
  }
  mask.allowed += allowed;
}


// ============================================================================


template <typename I>
I otswap::internal::vec2pix (const T_Healpix_Base<I>& base, const vec3& v)
{
  const double xl = 1./v.Length();
  const double phi = (v.x == 0. && v.y == 0.) ? 0.0 : det_atan2(v.y, v.x);
  const double nz = v.z*xl;
  if (std::abs(nz) > 0.99)
    return Locator<I>::pixel(base, nz, phi, std::sqrt(v.x*v.x + v.y*v.y)*xl, true);
  else
    return Locator<I>::pixel(base, nz, phi, 0, false);
}

template int otswap::internal::vec2pix<int> (const T_Healpix_Base<int>&, const vec3&);
template std::int64_t otswap::internal::vec2pix<std::int64_t> (const T_Healpix_Base<std::int64_t>&,
                                                               const vec3&);


// ============================================================================


template <typename I>
I otswap::internal::ang2pix (const T_Healpix_Base<I>& base, const double theta, const double phi)
{
  return ((theta < 0.01) || (theta > 3.14159-0.01)) ?
    Locator<I>::pixel(base, det_cos(theta), phi, det_sin(theta), true) :
    Locator<I>::pixel(base, det_cos(theta), phi, 0., false);
}

template int otswap::internal::ang2pix<int> (const T_Healpix_Base<int>&, double, double);
template std::int64_t otswap::internal::ang2pix<std::int64_t> (const T_Healpix_Base<std::int64_t>&,
                                                               double, double);


// ============================================================================


std::size_t otswap::internal::arc_unobserved_pixels (const PixelMask& mask,
                                                     const vec3& a, const vec3& b,
                                                     const unsigned limit,
                                                     std::vector<std::int64_t>& found)
{
  found.clear();

  const std::int64_t pixA = vec2pix(mask.base, a);
  const std::int64_t pixB = vec2pix(mask.base, b);
  const double pixelSize = std::sqrt(4. * kPi / (double)mask.base.Npix());
  const double arcFloor = kArcFloorInPixels * pixelSize;

  const ArcSearch search {mask, pixA, pixB, arcFloor * arcFloor, pixelSize * pixelSize,
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

  const internal::PixelMask& pixels = mask.m_impl->mask;

#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t i = 0; i < nObjects; ++i) {

    std::vector<std::int64_t> found;

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
          internal::arc_unobserved_pixels(pixels, a, b, maxForbiddenPixels, found) > maxForbiddenPixels)
        result.valid[d] = 0;
    }
  }

  internal::summarize(result);
  result.filteredNside = mask.nside();
}
