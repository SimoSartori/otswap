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
#include <array>
#include <cmath>
#include <exception>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <fitsio.h>

#include "arc.h"
#include "detmath.h"
#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;

  const double kFullSkyDeg2 = 4. * kPi * (180. / kPi) * (180. / kPi);

  /// Pixel values are read 1024 to a FITS row; the layout is checked
  /// against the pixel count before anything is read.
  constexpr LONGLONG kPixelsPerRow = 1024;

  /// Rows read at a time: the values of one block, 512 KiB of doubles, are
  /// turned into bytes before the next block is read.
  constexpr LONGLONG kRowsPerBlock = 64;

  /// A sub-arc shorter than this fraction of the pixel size is not split
  /// further.
  constexpr double kArcFloorInPixels = 1.e-6;

  using Pixel = std::int64_t;

  struct ArcSearch {
    const otswap::internal::PixelMask& mask;
    Pixel pixA, pixB;
    double floor2, pixel2;
    unsigned limit;
    std::vector<Pixel>& path;
    std::size_t runStart;
    std::size_t run;
    std::size_t longest;
  };

  bool counted (const ArcSearch& search, const Pixel pixel)
  {
    return pixel != search.pixA && pixel != search.pixB &&
           !otswap::internal::pixel_observed(search.mask, pixel);
  }

  /// True when the sub-arc from a point in pixP to a point in pixQ, shorter
  /// than a pixel, needs no further splitting: pixQ shares an edge with
  /// pixP, and no pixel they both neighbour could change the runs. When
  /// both count, any common neighbour could split their run; otherwise
  /// only a counted one could add to a run.
  bool edge_crossing_safe (const ArcSearch& search, const Pixel pixP, const Pixel pixQ)
  {
    std::array<Pixel, 8> nbP, nbQ;
    search.mask.base.neighbors(pixP, nbP);
    if (nbP[0] != pixQ && nbP[2] != pixQ && nbP[4] != pixQ && nbP[6] != pixQ) return false;

    const bool bothCounted = counted(search, pixP) && counted(search, pixQ);
    search.mask.base.neighbors(pixQ, nbQ);
    for (const Pixel r : nbP) {
      if (r < 0 || r == pixQ || (!bothCounted && !counted(search, r))) continue;
      for (const Pixel q : nbQ)
        if (q == r) return false;
    }
    return true;
  }

  bool meet (ArcSearch& search, const Pixel pixel)
  {
    if (pixel == search.path.back()) return false;
    if (!counted(search, pixel)) {
      search.path.push_back(pixel);
      search.runStart = search.path.size();
      search.run = 0;
      return false;
    }
    if (std::find(search.path.begin() + (std::ptrdiff_t)search.runStart, search.path.end(),
                  pixel) == search.path.end())
      search.longest = std::max(search.longest, ++search.run);
    search.path.push_back(pixel);
    return search.run > search.limit;
  }

  /// Examines the sub-arc from p to q, in order along the arc: the pixel of
  /// p has been met already, the pixel of q is met after. Returns true once
  /// a run holds more than search.limit pixels.
  bool segment_exceeds (ArcSearch& search,
                        const otswap::internal::Vec3& p, const Pixel pixP,
                        const otswap::internal::Vec3& q, const Pixel pixQ)
  {
    if (pixP == pixQ) return false;
    const double d2 = (p - q).squared_length();
    if (d2 < search.floor2) return false;
    if (d2 < search.pixel2 && edge_crossing_safe(search, pixP, pixQ)) return false;

    const otswap::internal::Vec3 m = (p + q).normalized();
    const Pixel pixM = otswap::internal::vec2pix(search.mask.base, m);

    return segment_exceeds(search, p, pixP, m, pixM) || meet(search, pixM) ||
           segment_exceeds(search, m, pixM, q, pixQ);
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

  const std::uint64_t perFace = (std::uint64_t)count / 12;
  std::uint64_t nside = (std::uint64_t)std::llround(std::sqrt((double)perFace));
  while (nside > 0 && nside * nside > perFace) --nside;
  while ((nside + 1) * (nside + 1) <= perFace) ++nside;

  if (count % 12 != 0 || nside * nside != perFace)
    throw Error(source + " is not 12 * NSIDE^2 for any integer NSIDE: it implies NSIDE = " +
                std::to_string(std::sqrt((double)count / 12.)) +
                "; a full-sky Healpix map holds 12 * NSIDE^2 values");

  auto impl = std::make_shared<Impl>();
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
  return kFullSkyDeg2 * (double)m_impl->mask.allowed / (double)m_impl->mask.base.npix();
}


// ============================================================================


int otswap::Mask::nside () const
{
  return (int)m_impl->mask.base.nside();
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
  mask.base.set(nside, nested);
  mask.observed.assign((std::size_t)mask.base.npix(), 0);
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


double otswap::internal::arc_pixel_size (const PixelMask& mask)
{
  return std::sqrt(4. * kPi / (double)mask.base.npix());
}


// ============================================================================


std::size_t otswap::internal::arc_longest_unobserved_run (const PixelMask& mask,
                                                         const Vec3& a, const std::int64_t pixA,
                                                         const Vec3& b, const double pixelSize,
                                                         const unsigned limit,
                                                         std::vector<std::int64_t>& path)
{
  const std::int64_t pixB = vec2pix(mask.base, b);
  const double arcFloor = kArcFloorInPixels * pixelSize;

  path.assign(1, pixA);
  ArcSearch search {mask, pixA, pixB, arcFloor * arcFloor, pixelSize * pixelSize,
                    limit, path, 1, 0, 0};
  if (!segment_exceeds(search, a, pixA, b, pixB)) meet(search, pixB);

  return search.longest;
}


// ============================================================================


std::size_t otswap::internal::arc_longest_unobserved_run (const PixelMask& mask,
                                                         const Vec3& a, const Vec3& b,
                                                         const unsigned limit,
                                                         std::vector<std::int64_t>& path)
{
  return arc_longest_unobserved_run(mask, a, vec2pix(mask.base, a), b, arc_pixel_size(mask),
                                    limit, path);
}


// ============================================================================


void otswap::rejectMaskCrossings (Result& result, const Mask& mask,
                                  const unsigned maxUnobservedPixelsCrossed)
{
  const std::size_t nObjects = result.nObjects;
  const unsigned nRealizations = result.nRealizations;
  const std::size_t nDisplacements = (std::size_t)nRealizations * nObjects;
  const std::size_t expected = 3 * nDisplacements;

  if (result.geometry != Geometry::Lightcone)
    throw Error("rejectMaskCrossings filters a lightcone result, and this one has box geometry; "
                "a box has no sky to filter against");

  if (nObjects == 0 || nRealizations == 0)
    throw Error("the result is empty");

  if (result.displacement.size() != expected)
    throw Error("the result is malformed: displacement holds " +
                std::to_string(result.displacement.size()) + " entries for " +
                std::to_string(nRealizations) + " realizations of " + std::to_string(nObjects) +
                " objects");

  if (result.tracers.size() != 3 * nObjects)
    throw Error("the result is malformed: tracers holds " + std::to_string(result.tracers.size()) +
                " entries for " + std::to_string(nObjects) + " objects; the arcs start at the "
                "tracers, which a result returned by reconstructLightcone carries");

  if (result.valid.size() != nDisplacements)
    throw Error("the result is malformed: valid holds " + std::to_string(result.valid.size()) +
                " entries against " + std::to_string(nDisplacements) + " displacements");

  for (std::size_t k = 0; k < nDisplacements; ++k)
    if (result.valid[k] > 1)
      throw Error("the result is malformed: valid holds " + std::to_string(result.valid[k]) +
                  " at entry " + std::to_string(k) + " (from 0); every entry must be 0 or 1");

  internal::check_flags(result);

  for (std::size_t i = 0; i < nObjects; ++i) {
    if (internal::excluded(result, i)) continue;
    for (std::size_t c = 0; c < 3; ++c)
      if (!std::isfinite(result.tracers[3*i+c]))
        throw Error("the result is malformed: tracers holds a non-finite coordinate at entry " +
                    std::to_string(3*i+c) + " (from 0)");
  }

  for (std::size_t k = 0; k < nDisplacements; ++k) {
    const std::size_t i = k % nObjects;
    if (internal::excluded(result, i)) {
      if (result.valid[k] != 0)
        throw Error("the result is malformed: object " + std::to_string(i) + " (from 0) took no part in the "
                    "reconstruction (outside the redshift or distance cut, the mask or the box limits) but "
                    "has a valid displacement at entry " + std::to_string(k) + " (from 0)");
      continue;
    }
    if (result.valid[k] == 0) continue;
    for (std::size_t c = 0; c < 3; ++c)
      if (!std::isfinite(result.displacement[3*k+c]))
        throw Error("the result is malformed: displacement holds a non-finite coordinate at entry " +
                    std::to_string(3*k+c) + " (from 0)");
  }

  if (result.filteredNside != 0 && result.filteredNside != mask.nside())
    throw Error("this result was already filtered against a mask of NSIDE " +
                std::to_string(result.filteredNside) + ", and this mask has NSIDE " +
                std::to_string(mask.nside()) + "; a pixel count is not comparable "
                "between the two, so filter an unfiltered result instead");

  const internal::PixelMask& pixels = mask.m_impl->mask;
  const double pixelSize = internal::arc_pixel_size(pixels);
  const std::ptrdiff_t nValidBefore = std::count(result.valid.begin(), result.valid.end(), 1);

  std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t i = 0; i < nObjects; ++i) {
    try {

      std::vector<std::int64_t> path;

      const double ex = result.tracers[3*i];
      const double ey = result.tracers[3*i+1];
      const double ez = result.tracers[3*i+2];

      // Left out objects may hold any coordinates and have no valid
      // displacement: the Eulerian direction is taken at the first valid one
      bool started = false;
      internal::Vec3 a;
      std::int64_t pixA = 0;

      for (unsigned rec = 0; rec < nRealizations; ++rec) {

        const std::size_t d = (std::size_t)rec * nObjects + i;
        if (!result.valid[d]) continue;

        if (!started) {
          started = true;
          const double r0 = std::sqrt(ex*ex + ey*ey + ez*ez);
          if (r0 == 0.) break;
          a = internal::Vec3(ex/r0, ey/r0, ez/r0);
          pixA = internal::vec2pix(pixels.base, a);
        }

        const std::size_t at = 3 * d;

        const double lx = ex + result.displacement[at];
        const double ly = ey + result.displacement[at+1];
        const double lz = ez + result.displacement[at+2];

        const double r1 = std::sqrt(lx*lx + ly*ly + lz*lz);
        if (r1 == 0.) continue;

        const internal::Vec3 b(lx/r1, ly/r1, lz/r1);

        if ((a + b).squared_length() == 0. ||
            internal::arc_longest_unobserved_run(pixels, a, pixA, b, pixelSize,
                                                 maxUnobservedPixelsCrossed, path) >
              maxUnobservedPixelsCrossed)
          result.valid[d] = 0;
      }

    }
    catch (...) {
#pragma omp critical
      if (!failure) failure = std::current_exception();
    }
  }
  if (failure) std::rethrow_exception(failure);

  const std::ptrdiff_t nValidAfter = std::count(result.valid.begin(), result.valid.end(), 1);
  SelectionCounts& counts = result.selection;
  if (!counts.maxUnobservedPixelsCrossed) counts.displacements = (std::size_t)nValidBefore;
  counts.maxUnobservedPixelsCrossed = maxUnobservedPixelsCrossed;
  counts.displacementsCrossingMask += (std::size_t)(nValidBefore - nValidAfter);

  const bool mean = result.nRealizationsMean > 1 && nRealizations == 1 &&
                    result.validRealizations.size() == nObjects;
  std::vector<unsigned> averaged;
  if (mean) averaged = result.validRealizations;
  internal::summarize(result);
  if (mean)
    for (std::size_t i = 0; i < nObjects; ++i)
      if (result.validRealizations[i] > 0) result.validRealizations[i] = averaged[i];
  result.filteredNside = mask.nside();
}
