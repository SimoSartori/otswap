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
 *  @file src/Selection.cpp
 *
 *  @brief The selections of a reconstruction, by redshift cut, mask or box
 *  limits, and the checks and expansion that follow them.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "internal.h"

namespace {

  const double kNaN = std::numeric_limits<double>::quiet_NaN();

  const char* const kAxes[] = {"x", "y", "z"};

}


// ============================================================================


otswap::internal::Selection otswap::internal::select (const std::vector<double>& sky,
                                                      const RedshiftCut& cut, const Mask* mask)
{
  const std::size_t n = sky.size() / 3;
  Selection s;
  s.outsideCut.assign(n, 0);
  s.outsideMask.assign(n, 0);
  s.outsideBox.assign(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const bool outsideCut = !(sky[3*i+2] >= cut.min && sky[3*i+2] <= cut.max);
    const bool outsideMask = mask != nullptr && !mask->allows(sky[3*i], sky[3*i+1]);
    s.outsideCut[i] = outsideCut ? 1 : 0;
    s.outsideMask[i] = outsideMask ? 1 : 0;
    s.nOutsideCut += outsideCut ? 1 : 0;
    s.nOutsideMask += outsideMask ? 1 : 0;
    s.nOutsideBoth += (outsideCut && outsideMask) ? 1 : 0;
    if (!outsideCut && !outsideMask) s.keep.push_back(i);
  }
  return s;
}


// ============================================================================


otswap::internal::Selection otswap::internal::select (const std::vector<double>& positions,
                                                      const BoxLimits& limits)
{
  const std::size_t n = positions.size() / 3;
  Selection s;
  s.outsideCut.assign(n, 0);
  s.outsideMask.assign(n, 0);
  s.outsideBox.assign(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    bool outside = false;
    for (std::size_t k = 0; k < 3; ++k)
      outside = outside || !(positions[3*i+k] >= limits.lower[k] && positions[3*i+k] <= limits.upper[k]);
    s.outsideBox[i] = outside ? 1 : 0;
    s.nOutsideBox += outside ? 1 : 0;
    if (!outside) s.keep.push_back(i);
  }
  return s;
}


// ============================================================================


void otswap::internal::check_box (const BoxLimits& limits)
{
  for (std::size_t k = 0; k < 3; ++k)
    for (const double bound : {limits.lower[k], limits.upper[k]})
      if (!std::isfinite(bound))
        throw Error("the box limits hold a bound of " + std::to_string(bound) + " on axis " +
                    kAxes[k] + "; every bound must be finite");

  for (std::size_t k = 0; k < 3; ++k)
    if (!(limits.lower[k] < limits.upper[k]))
      throw Error("the box limits are [" + std::to_string(limits.lower[k]) + ", " +
                  std::to_string(limits.upper[k]) + "] on axis " + kAxes[k] +
                  "; the lower bound must be below the upper one");
}


// ============================================================================


std::vector<double> otswap::internal::rows (const std::vector<double>& a,
                                            const std::vector<std::size_t>& keep)
{
  std::vector<double> out(3 * keep.size());
  for (std::size_t k = 0; k < keep.size(); ++k)
    for (std::size_t c = 0; c < 3; ++c) out[3*k+c] = a[3*keep[k]+c];
  return out;
}


// ============================================================================


std::string otswap::internal::selection_name (const SelectionCounts& counts)
{
  if (counts.boxLimits) {
    const BoxLimits& b = *counts.boxLimits;
    std::string name = "the box ";
    for (std::size_t k = 0; k < 3; ++k)
      name += (k ? " x [" : "[") + std::to_string(b.lower[k]) + ", " + std::to_string(b.upper[k]) + "]";
    return name;
  }
  const std::string cut = !counts.redshiftCut ? std::string()
    : "the redshift cut [" + std::to_string(counts.redshiftCut->min) + ", " +
      std::to_string(counts.redshiftCut->max) + "]";
  return !counts.maskApplied ? cut : counts.redshiftCut ? cut + " and the mask" : "the mask";
}


// ============================================================================


void otswap::internal::check_kept_tracers (const SelectionCounts& counts,
                                           const std::size_t keptTracers)
{
  const std::size_t nTracers = counts.tracers;
  if (keptTracers < min_objects())
    throw Error(selection_name(counts) + " keeps " + std::to_string(keptTracers) + " of " +
                std::to_string(nTracers) + " tracers (" + std::to_string(nTracers - keptTracers) +
                " dropped); at least " + std::to_string(min_objects()) + " are required");
}


// ============================================================================


void otswap::internal::check_kept (const SelectionCounts& counts, const std::size_t keptTracers,
                                   const std::size_t keptRandoms, const unsigned nRealizations)
{
  check_kept_tracers(counts, keptTracers);
  const std::size_t nRandoms = counts.randoms;
  const std::size_t needed = (std::size_t)nRealizations * keptTracers;
  if (keptRandoms < needed)
    throw Error(selection_name(counts) + " keeps " + std::to_string(keptRandoms) + " of " +
                std::to_string(nRandoms) + " randoms (" + std::to_string(nRandoms - keptRandoms) +
                " dropped), for " + std::to_string(keptTracers) + " tracers kept; " +
                std::to_string(nRealizations) + " realizations need " + std::to_string(needed));
}


// ============================================================================


otswap::SelectionCounts otswap::internal::tracer_counts (const RedshiftCut& cut, const Mask* mask,
                                                         const std::size_t nObjects,
                                                         const Selection& selection)
{
  SelectionCounts counts;
  if (std::isfinite(cut.min) || std::isfinite(cut.max)) counts.redshiftCut = cut;
  counts.maskApplied = mask != nullptr;
  counts.tracers = nObjects;
  counts.tracersOutsideRedshiftCut = selection.nOutsideCut;
  counts.tracersOutsideMask = selection.nOutsideMask;
  counts.tracersOutsideBoth = selection.nOutsideBoth;
  return counts;
}


// ============================================================================


otswap::Result otswap::internal::expand (const Result& kept, const Selection& selection,
                                         const std::size_t nObjects)
{
  const std::vector<std::size_t>& keep = selection.keep;
  const unsigned nRealizations = kept.nRealizations;

  Result full;
  full.nObjects = nObjects;
  full.nRealizations = nRealizations;
  full.displacement.assign(3 * (std::size_t)nRealizations * nObjects, kNaN);
  full.valid.assign((std::size_t)nRealizations * nObjects, 0);
  full.outsideRedshiftCut = selection.outsideCut;
  full.outsideMask = selection.outsideMask;
  full.outsideBox = selection.outsideBox;
  full.filteredNside = kept.filteredNside;
  full.config = kept.config;

  for (std::size_t k = 0; k < keep.size(); ++k) {
    const std::size_t i = keep[k];
    for (unsigned rec = 0; rec < nRealizations; ++rec) {
      const std::size_t from = (std::size_t)rec * keep.size() + k;
      const std::size_t to = (std::size_t)rec * nObjects + i;
      full.valid[to] = kept.valid[from];
      for (std::size_t c = 0; c < 3; ++c)
        full.displacement[3*to+c] = kept.displacement[3*from+c];
    }
  }

  return full;
}
