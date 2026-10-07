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
 *  @file src/RSD.cpp
 *
 *  @brief The redshift-space correction, declared in otswap/RSD.h.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <meshsearch/MeshGrid.h>

#include "detmath.h"
#include "internal.h"

namespace {

  // The gaussian is truncated at this many sigma.
  constexpr double kTruncationInSigma = 3.;

  // The grid cell, in mean separations of the positions averaged over.
  constexpr double kCellInSeparations = 4.;

  const double kNaN = std::numeric_limits<double>::quiet_NaN();

  // Number of objects in a flat array of 3*N entries that may hold NaN but
  // no infinity, such as a displacement.
  std::size_t check_nan_allowed (const std::vector<double>& a, const std::string& name)
  {
    if (a.size() % 3 != 0)
      throw otswap::Error(name + " holds " + std::to_string(a.size()) +
                          " entries, which is not a multiple of three");
    for (std::size_t i = 0; i < a.size(); ++i)
      if (std::isinf(a[i]))
        throw otswap::Error(name + " holds an infinite value at entry " + std::to_string(i) +
                            " (object " + std::to_string(i/3) + ")");
    return a.size() / 3;
  }

  // Positions that define a radial line of sight: finite and away from the
  // origin.
  std::size_t check_radial_positions (const std::vector<double>& positions)
  {
    const std::size_t n = otswap::internal::check_coordinates(positions, "the position array");
    for (std::size_t i = 0; i < n; ++i) {
      const double x = positions[3*i], y = positions[3*i+1], z = positions[3*i+2];
      const double r = std::sqrt(x*x + y*y + z*z);
      if (!(std::isfinite(r) && r > 0.))
        throw otswap::Error("object " + std::to_string(i) + " is at the origin, or too far from "
                            "it for its distance to be finite; its line of sight is undefined");
    }
    return n;
  }

  void check_axis (const unsigned axis)
  {
    if (axis > 2)
      throw otswap::Error("the line-of-sight axis is " + std::to_string(axis) +
                          "; it must be 0, 1 or 2");
  }

  void check_shift (const std::vector<double>& shift, const std::size_t n)
  {
    if (shift.size() != n)
      throw otswap::Error("the shift holds " + std::to_string(shift.size()) +
                          " entries for " + std::to_string(n) + " objects");
    for (std::size_t i = 0; i < n; ++i)
      if (std::isinf(shift[i]))
        throw otswap::Error("the shift of object " + std::to_string(i) + " is infinite");
  }

  // The fields of a Result the correction reads, checked against the
  // number of tracers given.
  void check_result (const otswap::Result& result, const std::size_t n)
  {
    if (result.nObjects != n)
      throw otswap::Error("the result describes " + std::to_string(result.nObjects) +
                          " objects and the tracer array " + std::to_string(n) +
                          "; pass the tracers the reconstruction was run on");
    if (result.meanDisplacement.size() != 3 * n || result.validRealizations.size() != n)
      throw otswap::Error("the result is malformed: meanDisplacement and validRealizations hold " +
                          std::to_string(result.meanDisplacement.size()) + " and " +
                          std::to_string(result.validRealizations.size()) + " entries for " +
                          std::to_string(n) + " objects");
    otswap::internal::check_flags(result);
    for (std::size_t i = 0; i < n; ++i) {
      if (otswap::internal::excluded(result, i) && result.validRealizations[i] != 0)
        throw otswap::Error("the result is malformed: object " + std::to_string(i) + " took no part "
                            "in the reconstruction (outside the redshift cut or the mask) but has "
                            "valid realizations");
      if (result.validRealizations[i] > 0)
        for (std::size_t c = 0; c < 3; ++c)
          if (!std::isfinite(result.meanDisplacement[3*i+c]))
            throw otswap::Error("the result is malformed: object " + std::to_string(i) +
                                " has valid realizations but a non-finite mean displacement");
    }
  }

  /// Indices, increasing, of the tracers the correction works on: all but
  /// those outside the result's redshift cut or mask.
  std::vector<std::size_t> kept_tracers (const otswap::Result& result, const std::size_t n)
  {
    std::vector<std::size_t> keep;
    keep.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
      if (!otswap::internal::excluded(result, i)) keep.push_back(i);
    return keep;
  }

  std::vector<double> rows_of (const std::vector<double>& a, const std::vector<std::size_t>& keep)
  {
    std::vector<double> out(3 * keep.size());
    for (std::size_t k = 0; k < keep.size(); ++k)
      for (std::size_t c = 0; c < 3; ++c) out[3*k+c] = a[3*keep[k]+c];
    return out;
  }

  // Bias at z: linear between the nodes, extrapolated from the end
  // segments; sets extrapolated when z lies outside the nodes.
  double bias_at (const std::vector<double>& zb, const std::vector<double>& b, const double z,
                  bool& extrapolated)
  {
    const std::size_t n = zb.size();
    const std::size_t k = (std::size_t)(std::lower_bound(zb.begin(), zb.end(), z) - zb.begin());
    extrapolated = false;

    if (k == 0) {
      if (z == zb[0]) return b[0];
      extrapolated = true;
      return b[0] + (b[1] - b[0]) * (z - zb[0]) / (zb[1] - zb[0]);
    }
    if (k == n) {
      extrapolated = true;
      return b[n-1] + (b[n-1] - b[n-2]) * (z - zb[n-1]) / (zb[n-1] - zb[n-2]);
    }
    if (z == zb[k]) return b[k];
    return b[k-1] + (b[k] - b[k-1]) * (z - zb[k-1]) / (zb[k] - zb[k-1]);
  }

  double growth_rate_of (const otswap::DistanceTable& distances, const double z,
                         const std::string& who)
  {
    try {
      return distances.growthRateAt(z);
    }
    catch (const otswap::Error& e) {
      throw otswap::Error(who + ": " + e.what());
    }
  }

  // f / (b + 3 f / 5) at z for object who; counts an extrapolation of b.
  double factor_at (const std::vector<double>& zb, const std::vector<double>& b, const double z,
                    const otswap::DistanceTable& distances, const std::string& who,
                    std::size_t& nExtrapolated)
  {
    if (!std::isfinite(z))
      throw otswap::Error("the redshift of " + who + " is not finite");

    bool extrapolated = false;
    const double bz = bias_at(zb, b, z, extrapolated);
    if (extrapolated) {
      ++nExtrapolated;
      if (!(bz > 0.))
        throw otswap::Error("the bias extrapolated to the redshift of " + who + ", " +
                            std::to_string(z) + ", is " + std::to_string(bz) +
                            "; it must be positive. The bias table covers [" +
                            std::to_string(zb.front()) + ", " + std::to_string(zb.back()) + "]");
    }

    const double f = growth_rate_of(distances, z, who);
    return f / (bz + 3. * f / 5.);
  }

  // An empty catalogue with n rows, every position NaN.
  otswap::RealSpaceCatalog empty_catalog (const std::size_t n)
  {
    otswap::RealSpaceCatalog catalog;
    catalog.nObjects = n;
    catalog.positions.assign(3 * n, kNaN);
    return catalog;
  }

}


// ============================================================================


std::vector<double> otswap::lineOfSightProjection (const std::vector<double>& positions,
                                                   const std::vector<double>& displacement)
{
  const std::size_t n = check_radial_positions(positions);
  if (check_nan_allowed(displacement, "the displacement array") != n || displacement.size() != positions.size())
    throw Error("the displacement array holds " + std::to_string(displacement.size()) +
                " entries for " + std::to_string(n) + " positions");

  std::vector<double> projection(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = positions[3*i], y = positions[3*i+1], z = positions[3*i+2];
    const double r = std::sqrt(x*x + y*y + z*z);
    projection[i] = (x*displacement[3*i] + y*displacement[3*i+1] + z*displacement[3*i+2]) / r;
  }
  return projection;
}


// ============================================================================


std::vector<double> otswap::lineOfSightProjection (const std::vector<double>& displacement,
                                                   const unsigned axis)
{
  check_axis(axis);
  const std::size_t n = check_nan_allowed(displacement, "the displacement array");

  std::vector<double> projection(n);
  for (std::size_t i = 0; i < n; ++i) projection[i] = displacement[3*i+axis];
  return projection;
}


// ============================================================================


double otswap::internal::neighbour_cell (const std::vector<double>& positions)
{
  const std::size_t n = positions.size() / 3;

  double side[3];
  for (int k = 0; k < 3; ++k) {
    double lo = positions[(std::size_t)k], hi = lo;
    for (std::size_t i = 1; i < n; ++i) {
      lo = std::min(lo, positions[3*i+k]);
      hi = std::max(hi, positions[3*i+k]);
    }
    side[k] = hi - lo;
  }

  const double volume = side[0] * side[1] * side[2];
  return kCellInSeparations * det_pow(volume / (double)n, 1./3.);
}


// ============================================================================


std::vector<double> otswap::internal::neighbour_average (const std::vector<double>& positions,
                                                         const std::vector<double>& values,
                                                         const std::vector<unsigned>& validRealizations,
                                                         const double sigma,
                                                         const bool weightByRealizations,
                                                         const double cellSize,
                                                         std::vector<unsigned>& nNeighbours,
                                                         std::vector<unsigned>& nRealizations)
{
  const std::size_t n = check_coordinates(positions, "the position array");

  if (values.size() != n)
    throw Error("values holds " + std::to_string(values.size()) + " entries for " +
                std::to_string(n) + " objects");
  if (validRealizations.size() != n)
    throw Error("validRealizations holds " + std::to_string(validRealizations.size()) +
                " entries for " + std::to_string(n) + " objects");
  if (!std::isfinite(sigma) || sigma < 0.)
    throw Error("sigma is " + std::to_string(sigma) + "; it must be finite and non-negative");

  for (std::size_t i = 0; i < n; ++i)
    if (validRealizations[i] > 0 && !std::isfinite(values[i]))
      throw Error("the value of object " + std::to_string(i) + " is not finite, and it has " +
                  std::to_string(validRealizations[i]) + " valid realizations");

  nNeighbours.assign(n, 0u);
  nRealizations.assign(n, 0u);

  if (sigma == 0.) {
    for (std::size_t i = 0; i < n; ++i)
      if (validRealizations[i] > 0) {
        nNeighbours[i] = 1;
        nRealizations[i] = validRealizations[i];
      }
    return values;
  }

  // The grid spans every position, so each query point lies inside it, and
  // holds the valid objects only; grid index g is object member[g].
  std::vector<std::vector<double>> limits(3, std::vector<double>(2));
  for (int k = 0; k < 3; ++k) {
    double lo = positions[(std::size_t)k], hi = lo;
    for (std::size_t i = 1; i < n; ++i) {
      lo = std::min(lo, positions[3*i+k]);
      hi = std::max(hi, positions[3*i+k]);
    }
    limits[(std::size_t)k] = {lo, hi};
  }

  std::vector<std::size_t> member;
  std::vector<double> X, Y, Z;
  for (std::size_t i = 0; i < n; ++i)
    if (validRealizations[i] > 0) {
      member.push_back(i);
      X.push_back(positions[3*i]);
      Y.push_back(positions[3*i+1]);
      Z.push_back(positions[3*i+2]);
    }

  std::vector<double> average(n, kNaN);
  if (member.empty()) return average;

  // Built before the parallel region: a cell size the grid refuses, as for
  // positions spanning no volume, raises meshsearch's error here.
  const meshsearch::MeshGrid grid(X, Y, Z, cellSize, limits);

  const double cutoff = kTruncationInSigma * sigma;
  const double inverse = 1. / (2. * sigma * sigma);

  // Every argument was checked above, so the search cannot fail but for
  // memory; whatever is thrown is kept and rethrown after the region.
  std::exception_ptr failure;

#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t i = 0; i < n; ++i) {
    try {
      const double x = positions[3*i], y = positions[3*i+1], z = positions[3*i+2];
      const std::vector<unsigned int> found = grid.closeObjects(x, y, z, cutoff);

      double sum = 0., norm = 0.;
      unsigned realizations = 0;
      for (const unsigned int g : found) {
        const std::size_t j = member[g];
        const double dx = positions[3*j] - x, dy = positions[3*j+1] - y, dz = positions[3*j+2] - z;
        double weight = weightByRealizations ? (double)validRealizations[j] : 1.;
        weight *= det_exp(-(dx*dx + dy*dy + dz*dz) * inverse);
        sum += weight * values[j];
        norm += weight;
        realizations += validRealizations[j];
      }

      if (!found.empty()) average[i] = sum / norm;
      nNeighbours[i] = (unsigned)found.size();
      nRealizations[i] = realizations;
    }
    catch (...) {
#pragma omp critical
      if (!failure) failure = std::current_exception();
    }
  }

  if (failure) std::rethrow_exception(failure);

  return average;
}


// ============================================================================


std::vector<double> otswap::neighbourAverage (const std::vector<double>& positions,
                                              const std::vector<double>& values,
                                              const std::vector<unsigned>& validRealizations,
                                              const double sigma,
                                              const bool weightByRealizations,
                                              std::vector<unsigned>& nNeighbours,
                                              std::vector<unsigned>& nRealizations)
{
  internal::check_coordinates(positions, "the position array");
  return internal::neighbour_average(positions, values, validRealizations, sigma,
                                     weightByRealizations, internal::neighbour_cell(positions),
                                     nNeighbours, nRealizations);
}


// ============================================================================


std::vector<double> otswap::neighbourAverage (const std::vector<double>& positions,
                                              const std::vector<double>& values,
                                              const std::vector<unsigned>& validRealizations,
                                              const double sigma,
                                              const bool weightByRealizations)
{
  std::vector<unsigned> nNeighbours, nRealizations;
  return neighbourAverage(positions, values, validRealizations, sigma, weightByRealizations,
                          nNeighbours, nRealizations);
}


// ============================================================================


void otswap::internal::check_bias_table (const std::vector<double>& redshift,
                                         const std::vector<double>& bias,
                                         const std::string& what,
                                         const std::string& entry,
                                         const std::size_t first)
{
  auto name = [&] (const std::size_t k) { return entry + " " + std::to_string(k + first); };

  if (redshift.size() != bias.size())
    throw Error(what + " holds " + std::to_string(redshift.size()) + " redshifts and " +
                std::to_string(bias.size()) + " bias values");
  if (redshift.size() < 2)
    throw Error(what + " holds " + std::to_string(redshift.size()) +
                " nodes; at least 2 are needed");

  for (std::size_t k = 0; k < redshift.size(); ++k) {
    if (!std::isfinite(redshift[k]))
      throw Error(what + ": the redshift of " + name(k) + " is not finite");
    if (k > 0 && !(redshift[k] > redshift[k-1]))
      throw Error(what + ": the redshift of " + name(k) + ", " +
                  std::to_string(redshift[k]) + ", is not greater than the one before, " +
                  std::to_string(redshift[k-1]) + "; the redshifts must increase strictly");
    if (!std::isfinite(bias[k]) || !(bias[k] > 0.))
      throw Error(what + ": the bias of " + name(k) + " is " +
                  std::to_string(bias[k]) + "; it must be finite and positive");
  }
}


// ============================================================================


std::vector<double> otswap::rsdFactor (const std::vector<double>& redshift,
                                       const DistanceTable& distances,
                                       const std::vector<double>& biasRedshift,
                                       const std::vector<double>& bias,
                                       std::size_t& nExtrapolated)
{
  internal::check_bias_table(biasRedshift, bias, "the bias table");
  if (!distances.hasGrowthRate())
    throw Error("the distance table carries no growth rate, which the correction needs");

  nExtrapolated = 0;
  std::vector<double> factor(redshift.size());

  for (std::size_t i = 0; i < redshift.size(); ++i)
    factor[i] = factor_at(biasRedshift, bias, redshift[i], distances,
                          "object " + std::to_string(i), nExtrapolated);

  return factor;
}


// ============================================================================


std::vector<double> otswap::rsdFactor (const std::vector<double>& redshift,
                                       const DistanceTable& distances,
                                       const std::vector<double>& biasRedshift,
                                       const std::vector<double>& bias)
{
  std::size_t nExtrapolated = 0;
  std::vector<double> factor = rsdFactor(redshift, distances, biasRedshift, bias, nExtrapolated);

  if (nExtrapolated > 0)
    std::clog << "otswap: b(z) extrapolated at " << nExtrapolated << " of " << redshift.size()
              << " redshifts, outside the bias table's range [" << biasRedshift.front() << ", "
              << biasRedshift.back() << "]" << std::endl;

  return factor;
}


// ============================================================================


double otswap::rsdFactorBox (const double redshift, const DistanceTable& distances,
                             const double bias)
{
  if (!std::isfinite(redshift))
    throw Error("the redshift is not finite");
  if (!std::isfinite(bias) || !(bias > 0.))
    throw Error("the bias is " + std::to_string(bias) + "; it must be finite and positive");
  if (!distances.hasGrowthRate())
    throw Error("the distance table carries no growth rate, which the correction needs");

  const double f = growth_rate_of(distances, redshift, "the box redshift");
  return f / (bias + 3. * f / 5.);
}


// ============================================================================


std::vector<double> otswap::shiftAlongLineOfSight (const std::vector<double>& positions,
                                                   const std::vector<double>& shift)
{
  const std::size_t n = check_radial_positions(positions);
  check_shift(shift, n);

  std::vector<double> shifted(3 * n);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = positions[3*i], y = positions[3*i+1], z = positions[3*i+2];
    const double scale = shift[i] / std::sqrt(x*x + y*y + z*z);
    shifted[3*i]   = x + scale * x;
    shifted[3*i+1] = y + scale * y;
    shifted[3*i+2] = z + scale * z;
  }
  return shifted;
}


// ============================================================================


std::vector<double> otswap::shiftAlongLineOfSight (const std::vector<double>& positions,
                                                   const std::vector<double>& shift,
                                                   const unsigned axis)
{
  check_axis(axis);
  const std::size_t n = internal::check_coordinates(positions, "the position array");
  check_shift(shift, n);

  std::vector<double> shifted = positions;
  for (std::size_t i = 0; i < n; ++i) shifted[3*i+axis] += shift[i];
  return shifted;
}


// ============================================================================


otswap::RealSpaceCatalog otswap::realSpaceLightcone (const Result& result,
                                                     const std::vector<double>& tracersSky,
                                                     const DistanceTable& distances,
                                                     const std::vector<double>& biasRedshift,
                                                     const std::vector<double>& bias,
                                                     const double sigma,
                                                     const bool weightByRealizations,
                                                     std::size_t& nExtrapolated)
{
  const std::size_t n = internal::check_sky(tracersSky, "the tracer sky array");
  check_result(result, n);
  internal::check_bias_table(biasRedshift, bias, "the bias table");
  if (!distances.hasGrowthRate())
    throw Error("the distance table carries no growth rate, which the correction needs");
  if (!std::isfinite(sigma) || sigma < 0.)
    throw Error("sigma is " + std::to_string(sigma) + "; it must be finite and non-negative");

  const std::vector<std::size_t> keep = kept_tracers(result, n);
  const std::size_t m = keep.size();

  RealSpaceCatalog catalog = empty_catalog(n);
  catalog.nNeighbours.assign(n, 0u);
  catalog.nRealizationsAveraged.assign(n, 0u);
  nExtrapolated = 0;
  if (m == 0) {
    for (std::size_t i = 0; i < n; ++i) catalog.uncorrected.push_back(i);
    return catalog;
  }

  std::vector<double> positions(3 * m), redshift(m), distance(m), factor(m);
  for (std::size_t k = 0; k < m; ++k) {
    const std::size_t i = keep[k];
    const std::string who = "object " + std::to_string(i);
    redshift[k] = tracersSky[3*i+2];
    try {
      distance[k] = distances.distanceAt(redshift[k]);
    }
    catch (const Error& e) {
      throw Error(who + " of the tracer sky array: " + e.what());
    }
    if (!(distance[k] > 0.))
      throw Error(who + " is at zero distance; its line of sight is undefined");
    internal::to_cartesian(tracersSky[3*i], tracersSky[3*i+1], distance[k],
                           positions[3*k], positions[3*k+1], positions[3*k+2]);
    factor[k] = factor_at(biasRedshift, bias, redshift[k], distances, who, nExtrapolated);
  }

  const std::vector<double> projection =
    lineOfSightProjection(positions, rows_of(result.meanDisplacement, keep));

  std::vector<unsigned> valid(m), nNeighbours, nRealizations;
  for (std::size_t k = 0; k < m; ++k) valid[k] = result.validRealizations[keep[k]];
  const std::vector<double> averaged =
    neighbourAverage(positions, projection, valid, sigma, weightByRealizations,
                     nNeighbours, nRealizations);

  const double dMin = distances.distanceAt(distances.minRedshift());
  const double dMax = distances.distanceAt(distances.maxRedshift());

  std::size_t next = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (next == m || keep[next] != i) {
      catalog.uncorrected.push_back(i);
      continue;
    }
    const std::size_t k = next++;
    catalog.nNeighbours[i] = nNeighbours[k];
    catalog.nRealizationsAveraged[i] = nRealizations[k];

    const double shift = factor[k] * averaged[k];
    if (!std::isfinite(shift)) {
      catalog.uncorrected.push_back(i);
      continue;
    }

    const double corrected = distance[k] + shift;
    double z = 0.;
    bool inside = corrected > 0.;
    if (inside) {
      try {
        z = distances.redshiftAt(corrected);
      }
      catch (const Error&) {
        inside = false;
      }
    }
    if (!inside)
      throw Error("object " + std::to_string(i) + ": its corrected comoving distance, " +
                  std::to_string(corrected) + " Mpc/h, " +
                  (corrected > 0. ? "lies outside" : "is not positive, and so outside") +
                  " the distance table's range [" + std::to_string(dMin) + ", " +
                  std::to_string(dMax) + "]; use a distance table covering a wider redshift "
                  "range than the catalogue's, [" + std::to_string(distances.minRedshift()) +
                  ", " + std::to_string(distances.maxRedshift()) + "] here");

    catalog.positions[3*i]   = tracersSky[3*i];
    catalog.positions[3*i+1] = tracersSky[3*i+1];
    catalog.positions[3*i+2] = z;
  }

  return catalog;
}


// ============================================================================


otswap::RealSpaceCatalog otswap::realSpaceLightcone (const Result& result,
                                                     const std::vector<double>& tracersSky,
                                                     const DistanceTable& distances,
                                                     const std::vector<double>& biasRedshift,
                                                     const std::vector<double>& bias,
                                                     const double sigma,
                                                     const bool weightByRealizations)
{
  std::size_t nExtrapolated = 0;
  RealSpaceCatalog catalog = realSpaceLightcone(result, tracersSky, distances, biasRedshift, bias,
                                                sigma, weightByRealizations, nExtrapolated);

  if (nExtrapolated > 0)
    std::clog << "otswap: b(z) extrapolated at " << nExtrapolated << " of "
              << kept_tracers(result, catalog.nObjects).size()
              << " redshifts, outside the bias table's range [" << biasRedshift.front() << ", "
              << biasRedshift.back() << "]" << std::endl;

  return catalog;
}


// ============================================================================


otswap::RealSpaceCatalog otswap::realSpaceBox (const Result& result,
                                               const std::vector<double>& tracers,
                                               const unsigned axis,
                                               const double redshift,
                                               const DistanceTable& distances,
                                               const double bias,
                                               const double sigma,
                                               const bool weightByRealizations)
{
  const std::size_t n = internal::check_coordinates(tracers, "the tracer array");
  check_axis(axis);
  check_result(result, n);
  const double factor = rsdFactorBox(redshift, distances, bias);

  const std::vector<std::size_t> keep = kept_tracers(result, n);
  const std::size_t m = keep.size();

  RealSpaceCatalog catalog = empty_catalog(n);
  catalog.nNeighbours.assign(n, 0u);
  catalog.nRealizationsAveraged.assign(n, 0u);

  std::vector<double> averaged;
  std::vector<unsigned> nNeighbours, nRealizations;
  if (m > 0) {
    std::vector<unsigned> valid(m);
    for (std::size_t k = 0; k < m; ++k) valid[k] = result.validRealizations[keep[k]];
    const std::vector<double> projection =
      lineOfSightProjection(rows_of(result.meanDisplacement, keep), axis);
    averaged = neighbourAverage(rows_of(tracers, keep), projection, valid, sigma,
                                weightByRealizations, nNeighbours, nRealizations);
  }

  std::size_t next = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (next == m || keep[next] != i) {
      catalog.uncorrected.push_back(i);
      continue;
    }
    const std::size_t k = next++;
    catalog.nNeighbours[i] = nNeighbours[k];
    catalog.nRealizationsAveraged[i] = nRealizations[k];

    const double shift = factor * averaged[k];
    if (!std::isfinite(shift)) {
      catalog.uncorrected.push_back(i);
      continue;
    }
    for (unsigned c = 0; c < 3; ++c) catalog.positions[3*i+c] = tracers[3*i+c];
    catalog.positions[3*i+axis] += shift;
  }

  return catalog;
}
