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
#include <chrono>
#include <cmath>
#include <exception>
#include <iostream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

#include <meshsearch/MeshGrid.h>

#include "detmath.h"
#include "internal.h"

namespace {

  /// The gaussian is truncated at this many sigma.
  constexpr double kTruncationInSigma = 3.;

  /// The grid cell, in mean separations of the positions averaged over.
  constexpr double kCellInSeparations = 4.;

  const double kNaN = std::numeric_limits<double>::quiet_NaN();

  /// Number of objects in a flat array of 3*N entries that may hold NaN but
  /// no infinity, such as a displacement.
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

  /// Positions that define a radial line of sight: finite and away from the
  /// origin.
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

  /// The fields of a Result the correction reads: its geometry, its
  /// tracers (and their sky coordinates in a lightcone), the mean
  /// displacements and the counts of valid realizations.
  void check_result (const otswap::Result& result, const otswap::Geometry geometry,
                     const std::string& function)
  {
    const std::size_t n = result.nObjects;
    const bool lightcone = geometry == otswap::Geometry::Lightcone;
    if (result.geometry != geometry)
      throw otswap::Error(function + " corrects a " + (lightcone ? "lightcone" : "box") +
                          " result, and this one was reconstructed in " +
                          (lightcone ? "box" : "lightcone") + " geometry");
    if (result.tracers.size() != 3 * n || (lightcone && result.tracersSky.size() != 3 * n))
      throw otswap::Error("the result does not carry its tracers: tracers" +
                          std::string(lightcone ? " and tracersSky hold " : " holds ") +
                          std::to_string(result.tracers.size()) +
                          (lightcone ? " and " + std::to_string(result.tracersSky.size()) : "") +
                          " entries for " + std::to_string(n) + " objects; a result returned by " +
                          (lightcone ? "reconstructLightcone" : "reconstructBox") + " carries them");
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

  /// Bias at z: linear between the nodes, extrapolated from the end
  /// segments; sets extrapolated when z lies outside the nodes.
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

  /// f / (b + 3 f / 5) at z for object who; counts an extrapolation of b.
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

  /// A catalogue of the objects of result, nothing moved: every position,
  /// shift and factor NaN, every diagnostic 0, the status LeftOut for the
  /// tracers left out of the reconstruction and NoValidNeighbour for the
  /// others.
  otswap::RealSpaceCatalog empty_catalog (const otswap::Result& result)
  {
    const std::size_t n = result.nObjects;
    otswap::RealSpaceCatalog catalog;
    catalog.nObjects = n;
    catalog.geometry = result.geometry;
    if (result.geometry == otswap::Geometry::Lightcone) catalog.sky.assign(3 * n, kNaN);
    catalog.cartesian.assign(3 * n, kNaN);
    catalog.shift.assign(n, kNaN);
    catalog.factor.assign(n, kNaN);
    catalog.status.assign(n, otswap::CorrectionStatus::NoValidNeighbour);
    for (std::size_t i = 0; i < n; ++i)
      if (otswap::internal::excluded(result, i)) catalog.status[i] = otswap::CorrectionStatus::LeftOut;
    catalog.validRealizations = result.validRealizations;
    catalog.nNeighbours.assign(n, 0u);
    catalog.nRealizationsAveraged.assign(n, 0u);
    return catalog;
  }

  /// Write the report of a correction to config.log: first, under Normal,
  /// the lines given (the extrapolation of b(z)); under Detailed the
  /// breakdown of the corrected and of the uncorrected; then the line of
  /// the call.
  void report (const std::string& name, const otswap::RealSpaceCatalog& catalog,
               const std::string& first, const otswap::CorrectionConfig& config)
  {
    using otswap::CorrectionStatus;
    using otswap::Verbosity;
    if (config.log == nullptr || config.verbosity == Verbosity::Silent) return;

    std::size_t count[4] = {0, 0, 0, 0};
    for (const CorrectionStatus s : catalog.status) ++count[(std::size_t)s];
    const std::size_t corrected = count[(std::size_t)CorrectionStatus::Corrected] +
                                  count[(std::size_t)CorrectionStatus::MovedByNeighbours];
    const std::size_t uncorrected = catalog.nObjects - corrected;

    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << first;
    if (config.verbosity == Verbosity::Detailed) {
      out << "otswap: corrected " << corrected << " tracers, "
          << count[(std::size_t)CorrectionStatus::MovedByNeighbours]
          << " of them moved with the average of their neighbours (no valid realization)\n";
      out << "otswap: left " << uncorrected << " uncorrected: "
          << count[(std::size_t)CorrectionStatus::LeftOut] << " left out of the reconstruction, "
          << count[(std::size_t)CorrectionStatus::NoValidNeighbour]
          << " with no valid tracer within 3 sigma\n";
    }
    out << "otswap: " << name << ": corrected " << corrected << " of " << catalog.nObjects
        << " tracers in " << std::fixed << std::setprecision(2) << catalog.elapsedSeconds << " s; "
        << uncorrected << " left uncorrected\n";
    *config.log << out.str() << std::flush;
  }

}


// ============================================================================


std::vector<double> otswap::radialProjection (const std::vector<double>& displacement,
                                              const std::vector<double>& positions)
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


std::vector<double> otswap::axisProjection (const std::vector<double>& displacement,
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
  return kCellInSeparations * bounding_box_separation(positions);
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

  const meshsearch::MeshGrid grid(X, Y, Z, cellSize, limits);

  const double cutoff = kTruncationInSigma * sigma;
  const double inverse = 1. / (2. * sigma * sigma);

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


otswap::NeighbourAverage otswap::neighbourAverage (const std::vector<double>& positions,
                                                   const std::vector<double>& values,
                                                   const std::vector<unsigned>& validRealizations,
                                                   const double sigma,
                                                   const bool weightByRealizations)
{
  internal::check_coordinates(positions, "the position array");
  NeighbourAverage average;
  average.values = internal::neighbour_average(positions, values, validRealizations, sigma,
                                               weightByRealizations, internal::neighbour_cell(positions),
                                               average.nNeighbours, average.nRealizationsAveraged);
  return average;
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
                                       const BiasTable& bias,
                                       std::size_t* nExtrapolated)
{
  internal::check_bias_table(bias.redshift, bias.bias, "the bias table");
  if (!distances.hasGrowthRate())
    throw Error("the distance table carries no growth rate, which the correction needs");

  std::size_t extrapolated = 0;
  std::vector<double> factor(redshift.size());

  for (std::size_t i = 0; i < redshift.size(); ++i)
    factor[i] = factor_at(bias.redshift, bias.bias, redshift[i], distances,
                          "object " + std::to_string(i), extrapolated);

  if (nExtrapolated != nullptr) *nExtrapolated = extrapolated;
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


std::vector<double> otswap::shiftRadially (const std::vector<double>& positions,
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


std::vector<double> otswap::shiftAlongAxis (const std::vector<double>& positions,
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
                                                     const DistanceTable& distances,
                                                     const BiasTable& bias,
                                                     const double sigma,
                                                     const CorrectionConfig& config)
{
  const auto start = std::chrono::steady_clock::now();
  check_result(result, Geometry::Lightcone, "realSpaceLightcone");
  internal::check_bias_table(bias.redshift, bias.bias, "the bias table");
  if (!distances.hasGrowthRate())
    throw Error("the distance table carries no growth rate, which the correction needs");
  if (!std::isfinite(sigma) || sigma < 0.)
    throw Error("sigma is " + std::to_string(sigma) + "; it must be finite and non-negative");

  const std::size_t n = result.nObjects;
  const std::vector<std::size_t> keep = kept_tracers(result, n);
  const std::size_t m = keep.size();

  RealSpaceCatalog catalog = empty_catalog(result);
  catalog.sigma = sigma;
  catalog.weightByRealizations = config.weightByRealizations;
  std::size_t nExtrapolated = 0;
  if (m == 0) {
    catalog.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    report("realSpaceLightcone", catalog, "", config);
    return catalog;
  }

  const std::vector<double>& sky = result.tracersSky;
  internal::check_sky(rows_of(sky, keep), "the tracer sky array of the result");

  const std::vector<double> positions = rows_of(result.tracers, keep);
  std::vector<double> distance(m), factor(m);
  for (std::size_t k = 0; k < m; ++k) {
    const std::size_t i = keep[k];
    const std::string who = "object " + std::to_string(i);
    try {
      distance[k] = distances.distanceAt(sky[3*i+2]);
    }
    catch (const Error& e) {
      throw Error(who + " of the tracer sky array of the result: " + e.what());
    }
    if (!(distance[k] > 0.))
      throw Error(who + " is at zero distance; its line of sight is undefined");
    factor[k] = factor_at(bias.redshift, bias.bias, sky[3*i+2], distances, who, nExtrapolated);
  }

  const std::vector<double> projection =
    radialProjection(rows_of(result.meanDisplacement, keep), positions);

  std::vector<unsigned> valid(m);
  for (std::size_t k = 0; k < m; ++k) valid[k] = result.validRealizations[keep[k]];
  const NeighbourAverage averaged =
    neighbourAverage(positions, projection, valid, sigma, config.weightByRealizations);

  const double dMin = distances.distanceAt(distances.minRedshift());
  const double dMax = distances.distanceAt(distances.maxRedshift());

  for (std::size_t k = 0; k < m; ++k) {
    const std::size_t i = keep[k];
    catalog.nNeighbours[i] = averaged.nNeighbours[k];
    catalog.nRealizationsAveraged[i] = averaged.nRealizationsAveraged[k];
    catalog.factor[i] = factor[k];

    const double shift = factor[k] * averaged.values[k];
    if (!std::isfinite(shift)) continue;

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

    catalog.shift[i] = shift;
    catalog.status[i] = valid[k] > 0 ? CorrectionStatus::Corrected : CorrectionStatus::MovedByNeighbours;
    catalog.sky[3*i]   = internal::fold_ra(sky[3*i]);
    catalog.sky[3*i+1] = sky[3*i+1];
    catalog.sky[3*i+2] = z;

    const double x = positions[3*k], y = positions[3*k+1], w = positions[3*k+2];
    const double scale = shift / std::sqrt(x*x + y*y + w*w);
    catalog.cartesian[3*i]   = x + scale * x;
    catalog.cartesian[3*i+1] = y + scale * y;
    catalog.cartesian[3*i+2] = w + scale * w;
  }

  catalog.nExtrapolated = nExtrapolated;
  catalog.elapsedSeconds =
    std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  std::ostringstream extrapolation;
  extrapolation.imbue(std::locale::classic());
  if (nExtrapolated > 0)
    extrapolation << "otswap: b(z) extrapolated at " << nExtrapolated << " of " << m
                  << " redshifts, outside the bias table's range [" << bias.redshift.front() << ", "
                  << bias.redshift.back() << "]\n";
  report("realSpaceLightcone", catalog, extrapolation.str(), config);

  return catalog;
}


// ============================================================================


otswap::RealSpaceCatalog otswap::realSpaceBox (const Result& result,
                                               const unsigned axis,
                                               const double redshift,
                                               const DistanceTable& distances,
                                               const double bias,
                                               const double sigma,
                                               const CorrectionConfig& config)
{
  const auto start = std::chrono::steady_clock::now();
  check_axis(axis);
  check_result(result, Geometry::Box, "realSpaceBox");
  const double factor = rsdFactorBox(redshift, distances, bias);

  const std::size_t n = result.nObjects;
  const std::vector<std::size_t> keep = kept_tracers(result, n);
  const std::size_t m = keep.size();

  RealSpaceCatalog catalog = empty_catalog(result);
  catalog.sigma = sigma;
  catalog.weightByRealizations = config.weightByRealizations;
  catalog.axis = axis;
  catalog.boxRedshift = redshift;
  catalog.boxBias = bias;

  if (m > 0) {
    const std::vector<double> positions = rows_of(result.tracers, keep);
    std::vector<unsigned> valid(m);
    for (std::size_t k = 0; k < m; ++k) valid[k] = result.validRealizations[keep[k]];
    const std::vector<double> projection = axisProjection(rows_of(result.meanDisplacement, keep), axis);
    const NeighbourAverage averaged =
      neighbourAverage(positions, projection, valid, sigma, config.weightByRealizations);

    for (std::size_t k = 0; k < m; ++k) {
      const std::size_t i = keep[k];
      catalog.nNeighbours[i] = averaged.nNeighbours[k];
      catalog.nRealizationsAveraged[i] = averaged.nRealizationsAveraged[k];
      catalog.factor[i] = factor;

      const double shift = factor * averaged.values[k];
      if (!std::isfinite(shift)) continue;
      catalog.shift[i] = shift;
      catalog.status[i] = valid[k] > 0 ? CorrectionStatus::Corrected : CorrectionStatus::MovedByNeighbours;
      for (unsigned c = 0; c < 3; ++c) catalog.cartesian[3*i+c] = positions[3*k+c];
      catalog.cartesian[3*i+axis] += shift;
    }
  }

  catalog.elapsedSeconds =
    std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  report("realSpaceBox", catalog, "", config);
  return catalog;
}
