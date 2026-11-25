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
 *  @file src/Reconstruct.cpp
 *
 *  @brief The reconstruction core, shared by both geometries.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <limits>
#include <locale>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <meshsearch/MeshGrid.h>

#include "otswap/Random.h"

#include "detmath.h"
#include "internal.h"

namespace {


  /// Radius, in units of the local mean particle separation, of the ball
  /// around a starter from which the seeding pass draws tracers to pair.
  constexpr double kSeedingRadiusInMps = 8.;

  /// Radius, in units of the local mean particle separation, of the
  /// neighbourhood the three swap partners are drawn from.
  constexpr double kNeighbourhoodRadiusInMps = 3.;

  /// Neighbours held per tracer: the count a sphere of
  /// kNeighbourhoodRadiusInMps mean separations holds at mean density,
  /// (4/3) * pi * 3^3 = 113.097..., truncated to 113. It is also the size
  /// of the list the swap partners are drawn from, and one less than the
  /// smallest catalog accepted.
  constexpr unsigned kNeighbours =
    (unsigned)((4. / 3.) * 3.14159265358979323846 *
               kNeighbourhoodRadiusInMps * kNeighbourhoodRadiusInMps *
               kNeighbourhoodRadiusInMps);

  /// Pairs the seeding pass makes around one starter before moving on.
  constexpr unsigned kMaxPairsPerStarter = 31;

  /// The 24 permutations of four elements. The order is load-bearing: the
  /// loop keeps the first permutation that strictly improves on the best so
  /// far, so a different order can pick a different permutation of equal
  /// cost.
  constexpr int kPermutations[24][4] = {
    {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1},
    {1,0,2,3}, {1,0,3,2}, {1,2,0,3}, {1,2,3,0}, {1,3,0,2}, {1,3,2,0},
    {2,0,1,3}, {2,0,3,1}, {2,1,0,3}, {2,1,3,0}, {2,3,0,1}, {2,3,1,0},
    {3,0,1,2}, {3,0,2,1}, {3,1,0,2}, {3,1,2,0}, {3,2,0,1}, {3,2,1,0}
  };

  /// The periods of a periodic box and their halves, per axis.
  struct Periods {
    double period[3];
    double half[3];
  };

  /// The squared distance from a tracer to a random; in a periodic box,
  /// of the minimum image, each difference of at most one period in
  /// magnitude brought into [-L/2, L/2).
  template <bool Periodic>
  inline double squared_distance (const double tx, const double ty, const double tz,
                                  const double rx, const double ry, const double rz,
                                  const Periods& periods)
  {
    double dx = tx - rx;
    double dy = ty - ry;
    double dz = tz - rz;
    if constexpr (Periodic) {
      double* d[3] = {&dx, &dy, &dz};
      for (int k = 0; k < 3; ++k) {
        if (*d[k] >= periods.half[k]) *d[k] -= periods.period[k];
        else if (*d[k] < -periods.half[k]) *d[k] += periods.period[k];
      }
    }
    return dx*dx + dy*dy + dz*dz;
  }

  /// The local search of one realization, from the first guess in
  /// partner: sweeps of 4-opt moves over every tracer until the fraction
  /// of successful moves no longer exceeds the convergence threshold.
  /// @return the mean squared cost after each sweep, when @p trace is set.
  template <bool Periodic>
  __attribute__((always_inline)) inline void swap_loop (const std::vector<double>& tracerX, const std::vector<double>& tracerY,
                  const std::vector<double>& tracerZ, const std::vector<double>& randomX,
                  const std::vector<double>& randomY, const std::vector<double>& randomZ,
                  const std::vector<unsigned>& neighbours, const Periods& periods,
                  const double convergence, std::mt19937 engine, std::vector<unsigned>& partner,
                  std::vector<double>* trace)
  {
    const std::size_t nObjects = tracerX.size();
    const double* const tX = tracerX.data();
    const double* const tY = tracerY.data();
    const double* const tZ = tracerZ.data();
    const double* const rX = randomX.data();
    const double* const rY = randomY.data();
    const double* const rZ = randomZ.data();
    const unsigned* const near = neighbours.data();
    unsigned* const match = partner.data();
    const Periods p = periods;

    std::vector<double> costs(nObjects);
    double* const cost = costs.data();
    for (std::size_t i = 0; i < nObjects; ++i)
      cost[i] = squared_distance<Periodic>(tX[i], tY[i], tZ[i], rX[match[i]], rY[match[i]],
                                           rZ[match[i]], p);

    std::vector<unsigned> order(nObjects);
    std::iota(order.begin(), order.end(), 0u);

    unsigned H[4], R[4];
    double ratio = 1.;

    while (ratio > convergence) {

      otswap::shuffle(order, engine);
      unsigned changed = 0;

      for (const unsigned i : order) {
        H[0] = i;
        unsigned r1 = otswap::uniformInt(engine, kNeighbours);
        unsigned r2 = otswap::uniformInt(engine, kNeighbours);
        unsigned r3 = otswap::uniformInt(engine, kNeighbours);
        while (r1 == r2) r2 = otswap::uniformInt(engine, kNeighbours);
        while (r3 == r1 || r3 == r2) r3 = otswap::uniformInt(engine, kNeighbours);

        const unsigned* const list = near + (std::size_t)i * kNeighbours;
        H[1] = list[r1];
        H[2] = list[r2];
        H[3] = list[r3];

        for (int k = 0; k < 4; ++k) R[k] = match[H[k]];

        double matrix[4][4];
        for (int k = 0; k < 4; ++k) {
          const double tx = tX[H[k]], ty = tY[H[k]], tz = tZ[H[k]];
          for (int m = 0; m < 4; ++m)
            matrix[k][m] = squared_distance<Periodic>(tx, ty, tz, rX[R[m]], rY[R[m]], rZ[R[m]], p);
        }

        const double current = cost[H[0]] + cost[H[1]] + cost[H[2]] + cost[H[3]];
        double best = current;
        int chosen = -1;

        for (int j = 0; j < 24; ++j) {
          if (kPermutations[j][0] == 0 && kPermutations[j][1] == 1 &&
              kPermutations[j][2] == 2 && kPermutations[j][3] == 3) continue;

          const double trial = matrix[0][kPermutations[j][0]] + matrix[1][kPermutations[j][1]] +
                               matrix[2][kPermutations[j][2]] + matrix[3][kPermutations[j][3]];

          if (trial < best) { best = trial; chosen = j; }
        }

        if (chosen >= 0) {
          ++changed;
          for (int k = 0; k < 4; ++k) {
            match[H[k]] = R[kPermutations[chosen][k]];
            cost[H[k]] = matrix[k][kPermutations[chosen][k]];
          }
        }
      }

      ratio = (double)changed / (double)nObjects;
      if (trace != nullptr)
        trace->push_back(std::accumulate(costs.begin(), costs.end(), 0.) / (double)nObjects);
    }
  }

}


// ============================================================================


std::size_t otswap::internal::min_objects ()
{
  return (std::size_t)kNeighbours + 1;
}


// ============================================================================


otswap::Result otswap::internal::reconstruct (const std::vector<double>& tracers,
                                              const std::vector<double>& randoms,
                                              const std::vector<double>& mps,
                                              const double representativeMps,
                                              const Config& config,
                                              const unsigned seed,
                                              const BoxLimits* box,
                                              std::vector<std::vector<double>>* sweepCosts,
                                              const bool periodic,
                                              const bool means)
{
  const std::size_t nObjects = tracers.size() / 3;

  if (periodic && box == nullptr)
    throw Error("a periodic reconstruction needs the limits of its box");
  const unsigned nRealizations = config.nRealizations;
  const bool generate = randoms.empty();

  if (nObjects < min_objects())
    throw Error("the catalog holds " + std::to_string(nObjects) + " objects; at least " +
                std::to_string(min_objects()) + " are required, because each tracer needs " +
                std::to_string(kNeighbours) + " distinct neighbours");

  if (!(representativeMps > 0.) || !std::isfinite(representativeMps))
    throw Error("the representative mean particle separation is " +
                std::to_string(representativeMps) + "; it must be positive");

  const double cellsize = config.cellSize * representativeMps;

  std::vector<double> tracerX(nObjects), tracerY(nObjects), tracerZ(nObjects);
  for (std::size_t i = 0; i < nObjects; ++i) {
    tracerX[i] = tracers[3*i];
    tracerY[i] = tracers[3*i+1];
    tracerZ[i] = tracers[3*i+2];
  }

  double minX, maxX, minY, maxY, minZ, maxZ;
  if (box == nullptr) {
    const auto xr = std::minmax_element(tracerX.begin(), tracerX.end());
    const auto yr = std::minmax_element(tracerY.begin(), tracerY.end());
    const auto zr = std::minmax_element(tracerZ.begin(), tracerZ.end());
    minX = *xr.first; maxX = *xr.second;
    minY = *yr.first; maxY = *yr.second;
    minZ = *zr.first; maxZ = *zr.second;
  }
  else {
    minX = box->lower[0]; maxX = box->upper[0];
    minY = box->lower[1]; maxY = box->upper[1];
    minZ = box->lower[2]; maxZ = box->upper[2];
  }

  std::array<std::array<double,2>,3> extent {{{minX, maxX}, {minY, maxY}, {minZ, maxZ}}};
  if (!periodic)
    for (std::size_t i = 0; i < randoms.size(); ++i) {
      std::array<double,2>& axis = extent[i % 3];
      axis[0] = std::min(axis[0], randoms[i]);
      axis[1] = std::max(axis[1], randoms[i]);
    }

  Periods periods {};
  if (periodic)
    for (int k = 0; k < 3; ++k) {
      periods.period[k] = box->upper[k] - box->lower[k];
      periods.half[k] = 0.5 * periods.period[k];
    }

  const meshsearch::MeshGrid tracerGrid =
    mesh_grid(tracerX, tracerY, tracerZ, cellsize, extent, periodic);
  const std::array<std::array<double,2>,3> lims = tracerGrid.get_lims();

  std::vector<unsigned> neighbours((std::size_t)nObjects * kNeighbours);
  {
    std::exception_ptr failure;
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nObjects; ++i) {
      try {
        const std::vector<unsigned> found = tracerGrid.nearestObjects(kNeighbours, (unsigned)i);
        if (found.size() != kNeighbours)
          throw Error("a tracer has " + std::to_string(found.size()) + " neighbours; " +
                      std::to_string(kNeighbours) + " are required");
        std::copy(found.begin(), found.end(), neighbours.begin() + (std::ptrdiff_t)(i * kNeighbours));
      }
      catch (...) {
#pragma omp critical
        if (!failure) failure = std::current_exception();
      }
    }
    if (failure) rethrow(failure);
  }

  std::vector<unsigned> supply;
  if (!generate) {
    const std::size_t nRandoms = randoms.size() / 3;
    supply.resize(nRandoms);
    std::iota(supply.begin(), supply.end(), 0u);
    std::mt19937 engine = stream(seed, 0, Stream::Supply);
    shuffle(supply, engine);
  }

  Result result;
  result.nObjects = nObjects;
  result.nRealizations = nRealizations;
  result.displacement.assign(3 * (std::size_t)nRealizations * nObjects, 0.);

  if (sweepCosts != nullptr) sweepCosts->assign(nRealizations, {});

  std::exception_ptr failure;

#pragma omp parallel for schedule(static)
  for (unsigned rec = 0; rec < nRealizations; ++rec) {

    try {

      std::vector<double> randomX(nObjects), randomY(nObjects), randomZ(nObjects);

      if (generate) {
        std::mt19937 engine = stream(seed, rec, Stream::Randoms);
        auto draw = [&engine] (const double lo, const double hi) {
          if (lo < hi) return uniformReal(engine, lo, hi);
          engine.discard(2);
          return lo;
        };
        for (std::size_t i = 0; i < nObjects; ++i) {
          randomX[i] = draw(minX, maxX);
          randomY[i] = draw(minY, maxY);
          randomZ[i] = draw(minZ, maxZ);
        }
      }
      else {
        for (std::size_t i = 0; i < nObjects; ++i) {
          const unsigned at = supply[(std::size_t)rec * nObjects + i];
          randomX[i] = randoms[3*(std::size_t)at];
          randomY[i] = randoms[3*(std::size_t)at+1];
          randomZ[i] = randoms[3*(std::size_t)at+2];
        }
      }

      meshsearch::MeshGrid::Options randomOptions;
      randomOptions.limits = lims;
      randomOptions.periodic = {periodic, periodic, periodic};
      meshsearch::MeshGrid randomGrid(randomX, randomY, randomZ, cellsize, randomOptions);

      std::mt19937 seeding = stream(seed, rec, Stream::Seeding);

      std::vector<unsigned> partner(nObjects);
      std::iota(partner.begin(), partner.end(), 0u);
      shuffle(partner, seeding);

      meshsearch::MeshGrid tracerCopy = tracerGrid;

      std::vector<unsigned> starters(nObjects);
      std::iota(starters.begin(), starters.end(), 0u);
      shuffle(starters, seeding);

      std::vector<std::uint8_t> paired(nObjects, 0);
      std::vector<unsigned> close;

      for (const unsigned i : starters) {
        if (paired[i]) continue;

        tracerCopy.closeObjects(tracerX[i], tracerY[i], tracerZ[i],
                                kSeedingRadiusInMps * mps[i], close);
        std::sort(close.begin(), close.end());

        unsigned toRemove = std::min((unsigned)kMaxPairsPerStarter, (unsigned)close.size());
        if (toRemove == 0) continue;

        const unsigned wanted =
          std::min((unsigned)close.size(), randomGrid.get_nObjects());
        std::vector<unsigned> closeRandom =
          randomGrid.nearestObjects(wanted, tracerX[i], tracerY[i], tracerZ[i]);

        const auto self = std::find(close.begin(), close.end(), i);
        if (self != close.end() && !closeRandom.empty()) {
          const unsigned r = closeRandom.front();

          partner[i] = r;
          paired[i] = 1;

          tracerCopy.removeObject(i);
          randomGrid.removeObject(r);

          close.erase(self);
          closeRandom.erase(closeRandom.begin());

          --toRemove;
        }

        while (toRemove > 0 && !close.empty() && !closeRandom.empty()) {
          const unsigned it = uniformInt(seeding, (unsigned)close.size());
          const unsigned ir = uniformInt(seeding, (unsigned)closeRandom.size());
          const unsigned t = close[it];
          const unsigned r = closeRandom[ir];

          partner[t] = r;
          paired[t] = 1;

          tracerCopy.removeObject(t);
          randomGrid.removeObject(r);

          close.erase(close.begin() + it);
          closeRandom.erase(closeRandom.begin() + ir);

          --toRemove;
        }
      }

      tracerCopy = {};
      randomGrid = {};

      const std::mt19937 engine = stream(seed, rec, Stream::Swap);
      std::vector<double> trace;
      std::vector<double>* traced = sweepCosts != nullptr ? &trace : nullptr;
      if (periodic)
        swap_loop<true>(tracerX, tracerY, tracerZ, randomX, randomY, randomZ, neighbours, periods,
                        config.convergence, engine, partner, traced);
      else
        swap_loop<false>(tracerX, tracerY, tracerZ, randomX, randomY, randomZ, neighbours, periods,
                         config.convergence, engine, partner, traced);

      if (sweepCosts != nullptr) (*sweepCosts)[rec] = trace;

      const std::size_t base = 3 * (std::size_t)rec * nObjects;
      for (std::size_t i = 0; i < nObjects; ++i) {
        const unsigned p = partner[i];
        result.displacement[base + 3*i]   = randomX[p] - tracerX[i];
        result.displacement[base + 3*i+1] = randomY[p] - tracerY[i];
        result.displacement[base + 3*i+2] = randomZ[p] - tracerZ[i];
      }
      if (periodic)
        for (std::size_t i = 0; i < 3 * nObjects; ++i)
          result.displacement[base + i] = minimum_image(result.displacement[base + i],
                                                        periods.period[i % 3]);

    }
    catch (...) {
#pragma omp critical
      if (!failure) failure = std::current_exception();
    }
  }

  if (failure) rethrow(failure);

  result.valid.assign((std::size_t)nRealizations * nObjects, 1);
  result.outsideRedshiftCut.assign(nObjects, 0);
  result.outsideMask.assign(nObjects, 0);
  result.outsideBox.assign(nObjects, 0);
  result.config = config;
  result.config.seed = seed;
  result.tracers = tracers;
  if (periodic) {
    result.box = *box;
    result.periodic = true;
  }
  if (means) summarize(result);

  return result;
}


// ============================================================================


std::string otswap::internal::reconstruction_line (const std::string& name, const Result& result)
{
  const std::size_t without =
    (std::size_t)std::count(result.validRealizations.begin(), result.validRealizations.end(), 0u);

  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << "otswap: " << name << ": " << result.nRealizations
      << (result.nRealizations == 1 ? " realization of " : " realizations of ") << result.nObjects
      << " tracers in " << std::fixed << std::setprecision(2) << result.elapsedSeconds << " s; "
      << without << " without a valid displacement\n";
  return out.str();
}


// ============================================================================


void otswap::internal::write_lines (std::ostream& log, const std::string& text)
{
  for (std::size_t start = 0, end; (end = text.find('\n', start)) != std::string::npos; start = end + 1)
    log << text.substr(start, end + 1 - start) << std::flush;
}


// ============================================================================


void otswap::internal::summarize (Result& result)
{
  const std::size_t nObjects = result.nObjects;
  const unsigned nRealizations = result.nRealizations;

  result.meanDisplacement.resize(3 * nObjects);
  result.validRealizations.resize(nObjects);
  std::fill(result.meanDisplacement.begin(), result.meanDisplacement.end(), 0.);
  std::fill(result.validRealizations.begin(), result.validRealizations.end(), 0u);

  for (std::size_t i = 0; i < nObjects; ++i) {
    unsigned kept = 0;
    for (unsigned rec = 0; rec < nRealizations; ++rec) {
      if (!result.valid[(std::size_t)rec * nObjects + i]) continue;
      const std::size_t at = 3 * ((std::size_t)rec * nObjects + i);
      result.meanDisplacement[3*i]   += result.displacement[at];
      result.meanDisplacement[3*i+1] += result.displacement[at+1];
      result.meanDisplacement[3*i+2] += result.displacement[at+2];
      ++kept;
    }

    result.validRealizations[i] = kept;

    if (kept == 0) {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      result.meanDisplacement[3*i] = result.meanDisplacement[3*i+1] =
        result.meanDisplacement[3*i+2] = nan;
    }
    else {
      result.meanDisplacement[3*i]   /= (double)kept;
      result.meanDisplacement[3*i+1] /= (double)kept;
      result.meanDisplacement[3*i+2] /= (double)kept;
    }
  }

  if (result.tracers.size() == 3 * nObjects) {
    result.lagrangian.resize(3 * nObjects);
    for (std::size_t k = 0; k < 3 * nObjects; ++k)
      result.lagrangian[k] = result.tracers[k] + result.meanDisplacement[k];
    if (result.periodic) wrap_into(result.lagrangian, *result.box);
  }
  else result.lagrangian.clear();

  if (result.geometry == Geometry::Lightcone && result.distances && !result.lagrangian.empty())
    sky_of_positions(result.lagrangian, *result.distances, result.lagrangianSky);
  else result.lagrangianSky.clear();
}


// ============================================================================


void otswap::recomputeMeans (Result& result)
{
  const std::size_t nDisplacements = (std::size_t)result.nRealizations * result.nObjects;

  if (result.displacement.size() != 3 * nDisplacements || result.valid.size() != nDisplacements)
    throw Error("the result is malformed: displacement and valid hold " +
                std::to_string(result.displacement.size()) + " and " +
                std::to_string(result.valid.size()) + " entries for " +
                std::to_string(result.nRealizations) + " realizations of " +
                std::to_string(result.nObjects) + " objects");

  for (std::size_t k = 0; k < nDisplacements; ++k)
    if (result.valid[k] > 1)
      throw Error("the result is malformed: valid holds " + std::to_string(result.valid[k]) +
                  " at entry " + std::to_string(k) + " (from 0); every entry must be 0 or 1");

  internal::check_flags(result);
  internal::summarize(result);
}


// ============================================================================


namespace {

  /// Both box overloads; a null mps is computed from the box: the limits
  /// when given, the bounding box of the tracers and the given randoms
  /// otherwise. With limits, the tracers and randoms outside them are left
  /// out before anything else is computed; in a periodic box they are
  /// wrapped into the limits instead, and every one is kept.
  otswap::Result box (const std::vector<double>& tracers, const std::vector<double>& randoms,
                      const double* givenMps, const otswap::Config& config,
                      const std::optional<otswap::BoxLimits>& limits, const bool periodic)
  {
    using namespace otswap;

    const auto start = std::chrono::steady_clock::now();
    const unsigned seed = internal::check_config(config);

    const std::size_t nObjects = internal::check_coordinates(tracers, "the tracer array");
    const std::size_t nRandoms = internal::check_coordinates(randoms, "the random array", true);
    if (nObjects > (std::size_t)std::numeric_limits<unsigned>::max())
      throw Error("the catalog holds " + std::to_string(nObjects) + " objects; it must hold at most " +
                  std::to_string(std::numeric_limits<unsigned>::max()));

    const bool reporting = config.log != nullptr && config.verbosity != Verbosity::Silent;

    if (periodic && !limits)
      throw Error("a periodic box needs its limits; give them");

    SelectionCounts counts;
    internal::Selection tracerSelection, randomSelection;
    std::vector<double> wrappedTracers, wrappedRandoms;
    if (periodic) {
      internal::check_box(*limits);
      internal::check_periodic_box(*limits);
      wrappedTracers = tracers;
      wrappedRandoms = randoms;
      counts.boxLimits = *limits;
      counts.periodic = true;
      counts.tracers = nObjects;
      counts.tracersWrapped = internal::wrap_into(wrappedTracers, *limits);
      counts.randoms = nRandoms;
      counts.randomsWrapped = internal::wrap_into(wrappedRandoms, *limits);
    }
    else if (limits) {
      internal::check_box(*limits);
      tracerSelection = internal::select(tracers, *limits);
      randomSelection = internal::select(randoms, *limits);
      counts.boxLimits = *limits;
      counts.tracers = nObjects;
      counts.tracersOutsideBox = tracerSelection.nOutsideBox;
      counts.randoms = nRandoms;
      counts.randomsOutsideBox = randomSelection.nOutsideBox;
      if (randoms.empty()) internal::check_kept_tracers(counts, tracerSelection.keep.size());
      else internal::check_kept(counts, tracerSelection.keep.size(), randomSelection.keep.size(),
                                config.nRealizations);
    }
    if (reporting) internal::write_lines(*config.log, counts.message());

    const bool dropTracers = limits && !periodic && tracerSelection.keep.size() != nObjects;
    const bool dropRandoms = limits && !periodic && randomSelection.keep.size() != nRandoms;
    const std::vector<double> keptTracers = dropTracers ? internal::rows(tracers, tracerSelection.keep)
                                                        : std::vector<double>();
    const std::vector<double> keptRandoms = dropRandoms ? internal::rows(randoms, randomSelection.keep)
                                                        : std::vector<double>();
    const std::vector<double>& t = periodic ? wrappedTracers : dropTracers ? keptTracers : tracers;
    const std::vector<double>& r = periodic ? wrappedRandoms : dropRandoms ? keptRandoms : randoms;
    const std::size_t nKept = t.size() / 3;

    BoxLimits span = limits ? *limits : internal::bounding_box(tracers);
    if (!limits && !randoms.empty()) {
      const BoxLimits given = internal::bounding_box(randoms);
      for (std::size_t k = 0; k < 3; ++k) {
        span.lower[k] = std::min(span.lower[k], given.lower[k]);
        span.upper[k] = std::max(span.upper[k], given.upper[k]);
      }
    }

    double mps;
    if (givenMps != nullptr) mps = *givenMps;
    else {
      const double volume = (span.upper[0] - span.lower[0]) * (span.upper[1] - span.lower[1]) *
                            (span.upper[2] - span.lower[2]);
      mps = internal::det_pow(volume / (double)nKept, 1./3.);
    }
    if (!limits && givenMps == nullptr && !(mps > 0.))
      throw Error(std::string(randoms.empty() ? "the tracers span"
                                              : "the tracers and the given randoms span") +
                  " no volume (all of them in one plane, on one line or at one point), so their mean "
                  "particle separation cannot be computed; give it");
    if (!std::isfinite(mps) || mps <= 0.)
      throw Error("the mean particle separation is " + std::to_string(mps) +
                  "; it must be positive");

    if ((!limits || periodic) && !randoms.empty())
      internal::check_random_supply(nRandoms, nObjects, config.nRealizations);

    if (reporting && config.verbosity == Verbosity::Detailed) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      out << "otswap: mean particle separation " << mps << " in coordinate units ("
          << (givenMps != nullptr ? "given"
              : periodic ? "from the periodic box"
              : limits ? "from the box limits"
              : randoms.empty() ? "from the bounding box of the tracers"
              : "from the bounding box of the tracers and the given randoms") << ")\n";
      internal::write_lines(*config.log, out.str());
    }

    const std::vector<double> perObject(nKept, mps);

    Result result = internal::reconstruct(t, r, perObject, mps, config, seed,
                                          limits ? &*limits : nullptr, nullptr, periodic, !dropTracers);
    if (dropTracers) {
      result = internal::expand(result, tracerSelection, nObjects);
      result.tracers = tracers;
      internal::summarize(result);
    }
    result.mps = mps;
    result.box = span;
    result.selection = counts;
    result.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (reporting) internal::write_lines(*config.log, internal::reconstruction_line("reconstructBox", result));

    return result;
  }

}


// ============================================================================


otswap::Result otswap::reconstructBox (const std::vector<double>& tracers,
                                       const std::vector<double>& randoms,
                                       const double mps,
                                       const Config& config,
                                       const std::optional<BoxLimits>& limits,
                                       const bool periodic)
{
  return box(tracers, randoms, &mps, config, limits, periodic);
}


// ============================================================================


otswap::Result otswap::reconstructBox (const std::vector<double>& tracers,
                                       const std::vector<double>& randoms,
                                       const Config& config,
                                       const std::optional<BoxLimits>& limits,
                                       const bool periodic)
{
  return box(tracers, randoms, nullptr, config, limits, periodic);
}
