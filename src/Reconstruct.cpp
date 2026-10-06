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
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <meshsearch/MeshGrid.h>

#include "internal.h"

namespace {

  // Properties of the algorithm, not knobs: none is exposed in the public
  // header.

  // Radius, in units of the local mean particle separation, of the ball
  // around a starter from which the seeding pass draws tracers to pair.
  constexpr double kSeedingRadiusInMps = 8.;

  // Radius, in units of the local mean particle separation, of the
  // neighbourhood the three swap partners are drawn from.
  constexpr double kNeighbourhoodRadiusInMps = 3.;

  // Neighbours held per tracer: the count a sphere of
  // kNeighbourhoodRadiusInMps mean separations holds at mean density,
  // (4/3) * pi * 3^3 = 113.097..., truncated to 113. It is also the size
  // of the list the swap partners are drawn from, and one less than the
  // smallest catalog accepted.
  constexpr unsigned kNeighbours =
    (unsigned)((4. / 3.) * 3.14159265358979323846 *
               kNeighbourhoodRadiusInMps * kNeighbourhoodRadiusInMps *
               kNeighbourhoodRadiusInMps);

  // Pairs the seeding pass makes around one starter before moving on.
  constexpr unsigned kMaxPairsPerStarter = 31;

  // The 24 permutations of four elements. The order is load-bearing: the
  // loop keeps the first permutation that strictly improves on the best so
  // far, so a different order can pick a different permutation of equal
  // cost.
  constexpr int kPermutations[24][4] = {
    {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1},
    {1,0,2,3}, {1,0,3,2}, {1,2,0,3}, {1,2,3,0}, {1,3,0,2}, {1,3,2,0},
    {2,0,1,3}, {2,0,3,1}, {2,1,0,3}, {2,1,3,0}, {2,3,0,1}, {2,3,1,0},
    {3,0,1,2}, {3,0,2,1}, {3,1,0,2}, {3,1,2,0}, {3,2,0,1}, {3,2,1,0}
  };

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
                                              std::vector<std::vector<double>>* sweepCosts)
{
  const std::size_t nObjects = tracers.size() / 3;
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

  const auto xr = std::minmax_element(tracerX.begin(), tracerX.end());
  const auto yr = std::minmax_element(tracerY.begin(), tracerY.end());
  const auto zr = std::minmax_element(tracerZ.begin(), tracerZ.end());
  const double minX = *xr.first, maxX = *xr.second;
  const double minY = *yr.first, maxY = *yr.second;
  const double minZ = *zr.first, maxZ = *zr.second;

  std::vector<std::vector<double>> extent {{minX, maxX}, {minY, maxY}, {minZ, maxZ}};
  for (std::size_t i = 0; i < randoms.size(); ++i) {
    std::vector<double>& axis = extent[i % 3];
    axis[0] = std::min(axis[0], randoms[i]);
    axis[1] = std::max(axis[1], randoms[i]);
  }

  const meshsearch::MeshGrid tracerGrid(tracerX, tracerY, tracerZ, cellsize, extent);
  const std::vector<std::vector<double>> lims = tracerGrid.get_lims();

  // An exception cannot leave a parallel region: the first one caught is
  // kept, whatever its type, and rethrown unchanged after the region. The
  // same holds for the realizations below.
  std::vector<std::vector<unsigned>> neighbours(nObjects);
  {
    std::exception_ptr failure;
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nObjects; ++i) {
      try {
        neighbours[i] = tracerGrid.nearestObjects(kNeighbours, (unsigned)i);
      }
      catch (...) {
#pragma omp critical
        if (!failure) failure = std::current_exception();
      }
    }
    if (failure) std::rethrow_exception(failure);
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
  result.matchedRandom.assign(3 * (std::size_t)nRealizations * nObjects, 0.);

  if (sweepCosts != nullptr) sweepCosts->assign(nRealizations, {});

  std::exception_ptr failure;

#pragma omp parallel for schedule(static)
  for (unsigned rec = 0; rec < nRealizations; ++rec) {

    try {

      std::vector<double> randomX(nObjects), randomY(nObjects), randomZ(nObjects);

      if (generate) {
        std::mt19937 engine = stream(seed, rec, Stream::Randoms);
        for (std::size_t i = 0; i < nObjects; ++i) {
          randomX[i] = uniform_real(engine, minX, maxX);
          randomY[i] = uniform_real(engine, minY, maxY);
          randomZ[i] = uniform_real(engine, minZ, maxZ);
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

      const meshsearch::MeshGrid randomGrid(randomX, randomY, randomZ, cellsize, lims);

      std::mt19937 seeding = stream(seed, rec, Stream::Seeding);

      std::vector<unsigned> partner(nObjects);
      std::iota(partner.begin(), partner.end(), 0u);
      shuffle(partner, seeding);

      meshsearch::MeshGrid tracerCopy = tracerGrid;
      meshsearch::MeshGrid randomCopy = randomGrid;

      std::vector<unsigned> starters(nObjects);
      std::iota(starters.begin(), starters.end(), 0u);
      shuffle(starters, seeding);

      std::vector<std::uint8_t> paired(nObjects, 0);

      for (const unsigned i : starters) {
        if (paired[i]) continue;

        // The grid returns the ball in cell order, which depends on the cell
        // size; sorted by index, the pairs drawn from it below do not. The
        // nearest-random list needs no sorting: it is ordered by distance,
        // nearest first, so only randoms at exactly equal distances could
        // come in an order that depends on the cell size.
        std::vector<unsigned> close =
          tracerCopy.closeObjects(tracerX[i], tracerY[i], tracerZ[i],
                                  kSeedingRadiusInMps * mps[i]);
        std::sort(close.begin(), close.end());

        unsigned toRemove = std::min((unsigned)kMaxPairsPerStarter, (unsigned)close.size());
        if (toRemove == 0) continue;

        const unsigned wanted =
          std::min((unsigned)close.size(), randomCopy.get_nObjects());
        std::vector<unsigned> closeRandom =
          randomCopy.nearestObjects(wanted, tracerX[i], tracerY[i], tracerZ[i]);

        // The starter is paired first, with its nearest unused random, and
        // the remaining pairs are drawn at random from the ball. Every
        // starter reached here is therefore paired, so every tracer leaves
        // the pass with a partner of its own and the matching is a
        // permutation, whatever the cap and the size of the ball.
        const auto self = std::find(close.begin(), close.end(), i);
        if (self != close.end() && !closeRandom.empty()) {
          const unsigned r = closeRandom.front();

          partner[i] = r;
          paired[i] = 1;

          tracerCopy.removeObject(i);
          randomCopy.removeObject(r);

          close.erase(self);
          closeRandom.erase(closeRandom.begin());

          --toRemove;
        }

        while (toRemove > 0 && !close.empty() && !closeRandom.empty()) {
          const unsigned it = uniform_int(seeding, (unsigned)close.size());
          const unsigned ir = uniform_int(seeding, (unsigned)closeRandom.size());
          const unsigned t = close[it];
          const unsigned r = closeRandom[ir];

          partner[t] = r;
          paired[t] = 1;

          tracerCopy.removeObject(t);
          randomCopy.removeObject(r);

          close.erase(close.begin() + it);
          closeRandom.erase(closeRandom.begin() + ir);

          --toRemove;
        }
      }

      tracerCopy = {};
      randomCopy = {};

      std::vector<double> cost(nObjects);
      for (std::size_t i = 0; i < nObjects; ++i) {
        const double dx = tracerX[i] - randomX[partner[i]];
        const double dy = tracerY[i] - randomY[partner[i]];
        const double dz = tracerZ[i] - randomZ[partner[i]];
        cost[i] = dx*dx + dy*dy + dz*dz;
      }

      std::mt19937 engine = stream(seed, rec, Stream::Swap);

      std::vector<unsigned> order(nObjects);
      std::iota(order.begin(), order.end(), 0u);

      unsigned H[4], R[4];
      double ratio = 1.;
      std::vector<double> trace;

      while (ratio > config.convergence) {

        shuffle(order, engine);
        unsigned changed = 0;

        for (const unsigned i : order) {
          H[0] = i;
          unsigned r1 = uniform_int(engine, kNeighbours);
          unsigned r2 = uniform_int(engine, kNeighbours);
          unsigned r3 = uniform_int(engine, kNeighbours);
          while (r1 == r2) r2 = uniform_int(engine, kNeighbours);
          while (r3 == r1 || r3 == r2) r3 = uniform_int(engine, kNeighbours);

          H[1] = neighbours[i][r1];
          H[2] = neighbours[i][r2];
          H[3] = neighbours[i][r3];

          for (int k = 0; k < 4; ++k) R[k] = partner[H[k]];

          double matrix[4][4];
          for (int k = 0; k < 4; ++k) {
            const double tx = tracerX[H[k]], ty = tracerY[H[k]], tz = tracerZ[H[k]];
            for (int m = 0; m < 4; ++m) {
              const double dx = tx - randomX[R[m]];
              const double dy = ty - randomY[R[m]];
              const double dz = tz - randomZ[R[m]];
              matrix[k][m] = dx*dx + dy*dy + dz*dz;
            }
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
              partner[H[k]] = R[kPermutations[chosen][k]];
              cost[H[k]] = matrix[k][kPermutations[chosen][k]];
            }
          }
        }

        ratio = (double)changed / (double)nObjects;
        if (sweepCosts != nullptr)
          trace.push_back(std::accumulate(cost.begin(), cost.end(), 0.) / (double)nObjects);
      }

      if (sweepCosts != nullptr) (*sweepCosts)[rec] = trace;

      const std::size_t base = 3 * (std::size_t)rec * nObjects;
      for (std::size_t i = 0; i < nObjects; ++i) {
        const unsigned p = partner[i];
        result.matchedRandom[base + 3*i]   = randomX[p];
        result.matchedRandom[base + 3*i+1] = randomY[p];
        result.matchedRandom[base + 3*i+2] = randomZ[p];
        result.displacement[base + 3*i]   = randomX[p] - tracerX[i];
        result.displacement[base + 3*i+1] = randomY[p] - tracerY[i];
        result.displacement[base + 3*i+2] = randomZ[p] - tracerZ[i];
      }

    }
    catch (...) {
#pragma omp critical
      if (!failure) failure = std::current_exception();
    }
  }

  if (failure) std::rethrow_exception(failure);

  result.valid.assign((std::size_t)nRealizations * nObjects, 1);
  summarize(result);

  return result;
}


// ============================================================================


void otswap::internal::summarize (Result& result)
{
  const std::size_t nObjects = result.nObjects;
  const unsigned nRealizations = result.nRealizations;

  // Resized only when the sizes differ, then overwritten in place: a filter
  // never moves the memory these arrays occupy.
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
}


// ============================================================================


otswap::Result otswap::reconstructBox (const std::vector<double>& tracers,
                                       const std::vector<double>& randoms,
                                       const double mps,
                                       const Config& config)
{
  const unsigned seed = internal::check_config(config);

  const std::size_t nObjects = internal::check_coordinates(tracers, "the tracer array");
  const std::size_t nRandoms = internal::check_coordinates(randoms, "the random array", true);

  if (!std::isfinite(mps) || mps <= 0.)
    throw Error("the mean particle separation is " + std::to_string(mps) +
                "; it must be positive");

  if (!randoms.empty())
    internal::check_random_supply(nRandoms, nObjects, config.nRealizations);

  const std::vector<double> perObject(nObjects, mps);

  return internal::reconstruct(tracers, randoms, perObject, mps, config, seed);
}
