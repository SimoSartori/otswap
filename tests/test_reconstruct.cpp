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
 *  @file tests/test_reconstruct.cpp
 *
 *  @brief The guarantees the reconstruction states: a valid permutation,
 *  a cost that never rises, a known optimum found, and a result that a
 *  seed pins down whatever the thread count.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "internal.h"
#include "check.h"

using namespace otswap;

namespace {

  // A cube of nSide^3 lattice sites, jittered so that no two tracers
  // coincide and the neighbour lookup has a well-defined answer.
  std::vector<double> lattice (const unsigned nSide, const double spacing,
                               const double jitter, const unsigned seed)
  {
    std::mt19937 rng(seed);
    std::vector<double> out;
    out.reserve(3ull * nSide * nSide * nSide);

    for (unsigned i = 0; i < nSide; ++i)
      for (unsigned j = 0; j < nSide; ++j)
        for (unsigned k = 0; k < nSide; ++k) {
          out.push_back(i * spacing + internal::uniform_real(rng, -jitter, jitter));
          out.push_back(j * spacing + internal::uniform_real(rng, -jitter, jitter));
          out.push_back(k * spacing + internal::uniform_real(rng, -jitter, jitter));
        }

    return out;
  }

  double total_cost (const Result& r, const std::vector<double>& tracers,
                     const unsigned rec)
  {
    const std::size_t base = 3 * (std::size_t)rec * r.nObjects;
    double sum = 0.;
    for (std::size_t i = 0; i < r.nObjects; ++i)
      for (int c = 0; c < 3; ++c) {
        const double d = r.matchedRandom[base + 3*i + c] - tracers[3*i + c];
        sum += d*d;
      }
    return sum;
  }

}

int main ()
{
  const unsigned nSide = 8;                       // 512 objects, above the 114 floor
  const double spacing = 10.;
  const std::vector<double> tracers = lattice(nSide, spacing, 1.0, 11);
  const std::size_t nObjects = tracers.size() / 3;

  group("a catalog below the neighbour floor is refused, loudly");
  {
    const std::vector<double> tiny(3 * 50, 0.);
    Config c;
    check_throws([&] { reconstructBox(tiny, {}, spacing, c); },
                 "fewer objects than the neighbourhood needs raises");
  }

  group("malformed input is refused before any work is done");
  {
    Config c;
    check_throws([&] { reconstructBox({}, {}, spacing, c); }, "an empty tracer array");
    check_throws([&] { reconstructBox({1., 2.}, {}, spacing, c); },
                 "a tracer array whose size is not a multiple of three");
    check_throws([&] { reconstructBox(tracers, {}, 0., c); },
                 "a non-positive mean particle separation");
    check_throws([&] { reconstructBox(tracers, {}, -1., c); }, "a negative separation");

    std::vector<double> nonFinite = tracers;
    nonFinite[7] = std::nan("");
    check_throws([&] { reconstructBox(nonFinite, {}, spacing, c); },
                 "a non-finite coordinate");

    Config two;
    two.nRealizations = 2;
    check_throws([&] { reconstructBox(tracers, tracers, spacing, two); },
                 "randoms too few for the realizations requested");
  }

  group("the matching is a permutation: every random is used at most once");
  {
    Config c;
    c.nRealizations = 3;
    c.seed = 2024;
    c.convergence = 1.e-2;

    // One lattice per realization, each with its own jitter, so that no
    // two randoms in the supply share a position and a repeated position
    // in the output can only mean a repeated match.
    std::vector<double> supply;
    for (unsigned r = 0; r < c.nRealizations; ++r) {
      const std::vector<double> block = lattice(nSide, spacing, 4.0, 99 + r);
      supply.insert(supply.end(), block.begin(), block.end());
    }

    const Result result = reconstructBox(tracers, supply, spacing, c);

    check(result.nObjects == nObjects, "the result reports the object count");
    check(result.nRealizations == 3, "the result reports the realization count");
    check(result.displacement.size() == 3 * 3 * nObjects, "displacement has the stated size");
    check(result.matchedRandom.size() == 3 * 3 * nObjects, "matchedRandom has the stated size");
    check(result.meanDisplacement.size() == 3 * nObjects,
          "meanDisplacement has the stated size");
    check(result.validRealizations.size() == nObjects,
          "validRealizations has the stated size");
    check(result.valid.size() == 3 * nObjects, "valid has the stated size");

    for (const std::uint8_t v : result.valid)
      check(v == 1, "every displacement is valid before a filter is applied");
    for (const unsigned v : result.validRealizations)
      check(v == 3, "every object counts every realization before a filter is applied");

    // Without a filter the mean is over every realization, summed in
    // realization order.
    bool plainMean = true;
    for (std::size_t i = 0; i < nObjects; ++i)
      for (int k = 0; k < 3; ++k) {
        double sum = 0.;
        for (unsigned rec = 0; rec < 3; ++rec)
          sum += result.displacement[3 * ((std::size_t)rec * nObjects + i) + k];
        if (result.meanDisplacement[3*i+k] != sum / 3.) plainMean = false;
      }
    check(plainMean, "the mean displacement is the plain mean over the realizations");

    for (unsigned rec = 0; rec < 3; ++rec) {
      const std::size_t base = 3 * (std::size_t)rec * nObjects;
      std::set<std::array<double, 3>> used;
      for (std::size_t i = 0; i < nObjects; ++i)
        used.insert({result.matchedRandom[base+3*i],
                     result.matchedRandom[base+3*i+1],
                     result.matchedRandom[base+3*i+2]});
      check(used.size() == nObjects,
            "no random is matched to two tracers within a realization");
    }
  }

  group("the Lagrangian position is the coordinate plus the displacement");
  {
    Config c;
    c.seed = 7;
    const Result result = reconstructBox(tracers, {}, spacing, c);

    for (std::size_t i = 0; i < nObjects; ++i)
      for (int k = 0; k < 3; ++k)
        check_close(result.matchedRandom[3*i+k],
                    tracers[3*i+k] + result.displacement[3*i+k], 1.e-12,
                    "matchedRandom equals the tracer plus its displacement");

    for (std::size_t i = 0; i < nObjects; ++i)
      for (int k = 0; k < 3; ++k)
        check_close(result.meanDisplacement[3*i+k], result.displacement[3*i+k], 1.e-12,
                    "with one realization the mean is that realization");
  }

  group("the total cost never rises across a sweep, and the loop terminates");
  {
    Config c;
    c.nRealizations = 2;
    c.seed = 31337;
    c.convergence = 1.e-3;

    const std::vector<double> mps(nObjects, spacing);
    std::vector<std::vector<double>> trace;

    internal::reconstruct(tracers, {}, mps, spacing, c,
                          internal::check_config(c), &trace);

    check(trace.size() == 2, "each realization reports its sweeps");

    for (const auto& perRealization : trace) {
      check(!perRealization.empty(), "at least one sweep ran");
      for (std::size_t s = 1; s < perRealization.size(); ++s)
        check(perRealization[s] <= perRealization[s-1] + 1.e-9 * perRealization[0],
              "the mean cost does not rise from one sweep to the next");
    }
  }

  group("a tighter convergence never leaves a worse matching");
  {
    Config loose, tight;
    loose.seed = tight.seed = 555;
    loose.convergence = 5.e-2;
    tight.convergence = 1.e-4;

    const Result a = reconstructBox(tracers, {}, spacing, loose);
    const Result b = reconstructBox(tracers, {}, spacing, tight);

    check(total_cost(b, tracers, 0) <= total_cost(a, tracers, 0),
          "the tighter run ends at or below the looser one");
  }

  group("a configuration with a known optimum is approached, to a local optimum");
  {
    // The randoms are the tracers themselves, so the optimal assignment
    // is the identity and its cost is exactly zero: the answer is known
    // in advance and the search can be measured against it.
    //
    // 4-opt local search need not reach it: the loop ends when a sweep
    // changes few enough pairs, which a matching that is not optimal can
    // satisfy. What is checked here is what the algorithm guarantees: a
    // large, monotone reduction towards the known answer, and a majority
    // of pairs recovered exactly.
    Config c;
    c.seed = 4;
    c.convergence = 1.e-4;

    const std::vector<double> mps(nObjects, spacing);
    std::vector<std::vector<double>> trace;
    const Result result =
      internal::reconstruct(tracers, tracers, mps, spacing, c, 4, &trace);

    const double seeded = trace[0].front();
    const double reached = total_cost(result, tracers, 0) / (double)nObjects;

    check(reached < 0.1 * seeded,
          "the search removes over nine tenths of the cost it started with");
    check(reached < 0.5 * spacing * spacing,
          "it ends well below the squared mean separation");
    check(trace[0].back() <= trace[0].front(),
          "and it never ends above where it started");

    std::size_t exact = 0;
    for (std::size_t i = 0; i < nObjects; ++i) {
      double d = 0.;
      for (int k = 0; k < 3; ++k) {
        const double e = result.matchedRandom[3*i+k] - tracers[3*i+k];
        d += e*e;
      }
      if (d == 0.) ++exact;
    }
    check(exact > nObjects / 2, "more than half the tracers find their own copy");
  }

  group("a seed pins the result down, and the thread count does not move it");
  {
    Config c;
    c.nRealizations = 4;
    c.seed = 98765;
    c.convergence = 1.e-3;

#ifdef _OPENMP
    const int restore = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    const Result one = reconstructBox(tracers, {}, spacing, c);

#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const Result many = reconstructBox(tracers, {}, spacing, c);

#ifdef _OPENMP
    omp_set_num_threads(restore);
#endif

    check(one.displacement == many.displacement,
          "the displacements are bit-identical at one and four threads");
    check(one.matchedRandom == many.matchedRandom,
          "the matched randoms are bit-identical at one and four threads");
    check(one.meanDisplacement == many.meanDisplacement,
          "the mean displacement is bit-identical at one and four threads");

    const Result again = reconstructBox(tracers, {}, spacing, c);
    check(one.displacement == again.displacement, "the same seed repeats exactly");

    // The same with randoms given, so that their assignment to
    // realizations is drawn too.
    std::vector<double> supply;
    for (unsigned r = 0; r < c.nRealizations + 1; ++r) {
      const std::vector<double> block = lattice(nSide, spacing, 4.0, 300 + r);
      supply.insert(supply.end(), block.begin(), block.end());
    }
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const Result givenOne = reconstructBox(tracers, supply, spacing, c);
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const Result givenMany = reconstructBox(tracers, supply, spacing, c);
#ifdef _OPENMP
    omp_set_num_threads(restore);
#endif
    check(givenOne.matchedRandom == givenMany.matchedRandom,
          "with randoms given, the result is bit-identical at one and four threads");

    Config other = c;
    other.seed = 98766;
    const Result different = reconstructBox(tracers, {}, spacing, other);
    check(one.displacement != different.displacement, "a different seed gives a different run");
  }

  group("the cell size affects speed only, never the result");
  {
    Config c;
    c.nRealizations = 2;
    c.seed = 4321;

    std::vector<double> supply;
    for (unsigned r = 0; r < c.nRealizations; ++r) {
      const std::vector<double> block = lattice(nSide, spacing, 4.0, 500 + r);
      supply.insert(supply.end(), block.begin(), block.end());
    }

    const DistanceTable distances(0.3, 0.7, -1., 0., 0., 1.5, 4000);
    std::mt19937 rng(77);
    std::vector<double> sky(3 * 1200), randomSky(3 * 2400);
    for (std::vector<double>* s : {&sky, &randomSky})
      for (std::size_t i = 0; i < s->size(); i += 3) {
        (*s)[i]   = internal::uniform_real(rng, 0.2, 0.6);
        (*s)[i+1] = internal::uniform_real(rng, -0.2, 0.2);
        (*s)[i+2] = internal::uniform_real(rng, 0.5, 0.9);
      }

    c.cellSize = 4.;
    const Result boxRef = reconstructBox(tracers, {}, spacing, c);
    const Result givenRef = reconstructBox(tracers, supply, spacing, c);
    const Result coneRef = reconstructLightcone(sky, randomSky, 1500., 3, distances, c);

    for (const double cellSize : {1., 2., 8.}) {
      c.cellSize = cellSize;
      const std::string at = " at cell size " + std::to_string((int)cellSize) + " as at 4";
      const Result box = reconstructBox(tracers, {}, spacing, c);
      const Result given = reconstructBox(tracers, supply, spacing, c);
      const Result cone = reconstructLightcone(sky, randomSky, 1500., 3, distances, c);
      check(box.matchedRandom == boxRef.matchedRandom &&
            box.displacement == boxRef.displacement,
            "box, generated randoms: the same result" + at);
      check(given.matchedRandom == givenRef.matchedRandom &&
            given.displacement == givenRef.displacement,
            "box, randoms given: the same result" + at);
      check(cone.matchedRandom == coneRef.matchedRandom &&
            cone.displacement == coneRef.displacement,
            "lightcone: the same result" + at);
    }
  }

  group("randoms outside the tracers' bounding box are accepted");
  {
    // The randoms reach a full spacing beyond the tracers on every side.
    Config c;
    c.seed = 8080;
    c.nRealizations = 2;
    std::vector<double> wide;
    std::mt19937 rng(17);
    const double lo = -1.5 * spacing, hi = (nSide + 0.5) * spacing;
    for (std::size_t i = 0; i < 2 * nObjects; ++i)
      for (int k = 0; k < 3; ++k) wide.push_back(internal::uniform_real(rng, lo, hi));

    Result r;
    bool ran = true;
    try {
      r = reconstructBox(tracers, wide, spacing, c);
    }
    catch (const Error&) {
      ran = false;
    }
    check(ran, "a box reconstruction with such randoms runs");

    if (ran)
      for (unsigned rec = 0; rec < 2; ++rec) {
        std::set<std::array<double, 3>> used;
        const std::size_t base = 3 * (std::size_t)rec * nObjects;
        for (std::size_t i = 0; i < nObjects; ++i)
          used.insert({r.matchedRandom[base+3*i], r.matchedRandom[base+3*i+1],
                       r.matchedRandom[base+3*i+2]});
        check(used.size() == nObjects, "and its matching is a permutation");
      }
  }

  group("mps(z) on a synthetic catalog recovers the density it was built from");
  {
    // A shell of constant comoving number density, so that every bin must
    // return the same mean particle separation, n^(-1/3).
    const DistanceTable distances(0.3, 0.7, -1., 0., 0., 1.5, 4000);

    const double skyAreaDeg2 = 2000.;
    const double fullSky = 4. * 3.14159265358979323846 * (180./3.14159265358979323846)
                              * (180./3.14159265358979323846);
    const double skyFraction = skyAreaDeg2 / fullSky;

    const unsigned nBins = 6;
    const double zLo = 0.4, zHi = 1.0;
    const double step = (zHi - zLo) / nBins;
    const double density = 2.e-3;                 // objects per (Mpc/h)^3

    std::vector<double> sky;
    std::mt19937 rng(5);

    for (unsigned b = 0; b < nBins; ++b) {
      const double lo = zLo + b*step, hi = zLo + (b+1)*step;
      const double dLo = distances.distanceAt(lo), dHi = distances.distanceAt(hi);
      const double volume = skyFraction * (4.*3.14159265358979323846/3.)
                            * (dHi*dHi*dHi - dLo*dLo*dLo);
      const unsigned count = (unsigned)(density * volume);

      for (unsigned i = 0; i < count; ++i) {
        const double t = (i + 0.5) / (double)count;
        sky.push_back(internal::uniform_real(rng, 0., 1.));
        sky.push_back(internal::uniform_real(rng, -0.1, 0.1));
        sky.push_back(lo + t * (hi - lo));
      }
    }

    // Pin the observed range so the measured bin edges are the intended
    // ones; without this the first and last bin would be shifted.
    sky[2] = zLo;
    sky[sky.size()-1] = zHi;

    const internal::MpsProfile profile =
      internal::mps_profile(sky, skyAreaDeg2, nBins, distances);

    const double expected = std::pow(density, -1./3.);
    check(profile.mps.size() == nBins, "one node per bin");
    for (unsigned b = 0; b < nBins; ++b)
      check_close(profile.mps[b], expected, 0.02 * expected,
                  "the recovered separation matches the density it was built from");

    const double rep = internal::representative(profile);
    check_close(rep, expected, 0.02 * expected,
                "the representative separation is that same value");
  }

  group("an under-populated redshift bin is refused, naming the bin and the count");
  {
    const DistanceTable distances(0.3, 0.7, -1., 0., 0., 2., 2000);

    std::vector<double> sky;
    for (unsigned i = 0; i < 400; ++i) {
      sky.push_back(0.1);
      sky.push_back(0.);
      sky.push_back(0.5 + 0.5 * (double)i / 400.);
    }

    // 400 tracers over 40 bins leaves 10 per bin, far below the floor.
    check_throws([&] { internal::mps_profile(sky, 1000., 40, distances); },
                 "a thinly populated bin raises");

    try {
      internal::mps_profile(sky, 1000., 40, distances);
    }
    catch (const Error& e) {
      const std::string message = e.what();
      check(message.find("bin 0") != std::string::npos, "the message names the bin");
      check(message.find("supports at most") != std::string::npos,
            "the message says how many bins the catalog supports");
    }

    check_throws([&] { internal::mps_profile(sky, 1000., 0, distances); },
                 "zero bins raise");
    check_throws([&] { internal::mps_profile(sky, -1., 4, distances); },
                 "a non-positive sky area raises");
    check_throws([&] { internal::mps_profile(sky, 1.e9, 4, distances); },
                 "a sky area beyond the whole sky raises");
  }

  group("the lightcone entry points run, and agree with each other");
  {
    const DistanceTable distances(0.3, 0.7, -1., 0., 0., 1.5, 4000);

    const double skyAreaDeg2 = 1500.;
    const unsigned nBins = 4;
    const double zLo = 0.5, zHi = 0.9;

    // 1600 tracers over four bins leaves 400 in each, above the floor.
    const std::size_t nLc = 1600;
    std::mt19937 rng(21);

    std::vector<double> sky(3 * nLc);
    for (std::size_t i = 0; i < nLc; ++i) {
      sky[3*i]   = internal::uniform_real(rng, 0.2, 0.6);
      sky[3*i+1] = internal::uniform_real(rng, -0.2, 0.2);
      sky[3*i+2] = zLo + (zHi - zLo) * (double)i / (double)(nLc - 1);
    }

    std::vector<double> randomSky(3 * nLc);
    for (std::size_t i = 0; i < nLc; ++i) {
      randomSky[3*i]   = internal::uniform_real(rng, 0.2, 0.6);
      randomSky[3*i+1] = internal::uniform_real(rng, -0.2, 0.2);
      randomSky[3*i+2] = internal::uniform_real(rng, zLo, zHi);
    }

    Config c;
    c.seed = 606;
    c.convergence = 1.e-2;

    const Result fromSky =
      reconstructLightcone(sky, randomSky, skyAreaDeg2, nBins, distances, c);

    check(fromSky.nObjects == nLc, "the sky overload reports the object count");
    check(fromSky.displacement.size() == 3 * nLc, "and returns a field of the right size");

    // The same run with the Cartesian coordinates supplied must give the
    // same answer: the second overload only skips the conversion, and
    // toCartesian is the conversion the first one performs.
    const std::vector<double> cart = toCartesian(sky, distances);
    const std::vector<double> randomCart = toCartesian(randomSky, distances);

    bool sameConversion = true;
    for (std::size_t i = 0; i < nLc; ++i) {
      double x, y, z;
      internal::to_cartesian(sky[3*i], sky[3*i+1], distances.distanceAt(sky[3*i+2]), x, y, z);
      if (x != cart[3*i] || y != cart[3*i+1] || z != cart[3*i+2]) sameConversion = false;
    }
    check(sameConversion, "toCartesian is the per-object conversion, bit for bit");

    const Result fromCartesian =
      reconstructLightcone(cart, randomCart, sky, randomSky, skyAreaDeg2, nBins,
                           distances, c);

    check(fromSky.displacement == fromCartesian.displacement,
          "both overloads produce the same field, bit for bit");

    // Exactly a permutation: the seeding pass pairs every tracer.
    std::set<std::array<double, 3>> used;
    for (std::size_t i = 0; i < nLc; ++i)
      used.insert({fromSky.matchedRandom[3*i], fromSky.matchedRandom[3*i+1],
                   fromSky.matchedRandom[3*i+2]});
    check(used.size() == nLc, "the lightcone matching is a permutation too");

    check_throws([&] {
      reconstructLightcone(cart, randomCart, sky, randomSky, skyAreaDeg2, 20,
                           distances, c);
    }, "too many bins for the catalog raises");

    // Randoms spread wider on the sky and deeper in redshift than the
    // tracers, so that many lie outside the tracers' bounding box.
    std::vector<double> wideSky(3 * nLc);
    for (std::size_t i = 0; i < nLc; ++i) {
      wideSky[3*i]   = internal::uniform_real(rng, 0.1, 0.7);
      wideSky[3*i+1] = internal::uniform_real(rng, -0.3, 0.3);
      wideSky[3*i+2] = internal::uniform_real(rng, zLo - 0.05, zHi + 0.05);
    }
    bool wideRan = true;
    try {
      const Result wide = reconstructLightcone(sky, wideSky, skyAreaDeg2, nBins, distances, c);
      std::set<std::array<double, 3>> wideUsed;
      for (std::size_t i = 0; i < nLc; ++i)
        wideUsed.insert({wide.matchedRandom[3*i], wide.matchedRandom[3*i+1],
                         wide.matchedRandom[3*i+2]});
      check(wideUsed.size() == nLc, "a lightcone with randoms outside the tracers' box matches a permutation");
    }
    catch (const Error&) {
      wideRan = false;
    }
    check(wideRan, "a lightcone with randoms outside the tracers' bounding box runs");

    std::vector<double> shortSky(sky.begin(), sky.end() - 3);
    check_throws([&] {
      reconstructLightcone(cart, randomCart, shortSky, randomSky, skyAreaDeg2, nBins,
                           distances, c);
    }, "Cartesian and sky arrays of different length raise");

    const DistanceTable narrow(0.3, 0.7, -1., 0., 0.6, 0.8, 500);
    check_throws([&] {
      reconstructLightcone(sky, randomSky, skyAreaDeg2, nBins, narrow, c);
    }, "redshifts outside the distance table raise");
  }

  return report("test_reconstruct");
}
