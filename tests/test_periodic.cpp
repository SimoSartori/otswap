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
 *  @file tests/test_periodic.cpp
 *
 *  @brief The guarantees of a periodic box: the checks of its limits, the
 *  wrap of tracers and randoms and its counts, a valid permutation with
 *  minimum-image displacements, a cost that never rises, a known optimum
 *  across the faces approached, the wrapped Lagrangian positions, and
 *  results that a seed pins down, whatever the thread count, to recorded
 *  hashes.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"
#include "otswap/Random.h"

#include "internal.h"
#include "check.h"

using namespace otswap;

namespace {

  // A cube of nSide^3 lattice sites at offset + i * spacing, jittered.
  std::vector<double> lattice (const unsigned nSide, const double spacing, const double offset,
                               const double jitter, const unsigned seed)
  {
    std::mt19937 rng(seed);
    std::vector<double> out;
    for (unsigned i = 0; i < nSide; ++i)
      for (unsigned j = 0; j < nSide; ++j)
        for (unsigned k = 0; k < nSide; ++k)
          for (const unsigned n : {i, j, k})
            out.push_back(offset + n * spacing + uniformReal(rng, -jitter, jitter));
    return out;
  }

  // FNV-1a, 64 bits, over the bytes of the doubles: a fingerprint of the
  // exact bits, the same on every little-endian platform.
  std::uint64_t fingerprint (const std::vector<double>& a)
  {
    std::uint64_t h = 14695981039346656037ull;
    for (const double x : a) {
      unsigned char bytes[sizeof(double)];
      std::memcpy(bytes, &x, sizeof(double));
      for (const unsigned char b : bytes) { h ^= b; h *= 1099511628211ull; }
    }
    return h;
  }

  void check_fingerprint (const std::vector<double>& a, const std::uint64_t expected,
                          const std::string& what)
  {
    const std::uint64_t got = fingerprint(a);
    char text[40];
    std::snprintf(text, sizeof text, "0x%016llx", (unsigned long long)got);
    check(got == expected, what + " (got " + text + ")");
  }

  bool inside (const std::vector<double>& a, const BoxLimits& b)
  {
    for (std::size_t i = 0; i < a.size(); ++i)
      if (!(a[i] >= b.lower[i % 3] && a[i] <= b.upper[i % 3])) return false;
    return true;
  }

  template <typename F>
  Result at_threads (const int threads, F&& f)
  {
#ifdef _OPENMP
    const int restore = omp_get_max_threads();
    omp_set_num_threads(threads);
    Result r = f();
    omp_set_num_threads(restore);
    return r;
#else
    (void)threads;
    return f();
#endif
  }

}

int main ()
{
  const unsigned nSide = 8;                       // 512 tracers
  const double spacing = 10.;
  BoxLimits cube;
  cube.lower = {0., 0., 0.};
  cube.upper = {80., 80., 80.};
  const std::vector<double> tracers = lattice(nSide, spacing, 5., 1., 11);
  const std::size_t nObjects = tracers.size() / 3;

  group("a periodic box needs limits whose sides are integer multiples of the shortest");
  {
    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 1;
    c.seed = 1;
    check(throws_naming([&] { reconstructBox(tracers, {}, spacing, c, std::nullopt, true); },
                        "a periodic box needs its limits"),
          "no limits");
    BoxLimits odd = cube;
    odd.upper[2] = 120.;
    check(throws_naming([&] { reconstructBox(tracers, {}, spacing, c, odd, true); },
                        "the side on axis z is not an integer multiple of the shortest one"),
          "80 x 80 x 120 is refused");
    check_throws([&] { reconstructBox(tracers, {}, c, odd, true); },
                 "also with the mps computed");
    BoxLimits flat = cube;
    flat.upper[1] = 0.;
    check(throws_naming([&] { reconstructBox(tracers, {}, spacing, c, flat, true); },
                        "the lower bound must be below the upper one"),
          "the limits are checked as without periodicity");
    BoxLimits tall = cube;
    tall.upper[2] = 160.;
    const Result r = reconstructBox(tracers, {}, spacing, c, tall, true);
    check(r.periodic && inside(r.tracers, tall) && inside(r.lagrangian, tall),
          "80 x 80 x 160 is accepted");
  }

  group("tracers and given randoms outside the limits are wrapped, counted and all kept");
  {
    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 2;
    c.seed = 7;
    std::vector<double> shifted = tracers;
    shifted[0] += 80.;                     // one period up on x
    shifted[4] -= 160.;                    // two periods down on y
    shifted[3*5+2] = 80.;                  // on the upper face: inside, kept as given
    std::vector<double> randoms = lattice(nSide, spacing, 2., 2., 21);
    const std::vector<double> more = lattice(nSide, spacing, 8., 2., 22);
    randoms.insert(randoms.end(), more.begin(), more.end());
    randoms[3] -= 80.;
    randoms[3*700+1] += 240.;

    std::ostringstream log;
    c.log = &log;
    c.verbosity = Verbosity::Normal;
    const Result r = reconstructBox(shifted, randoms, spacing, c, cube, true);
    check(r.selection.periodic && r.selection.boxLimits && r.selection.tracers == nObjects &&
          r.selection.tracersWrapped == 2 && r.selection.randoms == 2 * nObjects &&
          r.selection.randomsWrapped == 2 && r.selection.tracersOutsideBox == 0 &&
          r.selection.randomsOutsideBox == 0,
          "two tracers and two randoms wrapped, none outside");
    check(std::fabs(r.tracers[0] - tracers[0]) < 1.e-12 &&
          std::fabs(r.tracers[4] - tracers[4]) < 1.e-12 && r.tracers[3*5+2] == 80. &&
          r.tracers[1] == tracers[1],
          "Result::tracers holds the wrapped tracers; a point on a face is kept as given");
    check(r.periodic && r.box && r.box->upper == cube.upper && r.box->lower == cube.lower,
          "Result::periodic and Result::box");
    check(std::count(r.outsideBox.begin(), r.outsideBox.end(), 1) == 0 &&
          std::count(r.validRealizations.begin(), r.validRealizations.end(), 2u) == (long)nObjects,
          "no tracer left out: outsideBox all 0, every realization valid");
    check(log.str().find("otswap: wrapped 2 of 512 tracers into the periodic box [0, 80] x [0, 80] "
                         "x [0, 80]\notswap: wrapped 2 of 1024 randoms into the periodic box") == 0,
          "the report names what was wrapped, at Normal");
    check(log.str().find("0 without a valid displacement") != std::string::npos,
          "and the call line counts no loss");

    std::ostringstream none;
    c.log = &none;
    reconstructBox(tracers, {}, spacing, c, cube, true);
    check(none.str().find("otswap: wrapped 0 of 512 tracers into the periodic box") == 0 &&
          none.str().find("randoms into") == std::string::npos,
          "the tracer line is written with none wrapped; no random line with drawn randoms");
  }

  group("the matching is a permutation of the given randoms, by minimum-image displacements");
  {
    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 2;
    c.seed = 99;
    std::vector<double> randoms = lattice(nSide, spacing, 1., 1., 31);
    const std::vector<double> more = lattice(nSide, spacing, 9., 1., 32);
    randoms.insert(randoms.end(), more.begin(), more.end());
    const Result r = reconstructBox(tracers, randoms, spacing, c, cube, true);

    bool bounded = true;
    for (std::size_t k = 0; k < r.displacement.size(); ++k)
      bounded = bounded && r.displacement[k] >= -40. && r.displacement[k] < 40.;
    check(bounded, "every displacement component lies in [-L/2, L/2)");

    std::vector<int> used(2 * nObjects, 0);
    bool matched = true;
    for (std::size_t d = 0; d < 2 * nObjects; ++d) {
      const std::size_t i = d % nObjects;
      std::array<double, 3> p;
      for (int k = 0; k < 3; ++k)
        p[k] = internal::wrap_coordinate(r.tracers[3*i+k] + r.displacement[3*d+k], 0., 80.);
      std::size_t found = used.size();
      for (std::size_t j = 0; j < used.size() && found == used.size(); ++j) {
        double e = 0.;
        for (int k = 0; k < 3; ++k) {
          const double q = internal::minimum_image(p[k] - randoms[3*j+k], 80.);
          e += q*q;
        }
        if (e < 1.e-18) found = j;
      }
      matched = matched && found < used.size();
      if (found < used.size()) ++used[found];
    }
    check(matched, "each tracer plus its displacement, wrapped, is a given random");
    check(std::count(used.begin(), used.end(), 1) == (long)(2 * nObjects),
          "and each random is used once across the two realizations");
  }

  group("the total cost never rises across a sweep");
  {
    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 2;
    c.seed = 31337;
    const std::vector<double> mps(nObjects, spacing);
    std::vector<std::vector<double>> trace;
    internal::reconstruct(tracers, {}, mps, spacing, c, internal::check_config(c), &cube, &trace,
                          true);
    for (const auto& perRealization : trace) {
      check(!perRealization.empty(), "at least one sweep ran");
      for (std::size_t s = 1; s < perRealization.size(); ++s)
        check(perRealization[s] <= perRealization[s-1] + 1.e-9 * perRealization[0],
              "the mean cost does not rise from one sweep to the next");
    }
  }

  group("a known optimum across the faces is approached");
  {
    // The randoms are the tracers moved by +2 on x and wrapped, so the
    // tracers of the last layer find their copy across the face x = 80:
    // the optimum moves every tracer by (2, 0, 0), a cost reached only
    // through the minimum image.
    const std::vector<double> near = lattice(nSide, spacing, 9., 0.5, 41);
    std::vector<double> copies = near;
    for (std::size_t i = 0; i < nObjects; ++i)
      copies[3*i] = internal::wrap_coordinate(copies[3*i] + 2., 0., 80.);
    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 1;
    c.seed = 4;
    c.convergence = 1.e-4;
    const Result r = reconstructBox(near, copies, spacing, c, cube, true);
    std::size_t exact = 0, across = 0, acrossExact = 0;
    for (std::size_t i = 0; i < nObjects; ++i) {
      const bool hit = std::fabs(r.displacement[3*i] - 2.) < 1.e-9 &&
                       r.displacement[3*i+1] == 0. && r.displacement[3*i+2] == 0.;
      exact += hit ? 1 : 0;
      if (near[3*i] + 2. > 80.) { ++across; acrossExact += hit ? 1 : 0; }
    }
    check(across > 0, "some copies lie across the face");
    check(exact > nObjects / 2, "more than half the tracers find their own copy");
    check(acrossExact > across / 2, "and more than half of those whose copy is across the face");
  }

  group("the Lagrangian positions are wrapped into the box; the mps is that of the box volume");
  {
    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 3;
    c.seed = 5;
    const Result r = reconstructBox(tracers, {}, c, cube, true);
    check(r.mps == maths::pow(80. * 80. * 80. / (double)nObjects, 1. / 3.),
          "mps = (V/N)^(1/3) of the limits, with every tracer");
    bool wrapped = inside(r.lagrangian, cube);
    for (std::size_t k = 0; k < 3 * nObjects; ++k)
      wrapped = wrapped && r.lagrangian[k] == internal::wrap_coordinate(
                             r.tracers[k] + r.meanDisplacement[k], 0., 80.);
    check(wrapped, "lagrangian is tracers + meanDisplacement wrapped into the box");

    Result edited = r;
    edited.valid[3] = 0;
    recomputeMeans(edited);
    check(inside(edited.lagrangian, cube), "recomputeMeans wraps it too");
    edited.box.reset();
    check_throws([&] { recomputeMeans(edited); }, "a periodic result without a box is refused");
    edited.box = cube;
    edited.geometry = Geometry::Lightcone;
    check_throws([&] { recomputeMeans(edited); }, "a periodic lightcone result is refused");
  }

  group("the minimum image takes [-L/2, L/2): L/2 maps to -L/2");
  {
    check(internal::minimum_image(40., 80.) == -40. && internal::minimum_image(-40., 80.) == -40. &&
          internal::minimum_image(39.5, 80.) == 39.5 && internal::minimum_image(-70., 80.) == 10.,
          "the convention at and near the half period");
    check(internal::wrap_coordinate(-1.e-17, 0., 80.) == 0. ||
          internal::wrap_coordinate(-1.e-17, 0., 80.) == 80.,
          "a wrap that rounds past a bound lands on it");
  }

  group("neighbourAverage and realSpaceBox in a periodic box: across the faces, minimum-image, wrapped");
  {
    // Two valid tracers either side of the x faces of a cube of side 80: 1 from each face.
    const std::vector<double> pos = {1., 40., 40.,   79., 40., 40.,   40., 10., 70.};
    const std::vector<double> values = {1., 3., 0.};
    const std::vector<unsigned> valid = {1, 1, 0};
    const NeighbourAverage open = neighbourAverage(pos, values, valid, 2., false);
    const NeighbourAverage wrapped = neighbourAverage(pos, values, valid, 2., false, cube, true);
    check(open.nNeighbours[0] == 1 && wrapped.nNeighbours[0] == 2,
          "a tracer sees the one across the face only when periodic");
    const double w = std::exp(-0.5);   // the other tracer, at the minimum-image distance 2 = sigma
    check(std::fabs(wrapped.values[0] - (1. + 3. * w) / (1. + w)) < 1.e-12 &&
          std::fabs(wrapped.values[1] - (3. + w) / (1. + w)) < 1.e-12,
          "and each averages itself with the other at the minimum-image distance");
    check(std::isnan(wrapped.values[2]) && wrapped.nNeighbours[2] == 0, "the middle one has no neighbour");
    const NeighbourAverage same = neighbourAverage(pos, values, valid, 2., false, cube, false);
    check(same.values[0] == open.values[0] && same.nNeighbours[0] == open.nNeighbours[0],
          "periodic false is the first overload");
    std::vector<double> out = pos;
    out[0] = 90.;
    check_throws([&] { neighbourAverage(out, values, valid, 2., false, cube, true); },
                 "a position outside the periodic box is refused");
    BoxLimits odd = cube;
    odd.upper[2] = 120.;
    check_throws([&] { neighbourAverage(pos, values, valid, 2., false, odd, true); },
                 "a box that is not an integer multiple is refused");

    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 2;
    c.seed = 9;
    const Result r = reconstructBox(tracers, {}, c, cube, true);
    const DistanceTable table(0.3, 0.7, -1., 0., 0., 2., 21);
    CorrectionConfig quiet;
    quiet.verbosity = Verbosity::Silent;
    const RealSpaceCatalog moved = realSpaceBox(r, 0, 1., table, 1.5, 10., quiet);
    check(inside(moved.cartesian, cube), "every corrected coordinate lies inside the periodic box");
    std::size_t corrected = 0;
    for (const CorrectionStatus s : moved.status) corrected += s == CorrectionStatus::Corrected;
    check(corrected == nObjects, "every tracer is corrected");
    Result edited = r;
    edited.tracers[0] = -5.;
    check_throws([&] { realSpaceBox(edited, 0, 1., table, 1.5, 10., quiet); },
                 "a periodic result whose tracer lies outside its box is refused");
  }

  group("a seed pins a periodic result down at one and four threads, to the recorded bits");
  {
    Config c;
    c.verbosity = Verbosity::Silent;
    c.nRealizations = 3;
    c.seed = 2468;
    std::vector<double> randoms;
    for (unsigned r = 0; r < 4; ++r) {
      const std::vector<double> block = lattice(nSide, spacing, 1. + 2. * r, 3., 50 + r);
      randoms.insert(randoms.end(), block.begin(), block.end());
    }
    randoms[0] -= 80.;
    BoxLimits shifted = cube;
    for (int k = 0; k < 3; ++k) { shifted.lower[k] += 4.5; shifted.upper[k] += 4.5; }

    for (const int threads : {1, 4}) {
      const std::string at = " (" + std::to_string(threads) + " threads)";
      const Result drawn = at_threads(threads, [&] {
        return reconstructBox(tracers, {}, c, shifted, true);
      });
      const Result given = at_threads(threads, [&] {
        return reconstructBox(tracers, randoms, spacing, c, cube, true);
      });
      check_fingerprint(drawn.displacement, 0xf6714359252f51f0ull, "drawn randoms, displacement" + at);
      check_fingerprint(drawn.lagrangian, 0xe9556796d9e9d7acull, "drawn randoms, lagrangian" + at);
      check_fingerprint(drawn.tracers, 0x8de8bfd02ab4928bull, "drawn randoms, wrapped tracers" + at);
      check_fingerprint(given.displacement, 0x6a79ef17ede9e909ull, "given randoms, displacement" + at);
      check_fingerprint(given.meanDisplacement, 0x11519091b08faa46ull, "given randoms, mean displacement" + at);
    }
  }

  return report("test_periodic");
}
