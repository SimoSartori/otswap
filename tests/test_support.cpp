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
 *  @file tests/test_support.cpp
 *
 *  @brief The owned random-number routines, the input validation and the
 *  coordinate conversion.
 */

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "internal.h"
#include "check.h"

using namespace otswap;

namespace {

  // internal::stream(12345, 7, Stream::Swap)(), recorded once. A different
  // value on some toolchain means results are not portable there.
  constexpr std::uint32_t kFrozenStreamDraw = 1462087691u;

}

int main ()
{
  // The owned draws must be frozen: the whole point of not using the std::
  // distributions is that a library or compiler change cannot alter a
  // result silently. These expected values are the contract.
  group("uniform_int and shuffle are frozen against a fixed seed");
  {
    std::mt19937 rng(12345);
    const unsigned expected[8] = {
      internal::uniform_int(rng, 100), internal::uniform_int(rng, 100),
      internal::uniform_int(rng, 100), internal::uniform_int(rng, 100),
      internal::uniform_int(rng, 100), internal::uniform_int(rng, 100),
      internal::uniform_int(rng, 100), internal::uniform_int(rng, 100)
    };

    std::mt19937 again(12345);
    for (int i = 0; i < 8; ++i)
      check(internal::uniform_int(again, 100) == expected[i],
            "uniform_int repeats for a given seed");

    std::vector<unsigned> a(10), b(10);
    std::iota(a.begin(), a.end(), 0u);
    std::iota(b.begin(), b.end(), 0u);
    std::mt19937 r1(777), r2(777);
    internal::shuffle(a, r1);
    internal::shuffle(b, r2);
    check(a == b, "shuffle repeats for a given seed");

    std::vector<unsigned> sorted = a;
    std::sort(sorted.begin(), sorted.end());
    for (unsigned i = 0; i < 10; ++i)
      check(sorted[i] == i, "shuffle returns a permutation");
  }

  group("uniform_int stays inside its bound and covers it");
  {
    std::mt19937 rng(1);
    std::vector<int> seen(7, 0);
    for (int i = 0; i < 20000; ++i) {
      const unsigned v = internal::uniform_int(rng, 7);
      check(v < 7, "uniform_int(7) stays below 7");
      ++seen[v];
    }
    for (int i = 0; i < 7; ++i)
      check(seen[i] > 2000, "uniform_int(7) reaches every value");

    check(internal::uniform_int(rng, 1) == 0, "uniform_int(1) is always zero");
    check(internal::uniform_int(rng, 0) == 0, "uniform_int(0) is zero");
  }

  group("uniform_real stays inside its interval");
  {
    std::mt19937 rng(42);
    for (int i = 0; i < 20000; ++i) {
      const double v = internal::uniform_real(rng, -3., 5.);
      check(v >= -3. && v < 5., "uniform_real stays in [lo, hi)");
    }
  }

  group("the random streams are fixed by (seed, realization, stream) and are distinct");
  {
    using internal::Stream;

    // The first output of one stream, frozen: std::seed_seq and
    // std::mt19937 are specified exactly, so every toolchain must agree.
    std::mt19937 frozen = internal::stream(12345, 7, Stream::Swap);
    check(frozen() == kFrozenStreamDraw, "a stream's first draw is the same on every toolchain");

    std::mt19937 again = internal::stream(12345, 7, Stream::Swap);
    std::mt19937 same = internal::stream(12345, 7, Stream::Swap);
    for (int i = 0; i < 100; ++i)
      check(again() == same(), "a stream repeats exactly");

    // Distinct triples give distinct sequences, including the triples an
    // additive offset would have made coincide: seed + realization * stride
    // is equal for (s + stride, 0) and (s, 1).
    const Stream ids[4] = {Stream::Supply, Stream::Randoms, Stream::Seeding, Stream::Swap};
    std::set<std::vector<std::uint32_t>> seen;
    std::size_t made = 0;
    for (const unsigned seed : {1u, 2u, 100u, 100u + 19937u, 100u + 21317u, 4294967295u})
      for (const unsigned rec : {0u, 1u, 2u, 1000u})
        for (const Stream id : ids) {
          std::mt19937 e = internal::stream(seed, rec, id);
          std::vector<std::uint32_t> head(4);
          for (auto& h : head) h = (std::uint32_t)e();
          seen.insert(head);
          ++made;
        }
    check(seen.size() == made, "no two of the streams tried share their opening draws");
  }

  group("malformed coordinate arrays are rejected");
  {
    check_throws([] { internal::check_coordinates({}, "a"); }, "an empty array raises");
    check_throws([] { internal::check_coordinates({1., 2.}, "a"); },
                 "a size that is not a multiple of three raises");
    check_throws([] {
      internal::check_coordinates({1., 2., std::nan("")}, "a");
    }, "a non-finite entry raises");

    check(internal::check_coordinates({}, "a", true) == 0,
          "an empty array is accepted where the header allows it");
    check(internal::check_coordinates({1., 2., 3., 4., 5., 6.}, "a") == 2,
          "the object count is the size over three");
  }

  group("Config is validated and seed 0 is resolved");
  {
    Config bad;
    bad.nRealizations = 0;
    check_throws([&] { internal::check_config(bad); }, "zero realizations raise");

    Config c;
    c.convergence = 0.;
    check_throws([&] { internal::check_config(c); }, "a non-positive convergence raises");
    c.convergence = 1.;
    check_throws([&] { internal::check_config(c); }, "a convergence of one raises");

    Config d;
    d.cellSize = 0.;
    check_throws([&] { internal::check_config(d); }, "a non-positive cell size raises");

    Config e;
    e.seed = 4242;
    check(internal::check_config(e) == 4242u, "a given seed is kept");

    Config f;
    f.seed = 0;
    check(internal::check_config(f) != 0u, "seed 0 resolves to a drawn, non-zero seed");
  }

  group("a random supply too small for the realizations is refused");
  {
    check_throws([] { internal::check_random_supply(10, 10, 2); },
                 "fewer randoms than nRealizations * nObjects raises");
    internal::check_random_supply(20, 10, 2);
  }

  group("the sky and Cartesian conversions invert each other");
  {
    const double cases[4][3] = {
      {0.3, 0.2, 1000.}, {6.0, -0.9, 250.}, {3.14159, 0., 1.}, {0., 1.5, 4321.}
    };

    for (const auto& c : cases) {
      double x, y, z, ra, dec, d;
      internal::to_cartesian(c[0], c[1], c[2], x, y, z);
      internal::to_sky(x, y, z, ra, dec, d);
      check_close(ra, internal::normalize_ra(c[0]), 1.e-12, "right ascension round-trips");
      check_close(dec, c[1], 1.e-12, "declination round-trips");
      check_close(d, c[2], 1.e-9, "distance round-trips");
    }

    double ra, dec, d;
    internal::to_sky(0., 0., 0., ra, dec, d);
    check(d == 0. && ra == 0. && dec == 0., "the origin is defined, not a NaN");
  }

  group("a piecewise linear profile is extrapolated linearly beyond its ends");
  {
    const std::vector<double> x {1., 2., 4.};
    const std::vector<double> y {10., 20., 30.};
    check_close(internal::profile_at(x, y, 0.), 0., 1.e-12,
                "below the first node, along the first segment");
    check_close(internal::profile_at(x, y, 6.), 40., 1.e-12,
                "above the last node, along the last segment");
    check_close(internal::profile_at(x, y, 1.5), 15., 1.e-12, "between two nodes");
    check(internal::profile_at(x, y, 2.) == 20., "on a node, exactly the node's value");
    check(internal::profile_at(x, y, 4.) == 30., "on the last node, exactly its value");
    check_close(internal::profile_at({5.}, {7.}, 99.), 7., 1.e-12, "a single node is constant");
  }

  group("mps(z) is extrapolated linearly, and a non-positive value raises");
  {
    MpsProfile profile;
    profile.redshift = {0.5, 0.6, 0.7};
    profile.mps = {10., 12., 20.};
    profile.count = {300, 300, 300};

    check_close(profile.at(0.45), 9., 1.e-9, "below the first bin centre");
    check_close(profile.at(0.75), 24., 1.e-9, "above the last bin centre");
    check_close(profile.at(0.65), 16., 1.e-9, "between bin centres");

    // Along the first segment mps falls by 2 per 0.1 in z, so it reaches
    // zero at z = 0 and is negative below.
    check_throws([&] { profile.at(-0.01); },
                 "an extrapolation to a negative mps raises");
    check_throws([&] { profile.at(0.); },
                 "an extrapolation to exactly zero raises");

    try {
      profile.at(-0.05);
    }
    catch (const Error& e) {
      const std::string message = e.what();
      check(message.find("-0.05") != std::string::npos, "the message names the redshift");
      check(message.find("extrapolated") != std::string::npos,
            "and says the value came from an extrapolation");
    }
  }

  group("an Error records where it was constructed, and its message stays the message alone");
  {
    const int line = __LINE__ + 1;
    const Error here("a message");
    check(std::string(here.what()) == "a message", "what() holds the message alone");
    check(std::string(here.file()).find("test_support.cpp") != std::string::npos,
          "file() names this file");
    check(std::string(here.function()) == "main", "function() names this function");
    check(here.line() == line, "line() is the line of the construction");

    const Error elsewhere("message", "some/file.cpp", "f", 7);
    check(std::string(elsewhere.file()) == "some/file.cpp" && std::string(elsewhere.function()) == "f" &&
          elsewhere.line() == 7, "an explicit location is kept");

    try {
      internal::check_random_supply(1, 10, 1);
      check(false, "a supply too small raises");
    }
    catch (const Error& e) {
      check(std::string(e.file()).find("Support.cpp") != std::string::npos &&
            std::string(e.function()) == "check_random_supply" && e.line() > 0,
            "a library error carries the location of its throw site");
      check(std::string(e.what()).find("Support.cpp") == std::string::npos,
            "which its message does not mention");
    }
  }

  return report("test_support");
}
