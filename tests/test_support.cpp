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
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "otswap/Random.h"

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
  group("uniformInt is frozen: literal values for several bounds from a fixed seed");
  {
    std::mt19937 rng(12345);
    const unsigned bounds[6] = {1u, 2u, 7u, 100u, 1000003u, 4294967295u};
    const unsigned expected[6] = {0u, 1u, 2u, 13u, 183919u, 170765736u};
    for (int i = 0; i < 6; ++i)
      check(uniformInt(rng, bounds[i]) == expected[i],
            "uniformInt(" + std::to_string(bounds[i]) + ") gives its recorded value");

    std::mt19937 fresh(12345);
    fresh.discard(6);
    check(rng() == fresh(), "one engine output per accepted draw");
  }

  group("uniformInt draws again when the low half of the product is below 2^32 mod bound");
  {
    const unsigned bound = 3000000000u;
    const std::uint32_t threshold = (std::uint32_t)(0u - bound) % bound;

    std::mt19937 probe(3);
    const std::uint64_t first = (std::uint64_t)probe() * bound;
    const std::uint64_t second = (std::uint64_t)probe() * bound;
    check((std::uint32_t)first < threshold && (std::uint32_t)second >= threshold,
          "seed 3 rejects its first output for this bound and keeps the second");

    std::mt19937 rng(3);
    check(uniformInt(rng, bound) == 212174640u, "the value comes from the second output");

    std::mt19937 fresh(3);
    fresh.discard(2);
    check(rng() == fresh(), "exactly two outputs were drawn");
  }

  group("uniformInt with bound 0 returns 0 and leaves the engine untouched");
  {
    std::mt19937 rng(12345), fresh(12345);
    check(uniformInt(rng, 0) == 0u, "uniformInt(0) is zero");
    check(rng() == fresh(), "the next output is a fresh engine's first");
  }

  group("uniformInt stays inside its bound and covers it");
  {
    std::mt19937 rng(1);
    std::vector<int> seen(7, 0);
    for (int i = 0; i < 20000; ++i) {
      const unsigned v = uniformInt(rng, 7);
      check(v < 7, "uniformInt(7) stays below 7");
      ++seen[v];
    }
    for (int i = 0; i < 7; ++i)
      check(seen[i] > 2000, "uniformInt(7) reaches every value");

    check(uniformInt(rng, 1) == 0, "uniformInt(1) is always zero");
  }

  group("uniformReal is frozen: literal values from a fixed seed, two engine outputs per draw");
  {
    std::mt19937 rng(2026);
    check(uniformReal(rng, 0., 1.) == 0x1.c13848c4407d8p-3, "a draw in [0, 1)");
    check(uniformReal(rng, -3., 5.) == 0x1.376463d841cep-2, "a draw in [-3, 5)");
    check(uniformReal(rng, 100., 100.5) == 0x1.91f40990bd9fap+6, "a draw in [100, 100.5)");
    check(uniformReal(rng, -1.e-3, 0.) == -0x1.ddade73e35309p-11, "a draw in [-0.001, 0)");

    std::mt19937 fresh(2026);
    fresh.discard(8);
    check(rng() == fresh(), "four draws use eight outputs");
  }

  group("uniformReal never returns hi: the largest draw on [1, 2) rounds to 2 and gives the double below 2");
  {
    std::ostringstream words;
    for (int i = 0; i < 624; ++i) words << (i == 397 || i == 398 ? 316513203u : 0u) << ' ';
    words << 624;
    std::istringstream in(words.str());
    std::mt19937 rng;
    in >> rng;

    std::mt19937 copy = rng;
    check(copy() == 4294967295u && copy() == 4294967295u, "the engine's next two outputs are the largest");
    check(1. + (double)(((std::uint64_t)1 << 53) - 1) / 9007199254740992. == 2.,
          "the largest fraction added to 1 rounds to 2");
    check(uniformReal(rng, 1., 2.) == std::nextafter(2., 1.), "the result is the double just below 2");
  }

  group("uniformReal stays inside its interval");
  {
    std::mt19937 rng(42);
    for (int i = 0; i < 20000; ++i) {
      const double v = uniformReal(rng, -3., 5.);
      check(v >= -3. && v < 5., "uniformReal stays in [lo, hi)");
    }
  }

  group("uniformReal refuses a non-finite bound, lo not below hi and a width that overflows, before drawing");
  {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::mt19937 rng(5), fresh(5);

    throws_with([&] { uniformReal(rng, nan, 1.); },
                "the interval of a uniform draw has a bound of nan; both bounds must be finite", "a NaN lo");
    throws_with([&] { uniformReal(rng, 0., nan); },
                "the interval of a uniform draw has a bound of nan; both bounds must be finite", "a NaN hi");
    throws_with([&] { uniformReal(rng, -inf, 1.); },
                "the interval of a uniform draw has a bound of -inf; both bounds must be finite",
                "an infinite lo");
    throws_with([&] { uniformReal(rng, 0., inf); },
                "the interval of a uniform draw has a bound of inf; both bounds must be finite",
                "an infinite hi");
    throws_with([&] { uniformReal(rng, 2., 1.); },
                "the interval of a uniform draw is [2.000000, 1.000000); "
                "the lower bound must be below the upper one", "lo above hi");
    throws_with([&] { uniformReal(rng, 1., 1.); },
                "the interval of a uniform draw is [1.000000, 1.000000); "
                "the lower bound must be below the upper one", "lo equal to hi");
    throws_with([&] { uniformReal(rng, -1.e308, 1.e308); },
                "the interval of a uniform draw is [" + std::to_string(-1.e308) + ", " +
                std::to_string(1.e308) + "); its width, the upper bound minus the lower, must be finite",
                "a width that overflows");

    check(rng() == fresh(), "no refused call advanced the engine");
  }

  group("shuffle is frozen: literal permutations for two sizes and two seeds");
  {
    auto shuffled = [] (const std::size_t n, const unsigned seed) {
      std::vector<unsigned> v(n);
      std::iota(v.begin(), v.end(), 0u);
      std::mt19937 rng(seed);
      shuffle(v, rng);
      return v;
    };
    check(shuffled(5, 7) == std::vector<unsigned>{1, 3, 2, 4, 0}, "5 elements, seed 7");
    check(shuffled(5, 99) == std::vector<unsigned>{0, 4, 1, 2, 3}, "5 elements, seed 99");
    check(shuffled(12, 7) == std::vector<unsigned>{5, 1, 11, 9, 8, 4, 6, 3, 10, 7, 2, 0},
          "12 elements, seed 7");
    check(shuffled(12, 99) == std::vector<unsigned>{2, 1, 11, 9, 3, 0, 5, 6, 10, 4, 7, 8},
          "12 elements, seed 99");

    std::vector<double> one {3.5}, none;
    std::mt19937 rng(7), fresh(7);
    shuffle(one, rng);
    shuffle(none, rng);
    check(one == std::vector<double>{3.5} && rng() == fresh(),
          "one or no element: unchanged, and nothing drawn");
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
    throws_with([&] { internal::check_config(bad); },
                "the number of realizations is 0; it must be at least 1", "zero realizations raise");

    Config c;
    c.convergence = 0.;
    throws_with([&] { internal::check_config(c); },
                "the convergence threshold is 0.000000; it must lie strictly between 0 and 1",
                "a non-positive convergence raises");
    c.convergence = 1.;
    throws_with([&] { internal::check_config(c); },
                "the convergence threshold is 1.000000; it must lie strictly between 0 and 1",
                "a convergence of one raises");
    c.convergence = std::nan("");
    throws_with([&] { internal::check_config(c); },
                "the convergence threshold is nan; it must lie strictly between 0 and 1",
                "a NaN convergence raises");

    Config d;
    d.cellSize = 0.;
    throws_with([&] { internal::check_config(d); },
                "the grid cell size is 0.000000 mean particle separations; it must be positive",
                "a non-positive cell size raises");

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

  group("check_ids: empty or one per object, any 64-bit value, no value twice; the first repeat named");
  {
    auto message = [] (const std::vector<std::int64_t>& ids, const std::size_t n) {
      try {
        internal::check_ids(ids, n);
      }
      catch (const Error& e) {
        return std::string(e.what());
      }
      return std::string();
    };
    const std::int64_t lowest = std::numeric_limits<std::int64_t>::min();
    const std::int64_t highest = std::numeric_limits<std::int64_t>::max();
    check(message({}, 5).empty() && message({}, 0).empty(), "empty IDs are accepted for any count");
    check(message({lowest, -1, 0, highest, 9007199254740993}, 5).empty(),
          "negative and extreme values are accepted");
    check(message({1, 2, 3}, 4) == "ids holds 3 entries for 4 objects; give one per object, or none" &&
          message({1}, 0) == "ids holds 1 entries for 0 objects; give one per object, or none",
          "a wrong length is refused, with both counts");
    check(message({7, 3, 9, 3, 7, 3}, 6) == "ids holds 3 twice, at entries 1 and 3; each object needs its own ID",
          "the repeat with the smallest second entry is named, with the first entry holding its value");
    check(message({lowest, 5, 6, 5, lowest}, 5) ==
            "ids holds 5 twice, at entries 1 and 3; each object needs its own ID" &&
          message({highest, 2, highest}, 3) ==
            "ids holds 9223372036854775807 twice, at entries 0 and 2; each object needs its own ID",
          "the order of the values does not decide which repeat is named");
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

  group("the coordinate checks name the entry and the object from 0");
  {
    std::string m;
    try { internal::check_coordinates({1., 2., 3., 4., std::nan(""), 6.}, "the tracer array"); }
    catch (const Error& e) { m = e.what(); }
    check(m == "the tracer array holds a non-finite value at entry 4 (object 1, component 1), each from 0",
          "a non-finite coordinate: " + m);
    m.clear();
    try { internal::check_sky({0.1, 0.2, 0.3, 0.1, 2., 0.3}, "the sky array"); }
    catch (const Error& e) { m = e.what(); }
    check(m == "the sky array holds a declination of 2.000000 at object 1 (from 0), outside [-pi/2, pi/2]",
          "a declination outside the range: " + m);
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
