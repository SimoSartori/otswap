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
 *  @file tests/fixtures_v010.h
 *
 *  @brief The inputs of the otswap 0.1.0 files in tests/data: written by
 *  tests/data/write_fixtures_v010.cpp built against otswap 0.1.0, and
 *  rebuilt by the tests, which reconstruct them again and compare. Every
 *  value is an integer times a constant, the same double in both
 *  programs.
 */

#ifndef OTSWAP_TESTS_FIXTURES_V010_H
#define OTSWAP_TESTS_FIXTURES_V010_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "otswap/OT.h"

namespace fixtures_v010 {

  /// Objects of each file, and realizations.
  constexpr std::size_t kObjects = 500;
  constexpr unsigned kRealizations = 2;
  constexpr unsigned kSeed = 2026;

  /// Three integers of point i, each spread over [0, 1000) without
  /// repeating within the first thousand points.
  inline void spread (const std::size_t i, const std::size_t salt, long& a, long& b, long& c)
  {
    a = (long)((i * 7919 + salt * 31) % 1000);
    b = (long)((i * 104729 + salt * 17) % 997);
    c = (long)((i * 1299709 + salt * 13) % 991);
  }

  /// Box tracers in [0, 100) Mpc/h on each axis: three 20-bit integers of
  /// point i from splitmix64(i), times 100 / 2^20. No two neighbours of a
  /// tracer lie at exactly the same distance from it, so the result does
  /// not depend on how ties were ordered.
  inline std::vector<double> box_tracers ()
  {
    std::vector<double> t(3 * kObjects);
    for (std::size_t i = 0; i < kObjects; ++i) {
      std::uint64_t z = (std::uint64_t)i + 0x9E3779B97F4A7C15ULL;
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
      z ^= z >> 31;
      t[3*i] = (double)(z >> 44) * (100. / 1048576.);
      t[3*i+1] = (double)((z >> 24) & 0xFFFFF) * (100. / 1048576.);
      t[3*i+2] = (double)((z >> 4) & 0xFFFFF) * (100. / 1048576.);
    }
    return t;
  }

  /// Sky points, radians: right ascension in [0.2, 0.7), declination in
  /// [-0.25, 0.25), redshift in [0.2, 0.7).
  inline std::vector<double> sky (const std::size_t n, const std::size_t salt)
  {
    std::vector<double> s(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
      long a, b, c;
      spread(i, salt, a, b, c);
      s[3*i] = (double)(400 + a) * 0.0005;
      s[3*i+1] = (double)(b - 500) * 0.0005;
      s[3*i+2] = (double)(400 + c) * 0.0005;
    }
    return s;
  }

  inline std::vector<double> lightcone_tracers () { return sky(kObjects, 0); }
  inline std::vector<double> lightcone_randoms () { return sky(5 * kObjects, 1); }

  /// NSIDE 16, RING, every ninth pixel unobserved.
  inline otswap::Mask mask ()
  {
    std::vector<double> pixels(12 * 16 * 16, 1.);
    for (std::size_t p = 0; p < pixels.size(); p += 9) pixels[p] = 0.;
    return otswap::Mask(pixels, otswap::PixelOrdering::Ring);
  }

  inline otswap::DistanceTable table () { return otswap::DistanceTable(0.3, 0.7, -1., 0., 0., 1.5, 4000); }

  inline otswap::RedshiftCut cut () { return otswap::RedshiftCut{0.25, 0.65}; }

  inline otswap::Config config ()
  {
    otswap::Config c;
    c.nRealizations = kRealizations;
    c.seed = kSeed;
    c.verbosity = otswap::Verbosity::Silent;
    return c;
  }

}

#endif
