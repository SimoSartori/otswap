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
 *  @file include/otswap/Random.h
 *
 *  @brief Random draws that give the same values on every platform.
 *
 *  std::mt19937 is specified exactly by the C++ standard, but the standard
 *  distributions and std::shuffle are not: their algorithms are left to
 *  each library, so the same engine can give different values with
 *  different compilers. The functions here are written out over the
 *  engine's 32-bit outputs, so a given engine state gives the same values
 *  on every supported platform and compiler. The reconstructions draw
 *  through them; on an axis of zero width a draw takes its two outputs and
 *  gives the bound.
 *
 *  Their algorithms, described below down to the engine outputs each call
 *  consumes, are fixed: no later version of otswap changes them, so a
 *  program that draws with them keeps its values when otswap is updated.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_RANDOM_H
#define OTSWAP_RANDOM_H

#include <cstddef>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "otswap/OT.h"

namespace otswap {

  /// @brief A uniform integer in [0, bound).
  ///
  /// Lemire's multiply-and-reject method: an engine output x gives the
  /// 64-bit product x * bound; when its low 32 bits are below
  /// 2^32 mod bound, x is rejected and a new output drawn. The result is
  /// the high 32 bits of the first product kept. Every value in [0, bound)
  /// is exactly equally likely. A bound of 0 returns 0 and draws nothing.
  ///
  /// @param engine the engine, advanced by one output per attempt
  /// @param bound the number of possible values
  inline unsigned uniformInt (std::mt19937& engine, const unsigned bound)
  {
    if (bound == 0) return 0;

    const std::uint64_t range = bound;
    std::uint64_t product = (std::uint64_t)engine() * range;
    std::uint32_t low = (std::uint32_t)product;

    if (low < range) {
      const std::uint32_t threshold = (std::uint32_t)(0u - (std::uint32_t)range) % (std::uint32_t)range;
      while (low < threshold) {
        product = (std::uint64_t)engine() * range;
        low = (std::uint32_t)product;
      }
    }

    return (unsigned)(product >> 32);
  }

  /// @brief A uniform double in [lo, hi).
  ///
  /// Each call draws exactly two engine outputs, a then b, and forms
  /// u = ((a >> 5) * 2^26 + (b >> 6)) / 2^53, a multiple of 2^-53 in
  /// [0, 1). The result is lo + u * (hi - lo); when rounding makes it hi,
  /// the largest double below hi is returned instead, so the result is
  /// always below hi.
  ///
  /// @param engine the engine, advanced by two outputs
  /// @param lo,hi the interval's bounds: both finite, lo below hi, and
  /// hi - lo finite
  /// @throw Error when a bound is not finite, lo is not below hi, or
  /// hi - lo overflows; the engine is then not advanced
  double uniformReal (std::mt19937& engine, double lo, double hi);

  /// @brief Puts the elements of @p v in a uniformly random order.
  ///
  /// Fisher-Yates, from the end: for i = v.size() down to 2, the element at
  /// i - 1 is swapped with the one at uniformInt(engine, i). The engine is
  /// advanced only by those calls.
  ///
  /// @param v the elements, at most 2^32 - 1 of them
  /// @param engine the engine
  /// @throw Error when @p v holds more than 2^32 - 1 elements; @p v and the
  /// engine are then unchanged
  template <typename T>
  void shuffle (std::vector<T>& v, std::mt19937& engine)
  {
    if (v.size() > (std::size_t)std::numeric_limits<unsigned>::max())
      throw Error("the sequence to shuffle holds " + std::to_string(v.size()) +
                  " elements; it must hold at most " +
                  std::to_string(std::numeric_limits<unsigned>::max()));

    for (std::size_t i = v.size(); i > 1; --i)
      std::swap(v[i-1], v[uniformInt(engine, (unsigned)i)]);
  }

}

#endif
