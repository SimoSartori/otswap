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
 *  @file tests/test_detmath.cpp
 *
 *  @brief The deterministic elementary functions of src/detmath.h and the
 *  pixel lookup of src/arc.h, at their edge cases, against reference
 *  values computed at high precision (tests/detmath_reference.h, written
 *  by tools/determinism/edge_values.py).
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>

#include <healpix_base.h>
#include <pointing.h>
#include <vec3.h>

#include "arc.h"
#include "detmath.h"
#include "internal.h"
#include "check.h"
#include "detmath_reference.h"

using namespace otswap;

namespace {

  std::string hex (const double x)
  {
    char text[64];
    std::snprintf(text, sizeof text, "%a", x);
    return text;
  }

  // The doubles in increasing order mapped onto consecutive integers, so
  // that the difference of two is their distance in units in the last
  // place; -0 and +0 both map to 0.
  std::int64_t ordered (const double x)
  {
    std::int64_t bits;
    std::memcpy(&bits, &x, sizeof bits);
    return bits < 0 ? -(bits & 0x7fffffffffffffffLL) : bits;
  }

  bool same_bits (const double a, const double b)
  {
    return std::memcmp(&a, &b, sizeof a) == 0;
  }

  // Checks every case of a table, and returns how many matched the
  // correctly rounded value exactly (any NaN where NaN is expected).
  template <std::size_t N, typename F>
  std::size_t check_table (const std::string& name, const detref::Case (&cases)[N], F f,
                           const bool binary)
  {
    std::size_t exactMatches = 0;
    for (const detref::Case& c : cases) {
      const double got = f(c.x, c.y);
      const std::string call = name + "(" + hex(c.x) + (binary ? ", " + hex(c.y) : "") +
                               ") = " + hex(got) + ", expected " + hex(c.expected);

      if (std::isnan(c.expected)) {
        check(std::isnan(got), call);
        if (std::isnan(got)) ++exactMatches;
        continue;
      }
      if (c.special || c.exact) {
        check(same_bits(got, c.expected), call + ", exactly");
        if (same_bits(got, c.expected)) ++exactMatches;
        continue;
      }
      const std::int64_t d = ordered(got) - ordered(c.expected);
      check(!std::isnan(got) && d >= -1 && d <= 1, call + ", within one ulp");
      if (d == 0) ++exactMatches;
    }
    return exactMatches;
  }

  constexpr double kPi = 3.14159265358979323846;

}

int main ()
{
  group("sin and cos: signed zeros, NaN and infinities, near multiples of pi/2, "
        "and at large arguments, within one ulp of the correctly rounded value");
  {
    std::size_t n = 0, total = 0;
    n += check_table("det_sin", detref::k_sin, [] (double x, double) { return internal::det_sin(x); }, false);
    n += check_table("det_cos", detref::k_cos, [] (double x, double) { return internal::det_cos(x); }, false);
    total = sizeof(detref::k_sin) / sizeof(detref::Case) + sizeof(detref::k_cos) / sizeof(detref::Case);
    std::cout << "    " << n << " of " << total << " correctly rounded" << std::endl;
  }

  group("exp: at the overflow and underflow thresholds, and into the subnormals");
  {
    const std::size_t n = check_table("det_exp", detref::k_exp,
                                      [] (double x, double) { return internal::det_exp(x); }, false);
    std::cout << "    " << n << " of " << sizeof(detref::k_exp) / sizeof(detref::Case)
              << " correctly rounded" << std::endl;
  }

  group("log1p: near 0, near -1, and beyond its domain");
  {
    const std::size_t n = check_table("det_log1p", detref::k_log1p,
                                      [] (double x, double) { return internal::det_log1p(x); }, false);
    std::cout << "    " << n << " of " << sizeof(detref::k_log1p) / sizeof(detref::Case)
              << " correctly rounded" << std::endl;
  }

  group("asin: near 0 and +-1, and beyond its domain");
  {
    const std::size_t n = check_table("det_asin", detref::k_asin,
                                      [] (double x, double) { return internal::det_asin(x); }, false);
    std::cout << "    " << n << " of " << sizeof(detref::k_asin) / sizeof(detref::Case)
              << " correctly rounded" << std::endl;
  }

  group("atan2: signed zeros, infinities, NaN, and all four quadrants");
  {
    const std::size_t n = check_table("det_atan2", detref::k_atan2,
                                      [] (double y, double x) { return internal::det_atan2(y, x); }, true);
    std::cout << "    " << n << " of " << sizeof(detref::k_atan2) / sizeof(detref::Case)
              << " correctly rounded" << std::endl;
  }

  group("pow: the special cases of C99 Annex F, overflow, subnormal results, "
        "and the exponents otswap uses");
  {
    const std::size_t n = check_table("det_pow", detref::k_pow,
                                      [] (double x, double y) { return internal::det_pow(x, y); }, true);
    std::cout << "    " << n << " of " << sizeof(detref::k_pow) / sizeof(detref::Case)
              << " correctly rounded" << std::endl;
  }

  group("vec2pix and ang2pix: on pixel boundaries, at the poles, at RA 0 and 2 pi, "
        "the pixel loc2pix gives with correctly rounded atan2, cos and sin, "
        "on the 32-bit and the 64-bit base");
  {
    for (int k = 0; k < detref::kNNsides; ++k) {
      const int nside = detref::kNsides[k];
      const Healpix_Base ring(nside, RING, SET_NSIDE);
      const Healpix_Base nest(nside, NEST, SET_NSIDE);
      const Healpix_Base2 ring64(nside, RING, SET_NSIDE);
      const Healpix_Base2 nest64(nside, NEST, SET_NSIDE);

      for (const detref::VecCase& c : detref::kVectors) {
        const vec3 v(c.x, c.y, c.z);
        const std::string at = std::string(c.label) + ", NSIDE " + std::to_string(nside);
        check(internal::vec2pix(ring, v) == c.ring[k],
              "vec2pix, " + at + " RING: " + std::to_string(internal::vec2pix(ring, v)) +
              ", expected " + std::to_string(c.ring[k]));
        check(internal::vec2pix(nest, v) == c.nest[k],
              "vec2pix, " + at + " NESTED: " + std::to_string(internal::vec2pix(nest, v)) +
              ", expected " + std::to_string(c.nest[k]));
        check(internal::vec2pix(ring64, v) == c.ring[k],
              "vec2pix, 64-bit, " + at + " RING: " + std::to_string(internal::vec2pix(ring64, v)) +
              ", expected " + std::to_string(c.ring[k]));
        check(internal::vec2pix(nest64, v) == c.nest[k],
              "vec2pix, 64-bit, " + at + " NESTED: " + std::to_string(internal::vec2pix(nest64, v)) +
              ", expected " + std::to_string(c.nest[k]));
      }

      for (const detref::AngCase& c : detref::kAngles) {
        const std::string at = std::string(c.label) + ", NSIDE " + std::to_string(nside);
        check(internal::ang2pix(ring, c.theta, c.phi) == c.ring[k],
              "ang2pix, " + at + " RING: " + std::to_string(internal::ang2pix(ring, c.theta, c.phi)) +
              ", expected " + std::to_string(c.ring[k]));
        check(internal::ang2pix(nest, c.theta, c.phi) == c.nest[k],
              "ang2pix, " + at + " NESTED: " + std::to_string(internal::ang2pix(nest, c.theta, c.phi)) +
              ", expected " + std::to_string(c.nest[k]));
        check(internal::ang2pix(ring64, c.theta, c.phi) == c.ring[k],
              "ang2pix, 64-bit, " + at + " RING: " +
              std::to_string(internal::ang2pix(ring64, c.theta, c.phi)) +
              ", expected " + std::to_string(c.ring[k]));
        check(internal::ang2pix(nest64, c.theta, c.phi) == c.nest[k],
              "ang2pix, 64-bit, " + at + " NESTED: " +
              std::to_string(internal::ang2pix(nest64, c.theta, c.phi)) +
              ", expected " + std::to_string(c.nest[k]));
      }
    }
  }

  group("the 64-bit lookup at NSIDE 2^20 and 2^29, on the same points, with no map: "
        "the pixel loc2pix gives in exact integers");
  {
    for (int k = 0; k < detref::kNLargeNsides; ++k) {
      const std::int64_t nside = detref::kLargeNsides[k];
      const Healpix_Base2 ring(nside, RING, SET_NSIDE);
      const Healpix_Base2 nest(nside, NEST, SET_NSIDE);

      for (const detref::LargeVecCase& c : detref::kLargeVectors) {
        const vec3 v(c.x, c.y, c.z);
        const std::string at = std::string(c.label) + ", NSIDE " + std::to_string(nside);
        check(internal::vec2pix(ring, v) == c.ring[k],
              "vec2pix, " + at + " RING: " + std::to_string(internal::vec2pix(ring, v)) +
              ", expected " + std::to_string(c.ring[k]));
        check(internal::vec2pix(nest, v) == c.nest[k],
              "vec2pix, " + at + " NESTED: " + std::to_string(internal::vec2pix(nest, v)) +
              ", expected " + std::to_string(c.nest[k]));
      }

      for (const detref::LargeAngCase& c : detref::kLargeAngles) {
        const std::string at = std::string(c.label) + ", NSIDE " + std::to_string(nside);
        check(internal::ang2pix(ring, c.theta, c.phi) == c.ring[k],
              "ang2pix, " + at + " RING: " + std::to_string(internal::ang2pix(ring, c.theta, c.phi)) +
              ", expected " + std::to_string(c.ring[k]));
        check(internal::ang2pix(nest, c.theta, c.phi) == c.nest[k],
              "ang2pix, " + at + " NESTED: " + std::to_string(internal::ang2pix(nest, c.theta, c.phi)) +
              ", expected " + std::to_string(c.nest[k]));
      }
    }
  }

  group("up to NSIDE 8192 the 32-bit and the 64-bit base give the same pixel");
  {
    std::mt19937 rng(11);
    std::size_t differ = 0, total = 0;
    for (const int nside : {1, 64, 1024, 8192})
      for (const Healpix_Ordering_Scheme scheme : {RING, NEST}) {
        const Healpix_Base base(nside, scheme, SET_NSIDE);
        const Healpix_Base2 base64(nside, scheme, SET_NSIDE);
        for (int i = 0; i < 20000; ++i) {
          const double x = internal::uniform_real(rng, -1., 1.);
          const double y = internal::uniform_real(rng, -1., 1.);
          const double z = internal::uniform_real(rng, -1., 1.);
          if (x*x + y*y + z*z < 1.e-6) continue;
          const vec3 v(x, y, z);
          if ((std::int64_t)internal::vec2pix(base, v) != internal::vec2pix(base64, v)) ++differ;

          const double theta = internal::uniform_real(rng, 0., kPi);
          const double phi = internal::uniform_real(rng, 0., 2. * kPi);
          if ((std::int64_t)internal::ang2pix(base, theta, phi) !=
              internal::ang2pix(base64, theta, phi)) ++differ;
          total += 2;
        }
        for (const detref::VecCase& c : detref::kVectors) {
          const vec3 v(c.x, c.y, c.z);
          if ((std::int64_t)internal::vec2pix(base, v) != internal::vec2pix(base64, v)) ++differ;
          ++total;
        }
        for (const detref::AngCase& c : detref::kAngles) {
          if ((std::int64_t)internal::ang2pix(base, c.theta, c.phi) !=
              internal::ang2pix(base64, c.theta, c.phi)) ++differ;
          ++total;
        }
      }
    check(differ == 0, std::to_string(differ) + " of " + std::to_string(total) +
                       " directions land in a different pixel");
  }

  group("away from pixel boundaries, the pixel is the one Healpix's own vec2pix and ang2pix give, "
        "on the 32-bit base and, up to NSIDE 2^29, on the 64-bit one");
  {
    std::mt19937 rng(7);
    std::size_t differ = 0, total = 0;
    for (const int nside : {1, 64, 1024, 8192})
      for (const Healpix_Ordering_Scheme scheme : {RING, NEST}) {
        const Healpix_Base base(nside, scheme, SET_NSIDE);
        for (int i = 0; i < 20000; ++i) {
          const double x = internal::uniform_real(rng, -1., 1.);
          const double y = internal::uniform_real(rng, -1., 1.);
          const double z = internal::uniform_real(rng, -1., 1.);
          if (x*x + y*y + z*z < 1.e-6) continue;
          const vec3 v(x, y, z);
          if (internal::vec2pix(base, v) != base.vec2pix(v)) ++differ;

          const double theta = internal::uniform_real(rng, 0., kPi);
          const double phi = internal::uniform_real(rng, 0., 2. * kPi);
          if (internal::ang2pix(base, theta, phi) != base.ang2pix(pointing(theta, phi))) ++differ;
          total += 2;
        }
      }
    for (const std::int64_t nside : {std::int64_t(1) << 15, std::int64_t(1) << 20,
                                     std::int64_t(1) << 29})
      for (const Healpix_Ordering_Scheme scheme : {RING, NEST}) {
        const Healpix_Base2 base(nside, scheme, SET_NSIDE);
        for (int i = 0; i < 20000; ++i) {
          const double x = internal::uniform_real(rng, -1., 1.);
          const double y = internal::uniform_real(rng, -1., 1.);
          const double z = internal::uniform_real(rng, -1., 1.);
          if (x*x + y*y + z*z < 1.e-6) continue;
          const vec3 v(x, y, z);
          if (internal::vec2pix(base, v) != base.vec2pix(v)) ++differ;

          const double theta = internal::uniform_real(rng, 0., kPi);
          const double phi = internal::uniform_real(rng, 0., 2. * kPi);
          if (internal::ang2pix(base, theta, phi) != base.ang2pix(pointing(theta, phi))) ++differ;
          total += 2;
        }
      }
    check(differ == 0, std::to_string(differ) + " of " + std::to_string(total) +
                       " random directions land in a different pixel");
  }

  return report("test_detmath");
}
