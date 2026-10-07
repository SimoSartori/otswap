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
 *  @file tests/test_pixel.cpp
 *
 *  @brief otswap's Healpix geometry and neighbours, against Healpix's own
 *  64-bit base.
 */

#include <array>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include <healpix_base.h>

#include "internal.h"
#include "pixel.h"
#include "check.h"

using namespace otswap;

namespace {

  using Base64 = T_Healpix_Base<std::int64_t>;

  internal::HealpixBase base_of (const std::int64_t nside, const bool nested)
  {
    internal::HealpixBase base;
    base.set(nside, nested);
    return base;
  }

  // Pixels whose eight neighbours differ from Healpix's.
  std::size_t neighbour_disagreements (const internal::HealpixBase& ours, const Base64& healpix,
                                       const std::vector<std::int64_t>& pixels)
  {
    std::size_t differ = 0;
    std::array<std::int64_t, 8> mine;
    fix_arr<std::int64_t, 8> theirs;
    for (const std::int64_t p : pixels) {
      ours.neighbors(p, mine);
      healpix.neighbors(p, theirs);
      for (std::size_t m = 0; m < 8; ++m)
        if (mine[m] != theirs[m]) { ++differ; break; }
    }
    return differ;
  }

  // On every face, the pixels along its four edges at a few positions,
  // its corners included, and the ones next to them inside the face; then
  // count random pixels anywhere.
  std::vector<std::int64_t> edge_and_random_pixels (const Base64& healpix, const std::size_t count,
                                                    std::mt19937_64& rng)
  {
    const int n = (int)healpix.Nside();
    const std::vector<int> along = {0, 1, 2, n / 2, n - 3, n - 2, n - 1};
    std::vector<std::int64_t> pixels;
    for (int face = 0; face < 12; ++face)
      for (const int t : along)
        for (const int edge : {0, 1, n - 2, n - 1}) {
          pixels.push_back(healpix.xyf2pix(edge, t, face));
          pixels.push_back(healpix.xyf2pix(t, edge, face));
        }
    for (std::size_t i = 0; i < count; ++i)
      pixels.push_back((std::int64_t)(rng() % (std::uint64_t)healpix.Npix()));
    return pixels;
  }

}

int main ()
{
  group("the geometry: NSIDE, the pixel count and the scheme, from NSIDE 1 to 2^29");
  {
    for (const std::int64_t nside : {std::int64_t(1), std::int64_t(2), std::int64_t(3),
                                     std::int64_t(48), std::int64_t(8192),
                                     std::int64_t(1) << 20, std::int64_t(1) << 29})
      for (const bool nested : {false, true}) {
        if (nested && (nside & (nside - 1)) != 0) continue;
        const internal::HealpixBase ours = base_of(nside, nested);
        const Base64 healpix(nside, nested ? NEST : RING, SET_NSIDE);
        const std::string what = "NSIDE " + std::to_string(nside) + (nested ? " NESTED" : " RING");
        check(ours.nside() == healpix.Nside(), what + ": NSIDE");
        check(ours.npix() == healpix.Npix() && ours.npix() == 12 * nside * nside,
              what + ": 12 NSIDE^2 pixels");
        check(ours.nested() == nested, what + ": the scheme");
      }
  }

  group("the neighbours of every pixel are Healpix's, RING and NESTED from NSIDE 1 to 64, and "
        "RING at NSIDE 3, 5 and 48; face edges and corners included");
  {
    for (const int nside : {1, 2, 3, 4, 5, 8, 16, 32, 48, 64})
      for (const bool nested : {false, true}) {
        if (nested && (nside & (nside - 1)) != 0) continue;
        const internal::HealpixBase ours = base_of(nside, nested);
        const Base64 healpix(nside, nested ? NEST : RING, SET_NSIDE);
        std::vector<std::int64_t> pixels((std::size_t)healpix.Npix());
        for (std::size_t p = 0; p < pixels.size(); ++p) pixels[p] = (std::int64_t)p;
        const std::size_t differ = neighbour_disagreements(ours, healpix, pixels);
        check(differ == 0, "NSIDE " + std::to_string(nside) + (nested ? " NESTED: " : " RING: ") +
                           std::to_string(differ) + " of " + std::to_string(pixels.size()) +
                           " pixels have other neighbours");
      }
  }

  group("at NSIDE 2^15, 2^20 and 2^29, with no map, the neighbours of pixels on every face edge "
        "and corner and of random pixels are Healpix's, RING and NESTED");
  {
    std::mt19937_64 rng(29);
    for (const std::int64_t nside : {std::int64_t(1) << 15, std::int64_t(1) << 20,
                                     std::int64_t(1) << 29})
      for (const bool nested : {false, true}) {
        const internal::HealpixBase ours = base_of(nside, nested);
        const Base64 healpix(nside, nested ? NEST : RING, SET_NSIDE);
        const std::vector<std::int64_t> pixels = edge_and_random_pixels(healpix, 20000, rng);
        const std::size_t differ = neighbour_disagreements(ours, healpix, pixels);
        check(differ == 0, "NSIDE " + std::to_string(nside) + (nested ? " NESTED: " : " RING: ") +
                           std::to_string(differ) + " of " + std::to_string(pixels.size()) +
                           " pixels have other neighbours");
      }
  }

  return report("test_pixel");
}
