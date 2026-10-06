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
 *  @file src/arc.h
 *
 *  @brief Pixel lookup, and the search for unobserved pixels along a great
 *  circle arc. Internal, and the only internal header that needs Healpix.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_ARC_H
#define OTSWAP_ARC_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <healpix_base.h>
#include <vec3.h>

namespace otswap {

  namespace internal {

    /// The largest NSIDE a mask may have, 2^29: the limit of Healpix's
    /// 64-bit pixel index.
    constexpr std::int64_t kMaxNside = std::int64_t(1) << 29;

    /// A full-sky mask: the pixel geometry, and one byte per pixel, 1 when
    /// the pixel is observed.
    struct PixelMask {
      T_Healpix_Base<std::int64_t> base;
      std::vector<std::uint8_t> observed;
      std::int64_t allowed = 0;  ///< the number of observed pixels
    };

    /// Checks @p nside for a mask: in [1, kMaxNside], and a power of 2 when
    /// @p nested. @exception Error naming @p source otherwise.
    void check_nside (std::int64_t nside, bool nested, const std::string& source);

    /// Sets the geometry of @p mask, checked by check_nside, and sizes its
    /// pixels, all unobserved.
    void set_geometry (PixelMask& mask, std::int64_t nside, bool nested,
                       const std::string& source);

    /// Marks @p count pixels from @p first by the values at @p values: a
    /// pixel is observed when its value exceeds kMaskAllowedAbove, so NaN,
    /// -inf and Healpix's UNSEEN, -1.6375e30, are unobserved and +inf is
    /// observed. Updates the count of observed pixels.
    void mark_observed (PixelMask& mask, std::int64_t first, const double* values,
                        std::size_t count);

    /// Pixel holding the direction @p v, which need not be normalised.
    /// Healpix's T_Healpix_Base::vec2pix line for line, over its own
    /// loc2pix, with the atan2 of detmath.h in place of the system one, so
    /// that the pixel is the same on every platform. v must be finite and
    /// non-zero. Defined for the 32-bit (int) and 64-bit (std::int64_t)
    /// bases, which give the same pixel wherever the 32-bit one holds the
    /// pixel count.
    template <typename I>
    I vec2pix (const T_Healpix_Base<I>& base, const vec3& v);

    /// Pixel holding colatitude @p theta and longitude @p phi, radians.
    /// Healpix's T_Healpix_Base::ang2pix line for line, with the cos and
    /// sin of detmath.h, and without its assertion: theta must lie in
    /// [0, pi] and phi be finite, which the caller checks. Defined for the
    /// same two bases as vec2pix.
    template <typename I>
    I ang2pix (const T_Healpix_Base<I>& base, double theta, double phi);

    /// True when the pixel is observed.
    inline bool pixel_observed (const PixelMask& mask, const std::int64_t pixel)
    {
      return mask.observed[(std::size_t)pixel] != 0;
    }

    /// Distinct unobserved pixels met along the great circle arc from @p a
    /// to @p b, by the adaptive bisection rejectMaskCrossings documents.
    /// The pixels holding a and b are exempt.
    ///
    /// @param a,b unit vectors, not exactly opposite
    /// @param limit the search stops as soon as more than this many are found
    /// @param found receives the pixels, in the order met; cleared first
    /// @return found.size(), which exceeds @p limit only when the search
    ///   stopped early
    std::size_t arc_unobserved_pixels (const PixelMask& mask,
                                       const vec3& a, const vec3& b,
                                       unsigned limit, std::vector<std::int64_t>& found);

  }

}

#endif
