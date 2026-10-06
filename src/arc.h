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
#include <vector>

#include <healpix_map.h>
#include <vec3.h>

namespace otswap {

  namespace internal {

    /// Pixel holding the direction @p v, which need not be normalised.
    /// Healpix's T_Healpix_Base::vec2pix line for line, over its own
    /// loc2pix, with the atan2 of detmath.h in place of the system one, so
    /// that the pixel is the same on every platform. v must be finite and
    /// non-zero.
    int vec2pix (const Healpix_Base& base, const vec3& v);

    /// Pixel holding colatitude @p theta and longitude @p phi, radians.
    /// Healpix's T_Healpix_Base::ang2pix line for line, with the cos and
    /// sin of detmath.h, and without its assertion: theta must lie in
    /// [0, pi] and phi be finite, which the caller checks.
    int ang2pix (const Healpix_Base& base, double theta, double phi);

    /// True when the pixel is observed: its value exceeds kMaskAllowedAbove.
    /// NaN and Healpix's UNSEEN, -1.6375e30, are therefore unobserved.
    bool pixel_observed (const Healpix_Map<float>& map, int pixel);

    /// Distinct unobserved pixels met along the great circle arc from @p a
    /// to @p b, by the adaptive bisection rejectMaskCrossings documents.
    /// The pixels holding a and b are exempt.
    ///
    /// @param a,b unit vectors, not exactly opposite
    /// @param limit the search stops as soon as more than this many are found
    /// @param found receives the pixels, in the order met; cleared first
    /// @return found.size(), which exceeds @p limit only when the search
    ///   stopped early
    std::size_t arc_unobserved_pixels (const Healpix_Map<float>& map,
                                       const vec3& a, const vec3& b,
                                       unsigned limit, std::vector<int>& found);

  }

}

#endif
