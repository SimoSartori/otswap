/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
 ********************************************************************/

/**
 *  @file src/arc.h
 *
 *  @brief The search for unobserved pixels along a great circle arc.
 *  Internal, and the only internal header that needs Healpix.
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

    /// True when the pixel is observed: its value exceeds kMaskAllowedAbove.
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
