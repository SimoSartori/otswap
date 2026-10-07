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
 *  @file src/pixel.h
 *
 *  @brief The Healpix pixelisation: geometry, pixel lookup and neighbours,
 *  on a 64-bit pixel index. Internal.
 *
 *  Derived from HEALPix 3.82 (Healpix_cxx/healpix_base.cc,
 *  healpix_tables.cc, cxxsupport/math_utils.h and vec3.h), Copyright (C)
 *  2003-2016 Max-Planck-Society, author Martin Reinecke, distributed under
 *  the GNU General Public License, version 2 or later. The floating-point
 *  operations are those of Healpix, in the same order, so every pixel is
 *  the one Healpix computes from the same inputs.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_PIXEL_H
#define OTSWAP_PIXEL_H

#include <array>
#include <cmath>
#include <cstdint>

namespace otswap {

  namespace internal {

    /// The largest NSIDE a mask may have, 2^29: the limit of Healpix's
    /// 64-bit pixel index.
    constexpr std::int64_t kMaxNside = std::int64_t(1) << 29;

    /// A Cartesian vector, with the operations of Healpix's vec3 that the
    /// arc search uses, computed as Healpix computes them.
    struct Vec3 {
      double x = 0., y = 0., z = 0.;

      Vec3 () = default;
      Vec3 (const double xc, const double yc, const double zc) : x(xc), y(yc), z(zc) {}

      double length () const { return std::sqrt(x*x + y*y + z*z); }

      double squared_length () const { return x*x + y*y + z*z; }

      /// The vector scaled to length 1, by the reciprocal of its length.
      Vec3 normalized () const
      {
        const double l = 1./std::sqrt(x*x + y*y + z*z);
        return Vec3(x*l, y*l, z*l);
      }

      Vec3 operator+ (const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
      Vec3 operator- (const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    };

    /// The Healpix pixelisation of the sphere at one NSIDE, in the RING or
    /// the NESTED scheme, with 64-bit pixel indices.
    class HealpixBase {

    public:

      /// Sets the geometry. @p nside must lie in [1, kMaxNside], and be a
      /// power of 2 when @p nested; check_nside checks both.
      void set (std::int64_t nside, bool nested);

      std::int64_t nside () const { return m_nside; }

      /// The number of pixels, 12 NSIDE^2.
      std::int64_t npix () const { return m_npix; }

      bool nested () const { return m_nested; }

      /// Pixel at z = cos(theta) and longitude @p phi, radians; @p sth is
      /// sin(theta), used near the poles when @p haveSth.
      /// @exception Error if the in-ring index comes out of range, which
      /// does not happen for finite inputs with |z| <= 1.
      std::int64_t loc2pix (double z, double phi, double sth, bool haveSth) const;

      /// The eight neighbours of @p pix, in the order SW, W, NW, N, NE, E,
      /// SE, S; -1 where a pixel at a face corner has only seven.
      /// @exception Error if a neighbour's ring position comes out of range,
      /// which does not happen for a pixel of this geometry.
      void neighbors (std::int64_t pix, std::array<std::int64_t, 8>& result) const;

    private:

      int m_order = -1;            ///< log2(NSIDE), or -1 when NSIDE is not a power of 2
      std::int64_t m_nside = 0;
      std::int64_t m_npface = 0;   ///< pixels per face, NSIDE^2
      std::int64_t m_ncap = 0;     ///< pixels in each polar cap
      std::int64_t m_npix = 0;
      bool m_nested = false;

      std::int64_t xyf2nest (int ix, int iy, int face) const;
      void nest2xyf (std::int64_t pix, int& ix, int& iy, int& face) const;
      std::int64_t xyf2ring (int ix, int iy, int face) const;
      void ring2xyf (std::int64_t pix, int& ix, int& iy, int& face) const;
      void ring_info_small (std::int64_t ring, std::int64_t& startpix, std::int64_t& ringpix,
                            bool& shifted) const;
    };

    /// Pixel holding the direction @p v, which need not be normalised:
    /// Healpix's vec2pix, with the atan2 of detmath.h in place of the
    /// system one, so that the pixel is the same on every platform. v must
    /// be finite and non-zero.
    std::int64_t vec2pix (const HealpixBase& base, const Vec3& v);

    /// Pixel holding colatitude @p theta and longitude @p phi, radians:
    /// Healpix's ang2pix, with the cos and sin of detmath.h, and without
    /// its assertion: theta must lie in [0, pi] and phi be finite, which
    /// the caller checks.
    std::int64_t ang2pix (const HealpixBase& base, double theta, double phi);

  }

}

#endif
