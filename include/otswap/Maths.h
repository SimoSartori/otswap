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
 *  @file include/otswap/Maths.h
 *
 *  @brief Elementary functions that return the same bits on every
 *  platform.
 *
 *  Each returns the same bits on every supported platform and compiler,
 *  with any number of threads, whatever the system's maths library, and is
 *  within one unit in the last place of the exact value. Use them in place
 *  of the standard library's where a result must be reproducible bit for
 *  bit. sqrt, fabs, floor and fmod need no counterpart here: IEEE 754
 *  already makes them exact or correctly rounded everywhere.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_MATHS_H
#define OTSWAP_MATHS_H

namespace otswap {

  namespace maths {

    /// @brief The sine of x, in radians.
    double sin (double x);

    /// @brief The cosine of x, in radians.
    double cos (double x);

    /// @brief e raised to the power x.
    double exp (double x);

    /// @brief The natural logarithm of 1 + x, accurate also for x near 0.
    double log1p (double x);

    /// @brief The cube root of x, negative for negative x.
    double cbrt (double x);

    /// @brief The arc sine of x, in radians, in [-pi/2, pi/2].
    double asin (double x);

    /// @brief x raised to the power y.
    double pow (double x, double y);

    /// @brief The angle of the point (x, y) from the x axis, in radians, in
    /// [-pi, pi].
    double atan2 (double y, double x);

  }

}

#endif
