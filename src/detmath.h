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
 *  @file src/detmath.h
 *
 *  @brief The elementary functions on otswap's result path, computed by
 *  the same code on every platform. Nothing here is part of the public
 *  API.
 *
 *  They come from FreeBSD's msun, vendored unmodified in extern/fdlibm and
 *  compiled with contraction off, so they return the same bits wherever
 *  otswap is built, whatever the system libm. Their error is below one
 *  unit in the last place. Every transcendental function otswap evaluates
 *  on a value that reaches a result goes through here; sqrt, fabs, fmod and
 *  floor need not, since IEEE 754 makes them exact or correctly rounded.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_DETMATH_H
#define OTSWAP_DETMATH_H

extern "C" {
  double otswap_fdlibm_sin   (double);
  double otswap_fdlibm_cos   (double);
  double otswap_fdlibm_exp   (double);
  double otswap_fdlibm_log1p (double);
  double otswap_fdlibm_pow   (double, double);
  double otswap_fdlibm_atan2 (double, double);
  double otswap_fdlibm_asin  (double);
}

namespace otswap {

  namespace internal {

    inline double det_sin   (const double x)                 { return otswap_fdlibm_sin(x); }
    inline double det_cos   (const double x)                 { return otswap_fdlibm_cos(x); }
    inline double det_exp   (const double x)                 { return otswap_fdlibm_exp(x); }
    inline double det_log1p (const double x)                 { return otswap_fdlibm_log1p(x); }
    inline double det_pow   (const double x, const double y) { return otswap_fdlibm_pow(x, y); }
    inline double det_atan2 (const double y, const double x) { return otswap_fdlibm_atan2(y, x); }
    inline double det_asin  (const double x)                 { return otswap_fdlibm_asin(x); }

  }

}

#endif
