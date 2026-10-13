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
 *  @file src/Maths.cpp
 *
 *  @brief The functions of otswap/Maths.h, from the vendored FreeBSD msun
 *  sources in extern/fdlibm.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include "otswap/Maths.h"

extern "C" {
  double otswap_fdlibm_sin   (double);
  double otswap_fdlibm_cos   (double);
  double otswap_fdlibm_exp   (double);
  double otswap_fdlibm_log1p (double);
  double otswap_fdlibm_cbrt  (double);
  double otswap_fdlibm_asin  (double);
  double otswap_fdlibm_pow   (double, double);
  double otswap_fdlibm_atan2 (double, double);
}


// ============================================================================


double otswap::maths::sin (const double x)
{
  return otswap_fdlibm_sin(x);
}


// ============================================================================


double otswap::maths::cos (const double x)
{
  return otswap_fdlibm_cos(x);
}


// ============================================================================


double otswap::maths::exp (const double x)
{
  return otswap_fdlibm_exp(x);
}


// ============================================================================


double otswap::maths::log1p (const double x)
{
  return otswap_fdlibm_log1p(x);
}


// ============================================================================


double otswap::maths::cbrt (const double x)
{
  return otswap_fdlibm_cbrt(x);
}


// ============================================================================


double otswap::maths::asin (const double x)
{
  return otswap_fdlibm_asin(x);
}


// ============================================================================


double otswap::maths::pow (const double x, const double y)
{
  return otswap_fdlibm_pow(x, y);
}


// ============================================================================


double otswap::maths::atan2 (const double y, const double x)
{
  return otswap_fdlibm_atan2(y, x);
}
