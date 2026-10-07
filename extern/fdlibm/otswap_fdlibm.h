/* SPDX-License-Identifier: GPL-2.0-or-later */
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

/*
 * Forced-include header of the vendored FreeBSD msun sources: the build
 * passes it with -include to every file of msun/src, which stay unmodified.
 *
 * It renames every external symbol those files define, so that nothing
 * clashes with, or is resolved against, the system libm; src/detmath.h
 * declares the renamed functions. It also supplies the few FreeBSD
 * definitions the files need and that glibc and macOS lack, and it
 * disables the long double aliases, which would otherwise export
 * unrenamed libm names where long double is double.
 *
 * Not part of the public API.
 */

#ifndef OTSWAP_FDLIBM_H
#define OTSWAP_FDLIBM_H

#include <stdint.h>
#include <sys/types.h>
#include <sys/cdefs.h>

/* The sources assume that double arithmetic is evaluated in double: so it
   is with FLT_EVAL_METHOD 0, and with 16, which GCC 13 and 14 report on
   AArch64 and which widens only types narrower than _Float16. */
#if defined(__FLT_EVAL_METHOD__) && __FLT_EVAL_METHOD__ != 0 && __FLT_EVAL_METHOD__ != 16
#error "the vendored msun sources need FLT_EVAL_METHOD 0 or 16"
#endif

#undef  __weak_reference
#define __weak_reference(sym, alias)

/* FreeBSD's __CONCAT expands its arguments before pasting them, which
   math_private.h relies on (0x1.8p with LDBL_MANT_DIG); the glibc and macOS
   ones paste them as written. */
#undef  __CONCAT
#define OTSWAP_FDLIBM_CONCAT1(x, y) x ## y
#define __CONCAT(x, y)              OTSWAP_FDLIBM_CONCAT1(x, y)

#ifndef __always_inline
#define __always_inline inline __attribute__((__always_inline__))
#endif

#define __double_t double
#define __float_t  float

#define sin                 otswap_fdlibm_sin
#define cos                 otswap_fdlibm_cos
#define exp                 otswap_fdlibm_exp
#define log1p               otswap_fdlibm_log1p
#define pow                 otswap_fdlibm_pow
#define atan                otswap_fdlibm_atan
#define atan2               otswap_fdlibm_atan2
#define asin                otswap_fdlibm_asin
#define __kernel_sin        otswap_fdlibm_kernel_sin
#define __kernel_cos        otswap_fdlibm_kernel_cos
#define __kernel_rem_pio2   otswap_fdlibm_kernel_rem_pio2
#define __ieee754_rem_pio2  otswap_fdlibm_ieee754_rem_pio2

#endif
