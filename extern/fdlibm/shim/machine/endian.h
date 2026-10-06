/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it */

/*
 * The vendored msun/src/math_private.h includes <machine/endian.h> for
 * BYTE_ORDER, LITTLE_ENDIAN and BIG_ENDIAN. glibc has no such header and
 * defines them in <endian.h>. Added to the include path on systems other
 * than macOS, whose own <machine/endian.h> is used there.
 */

#ifndef OTSWAP_FDLIBM_MACHINE_ENDIAN_H
#define OTSWAP_FDLIBM_MACHINE_ENDIAN_H

#include <endian.h>

#if !defined(BYTE_ORDER) || !defined(LITTLE_ENDIAN) || !defined(BIG_ENDIAN)
#error "<endian.h> does not define BYTE_ORDER, LITTLE_ENDIAN and BIG_ENDIAN"
#endif

#endif
