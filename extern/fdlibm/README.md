# Vendored deterministic math: FreeBSD msun (fdlibm lineage)

The 14 files under `msun/src/` are copied from the **FreeBSD source tree**,
`lib/msun/src/`, byte for byte and **unmodified**:

- Repository: <https://github.com/freebsd/freebsd-src>
- Commit: `481447ad9be2e60e0bb3744ab199ce121baadf84` (branch `main`,
  committed 2026-10-05T20:22:44Z)

| File | Provides |
|---|---|
| `s_sin.c` | `sin` |
| `s_cos.c` | `cos` |
| `k_sin.c` | `__kernel_sin` |
| `k_cos.c` | `__kernel_cos` |
| `e_rem_pio2.c` | `__ieee754_rem_pio2`, included by `s_sin.c` and `s_cos.c`, not compiled alone |
| `k_rem_pio2.c` | `__kernel_rem_pio2` |
| `e_exp.c` | `exp` |
| `s_log1p.c` | `log1p` |
| `e_pow.c` | `pow` |
| `e_atan2.c` | `atan2` |
| `s_atan.c` | `atan` |
| `e_asin.c` | `asin` |
| `s_cbrt.c` | `cbrt` |
| `math_private.h` | word-access macros, `rnint`, kernel prototypes |

FreeBSD's own `math.h` is not vendored: it depends on FreeBSD's system
headers. The files' `#include "math.h"` resolves to the system `<math.h>`,
used only for declarations and `double_t`.

From the system libm the files still call `sqrt`, `fabs`, `floor` and
`scalbn`. IEEE 754 requires all four to be exact or correctly rounded, so
they give the same result on every conforming platform.

## Licence

Every vendored file carries the Sun Microsystems notice below, and no
other licence or copyright line. The 1993 files read:

```
Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.

Developed at SunSoft, a Sun Microsystems, Inc. business.
Permission to use, copy, modify, and distribute this
software is freely granted, provided that this notice
is preserved.
```

(`SunPro` in place of `SunSoft` in `s_sin.c`, `s_cos.c`, `s_atan.c`,
`s_log1p.c`, `s_cbrt.c` and `math_private.h`; `e_rem_pio2.c` and `s_cbrt.c`
add "Optimized by Bruce D. Evans.") `e_exp.c` and `e_pow.c` read:

```
Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.

Permission to use, copy, modify, and distribute this
software is freely granted, provided that this notice
is preserved.
```

This is a permissive licence whose only condition is that the notice be
preserved, which it is, in each file and here. It is compatible with
otswap's GPL-2.0-or-later.

## How otswap builds them

- As the internal object library `otswap_fdlibm`, in C, with warnings
  suppressed for these files only, and with `-ffp-contract=off` like every
  other target. Its objects go into `libotswap.a`; nothing of it is
  installed or exported on its own.
- `otswap_fdlibm.h`, otswap's own file, is force-included in every one of
  them. It renames each external symbol to `otswap_fdlibm_<name>`, so
  nothing clashes with the system libm, and supplies the FreeBSD
  definitions glibc and macOS lack or define differently (`__double_t`,
  `__float_t`, `__always_inline`, a `__CONCAT` that expands its arguments),
  and an empty `__weak_reference`, which would
  otherwise export unrenamed long double aliases where long double is
  double.
- `shim/machine/endian.h`, also otswap's, stands in for FreeBSD's
  `<machine/endian.h>` where the system has none (glibc).
- otswap calls them only through `src/detmath.h`, for its own results,
  and `src/Maths.cpp`, which defines the public `otswap/Maths.h`.

They are third-party code: do not edit them. Update by copying the same
files from a newer commit, and record it here.
