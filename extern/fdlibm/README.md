# Vendored deterministic math: FreeBSD msun (fdlibm lineage)

The 13 files under `msun/src/` are copied from the **FreeBSD source tree**,
`lib/msun/src/`, byte for byte and **unmodified**:

- Repository: <https://github.com/freebsd/freebsd-src>
- Commit: `481447ad9be2e60e0bb3744ab199ce121baadf84` (branch `main`,
  committed 2026-10-05T20:22:44Z)
- Downloaded from `https://raw.githubusercontent.com/freebsd/freebsd-src/<commit>/lib/msun/src/<file>`

| File | Provides | SHA-256 |
|---|---|---|
| `s_sin.c` | `sin` | `acb2837b21fd365ff7cb21d3749321e9bd19f50425c69f79e33064ec9fd12987` |
| `s_cos.c` | `cos` | `b9ebcab1e4d5813404d61f8d2f563e1d3a1793d45d481a064762aa19c4ddd818` |
| `k_sin.c` | `__kernel_sin` | `28388c44b0b3f59c72948452b45c8c60aded8d658be27bf1a9568c09ba63de09` |
| `k_cos.c` | `__kernel_cos` | `e1d3c819b635a744fcf09d082ebffc26a8464474ce4da567cd36464bc8772ddf` |
| `e_rem_pio2.c` | `__ieee754_rem_pio2`, included by `s_sin.c` and `s_cos.c`, not compiled alone | `279da18436f04c773b200a5eaa9fd2facf83b8a0e3873afe754f24f484a0d627` |
| `k_rem_pio2.c` | `__kernel_rem_pio2` | `1ce6c260dd5327b96fdc8dc211a41651fbb71e36fa6e529170ab4da258d1b910` |
| `e_exp.c` | `exp` | `902371ecd772951c2604bd9fcdf121247debb0c59c09d119817c77fc0474c5c6` |
| `s_log1p.c` | `log1p` | `6fb0248e5ade56613bf1ea3fe392a7e255383b953e9156717048a87340fd4da6` |
| `e_pow.c` | `pow` | `b250a94b6816382dceadc8890c3a4990782dd8a1bdc8bbd6192e8f6718106d58` |
| `e_atan2.c` | `atan2` | `cae6f5c5ddf95616c4d6bbfdcc081ac5898d18c3e79ed944a777c8a687540c1d` |
| `s_atan.c` | `atan` | `562febb0fe60383526f7fda1a668f1186b2dcd2e0773b738d57f2ed30bdb6ebb` |
| `e_asin.c` | `asin` | `4340fb7ca58f3b84a8e8c395b2fa0428ff24baa4a313302ea0d6ffc984ba775a` |
| `math_private.h` | word-access macros, `rnint`, kernel prototypes | `6a90b768df36af4f83f4b51e9a94d8d289f1560365371d3a5c7282ef4dd57cad` |

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
`s_log1p.c` and `math_private.h`; `e_rem_pio2.c` adds "Optimized by Bruce
D. Evans.") `e_exp.c` and `e_pow.c` read:

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

- As the static internal library `otswap_fdlibm`, in C, with warnings
  suppressed for these files only, and with `-ffp-contract=off` like every
  other target.
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
- otswap calls them only through `src/detmath.h`.

They are third-party code: do not edit them. Update by copying the same
files from a newer commit, and record it here.
