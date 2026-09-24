# Vendored Healpix C++ sources

The 21 files under `Healpix_cxx/` and `cxxsupport/` are copied from
**HEALPix 3.82**, `src/cxx/Healpix_cxx/` and `src/cxx/cxxsupport/`, byte for
byte and **unmodified**. They are the subset otswap needs: `T_Healpix_Base`,
`Healpix_Map`, `pointing`, `vec3`, `fix_arr` and `rangeset`, with the
out-of-line code those depend on.

- Source: HEALPix 3.82, <https://healpix.sourceforge.io>
- License: GPL-2.0-or-later, as stated in each file; `COPYING` is the copy
  shipped with Healpix's C++ package.
- Copyright: Max-Planck-Society (see each file's header).

They are third-party code: do not edit them. Update by copying the same files
from a newer HEALPix release. otswap builds them as the static internal
library `otswap_healpix`, with warnings suppressed for these files only; they
need nothing besides the C++ standard library.
