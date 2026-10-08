# Changelog

## 0.1.0 — 2026-10-08

First release.

### Compared with CosmoBolognaLib's OTreconstruction

- The same algorithm and inputs: seeding on a mesh, then local search that
  tries every permutation of the partners of four neighbouring tracers.
- With a fixed seed, the result does not depend on the number of threads, and
  is the same bit for bit on every supported platform and compiler: the
  assignment of the randoms, their draw, the seeding and the swap loop each
  draw from their own random stream.
- Every tracer is examined in every sweep.
- The seeding radius (8 mean particle separations) and the number of swap
  neighbours (113) are named constants, not parameters.
- Mask crossings: the adaptive search along the arc, with a limit on the
  distinct unobserved pixels crossed (0 rejects at the first), and validity
  kept per realization.
- Cosmology enters as a table of redshift, comoving distance and growth rate,
  computed from a flat w0-wa cosmology or given by the caller; no GSL.
- The redshift-space correction is a library function for each geometry, with
  its four steps also available on their own.
- Table writers for every product, ASCII and FITS, with column descriptions,
  units and header keywords.
- A Python package on the same code.
