# Changelog

## Unreleased

### Added

- `otswap::io::ParameterFile`, in the new header `otswap/ParameterFile.h`: a
  reader of `key = value` parameter files, with typed reads of strings,
  booleans, integers, doubles, characters and lists of each, errors naming
  the file, line and key, and the list of the keys never read. The syntax is
  on the "File formats" page. C++ only.

### Changed

- Requires meshsearch 0.2.1 or a later 0.2.x; 0.2.1 is fetched when none is
  installed.
- Tracers at exactly the same distance from a tracer, or randoms at exactly
  the same distance from a first-guess starter, are now ordered by ascending
  index. Their order used to come from the standard library, so a result
  could differ between compilers when distances were tied, as on a lattice.
  Results for inputs with exactly tied distances change on every compiler;
  results without ties are unchanged.
- The error messages for a wrong value of a configuration field, of the
  number of redshift bins, of a cosmological parameter, redshift range or
  sample count of a distance table, or of the smoothing width sigma now
  describe the quantity instead of naming the C++ field or argument, so they
  read the same from C++ and Python.
- Python: the parameters of a reconstruction are an `otswap.Config`, passed
  as `config=` to `reconstruct_box` and `reconstruct_lightcone`, and the
  options of the correction an `otswap.CorrectionConfig`, passed as
  `correction=` to `real_space_lightcone` and `real_space_box`; both are the
  C++ structs, with every field a keyword of the constructor and a
  read/write attribute. The keywords they replace (`n_realizations`,
  `convergence`, `seed`, `cell_size`, `reject_crossings`,
  `max_unobserved_pixels_crossed`, `verbosity`, `weight_by_realizations`) are
  removed. `measure_mps_profile` keeps `verbosity`, `reject_mask_crossings`
  and `neighbour_average` their own arguments.
- Python: `otswap.RedshiftCut` is a class; `redshift_cut=` takes one, a
  `(min, max)` pair or None, and `SelectionCounts.redshift_cut` returns one.
  `BoxLimits` and `RedshiftCut` compare equal by value.
- Python: `Result.config` and `Result.distances` return copies of the
  result's configuration and distance table; `Result.seed` is removed, its
  value is `Result.config.seed`. `Result.from_arrays` takes `config=`, in
  place of `seed=`.
- Python: `Result` and `RealSpaceCatalog` are the C++ structs too, so
  another nanobind module can take and return them, as it can `Config`,
  `CorrectionConfig`, `RedshiftCut`, `BoxLimits`, `Mask`, `DistanceTable`,
  `MpsProfile` and `BiasTable` (conditions in the "Names side by side" page).
  One made by another module is in radians, its `angle_unit` `"rad"`; its
  array properties raise `otswap.Error` when the C++ arrays do not have the
  sizes the object declares.

## 0.1.0 — 2026-10-08

First release.

### Compared with CosmoBolognaLib's OTreconstruction

- The same algorithm and inputs: a local first guess, then local search that
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
