# Changelog

## Unreleased

### Added

- `io::readKeywords` (Python `io.read_keywords`) reads the keywords of a
  table, ASCII or FITS, from the first extension or one named by EXTNAME.

### Changed

- Every translation unit is compiled with `-fno-math-errno` as well as
  `-ffp-contract=off` (except with MSVC). The functions it lets the compiler
  inline are exact, so results are unchanged, bit for bit.

## 0.2.0 — 2026-10-23

### Added

- The `otswap` executable: a reconstruction, and optionally its
  redshift-space correction, run from a parameter file, with every key and
  the output files on the new "The executable" page. Built and installed into
  `bin/` by a standalone build of the C++ library (`OTSWAP_BUILD_EXECUTABLE`,
  off when otswap is a subproject); not in the Python package.
  `examples/executable/` has a commented parameter file per geometry.
- `otswap::io::ParameterFile`, in the new header `otswap/ParameterFile.h`: a
  reader of `key = value` parameter files, with typed reads of strings,
  booleans, integers, doubles, characters and lists of each, errors naming
  the file, line and key, the line that sets a key, and the list of the keys
  never read. The syntax is on the "File formats" page. C++ only.
- The cosmology constructor of `DistanceTable` defaults `w0` to -1 and `wa`
  to 0, flat LCDM; Python's `DistanceTable.flat` already did.
- The header `otswap/Maths.h`: `otswap::maths::sin`, `cos`, `exp`, `log1p`,
  `cbrt`, `asin`, `pow` and `atan2`, which return the same bits on every
  supported platform and compiler, within one unit in the last place.
  `cbrt` is vendored from the same FreeBSD commit as the others. C++ only.
- `measureMpsProfile` (Python: `measure_mps_profile`): the mps(z) profile
  `reconstructLightcone` would use, measured from the tracers alone through
  the same redshift cut, mask and checks, bit for bit. In Python
  `angle_unit` is required.
- Optional box limits in `reconstructBox`: a `BoxLimits` (Python: `box=`, a
  `BoxLimits` or a `(lower, upper)` pair), closed on every face. Tracers
  outside are flagged in `Result::outsideBox` (`outside_box`) and left out,
  given randoms outside are dropped, and drawn randoms and the mean particle
  separation come from the limits. `Result::box` (`box`) records the box
  used, the limits or the tracers' bounding box. `SelectionCounts` gains
  `boxLimits`, `tracersOutsideBox` and `randomsOutsideBox`, and its report a
  line for each. `realSpaceBox` leaves the tracers outside with status
  `LeftOut`. The mean displacement file of a box with limits adds the
  `OutsideBox` group (`outside_box`) to its default, 11 columns. Without
  limits, results are unchanged.
- Object IDs: `Result::ids` and `RealSpaceCatalog::ids` (Python:
  `Result.ids`, settable and checked on assignment, and
  `RealSpaceCatalog.ids`), 64-bit integers, one per object or empty for
  0..N-1. `realSpaceBox` and `realSpaceLightcone` copy them into the
  catalogue, and the writers write them as the `id` column (see Changed).
  They enter no computation. `Result.from_arrays` takes `ids=`.
- `DistanceTable::redshift`, `distance` and `growthRate` (Python: the
  properties `redshift`, `distance` and `growth_rate`) return the table's
  arrays.
- The complete displacement file is self-contained. It adds, for a
  lightcone, the tracers' `tracRA`, `tracDec`, `tracRed` (right ascension as
  given, not folded into [0, 360)), the mps(z) profile and the distance
  table; and `outsideBox` in both geometries. In FITS the profile and the
  table are the HDUs `MPSPROF` and `DISTANCE` after the main table
  `DISPLACEMENTS`; in ASCII they are twin files, `run.dat` giving
  `run.mps.dat` and `run.distances.dat`, each with its own `PRODUCT`,
  declared by the keywords `MPSPART` and `DISTPART`. `lightconeTables`
  (`lightcone_tables`) false leaves them out.
- Both displacement files record the whole configuration (`CELLSIZE`,
  `REJCROSS`, `CFMAXPIX`, `VERBOSE` added to `NREC`, `SEED`, `CONVERG`), the
  box (`BOXXLO` to `BOXZHI`, `BOXLIM`) and the selection counts (`NTRAC` to
  `NDISPCR`, each when not 0). In FITS their main table is named
  `DISPLACEMENTS` (`EXTNAME`). Both writers refuse a result whose
  `nRealizations` and `config.nRealizations` differ, or whose box limits and
  `box` differ.
- `io::readDisplacementComplete` (Python: `io.read_displacement_complete`,
  with a required `angle_unit`): rebuilds the `Result` from a complete
  displacement file, of this version or of otswap 0.1. Every double reads
  back identical, NaN as NaN (with its sign and payload in FITS only), but
  for `tracersSky`, equal within the rounding of the conversion to degrees
  and back. In a file of otswap 0.1 the `index` column is read as the IDs.

### Changed

- Requires Python 3.10 or later, was 3.9: nanobind 3 does not support
  3.9. No wheel is built for CPython 3.9. Building the Python package
  requires nanobind 3.1.0 exactly, was 2.0 or later, the version the wheels
  are built with.
- `Config::nRealizations` (`n_realizations`) defaults to 30, was 1.
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
  removed. `measure_mps_profile` keeps its own `verbosity`,
  `reject_mask_crossings` its `max_unobserved_pixels_crossed` and
  `neighbour_average` its `weight_by_realizations`.
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
- An error of the neighbour search reaches the caller as an
  `otswap::Error` (`otswap.Error`) whose message starts with `meshsearch: `,
  the original nested in it (`std::rethrow_if_nested`). C++ code that caught
  `meshsearch::Error`, and Python code that caught the `RuntimeError` of
  `neighbour_average` and the real-space functions, must catch it instead.
- `rejectMaskCrossings` (`reject_mask_crossings`) refuses a box result, one
  built by hand included: a `Result` defaults to the box geometry, and
  `Result.from_arrays` without `tracers_sky` gives one.
- The `index` column is now `id`: the object's ID, or its row when no IDs
  are set, so the values are unchanged when none are. The groups
  `DisplacementGroup::Index` and `CatalogGroup::Index` are renamed `Id`
  (Python: `"index"` becomes `"id"`), with no alias: code that names them no
  longer compiles, and Python raises `otswap.Error` for the unknown group.
  The default columns start with `id` only when IDs are set, so the default
  files of 0.1 are unchanged.
- Descriptions in every file: `status` of the real-space catalogue now
  reads "0 corrected, 1 moved by neighbours, 2 no valid neighbour, 3 left
  out", short enough for a FITS card, and `NOBJECTS` "number of input
  objects". Values and columns are unchanged.
- `io::writeDisplacements` and `io::writeDisplacementField` are renamed
  `writeDisplacementMean` and `writeDisplacementComplete` (Python:
  `write_displacement_mean`, `write_displacement_complete`). The complete
  file's `PRODUCT` is `complete displacements`, was `displacement field`.

### Removed

- `Result::matchedRandom` (Python: `Result.matched_random`), the random
  matched to each tracer in each realization; the tracer plus its
  displacement gives it within rounding. A result with the default 30
  realizations takes about 45% less memory.
  `Result.from_arrays` no longer takes it: its positional arguments are now
  `displacement`, `valid`, `tracers`. `rejectMaskCrossings` takes the arcs
  from the tracers (`Result::tracers`, which it now requires) to the tracer
  plus its displacement.
- The columns lagrX, lagrY, lagrZ of the complete displacement file, which
  held the matched randoms. `readDisplacementComplete` ignores them in a
  file of otswap 0.1.

### Fixed

- Python: `io.write` reports an out-of-range precision with the value
  given; a value above 1000 was reported as 1000.
- The documentation of `reconstructLightcone` stated the random supply rule
  before the selection: the randoms kept by the redshift cut and the mask
  must number at least `nRealizations` times the tracers kept.

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
