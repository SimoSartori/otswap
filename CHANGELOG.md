# Changelog

## 0.4.0 — 2026-11-26

Changes marked **Breaking** require changes to code or builds that worked
with 0.3.

### Added

- Periodic boxes. `reconstructBox` takes `periodic` after the limits
  (Python `reconstruct_box(..., periodic=True)`), which it then requires;
  each side must be an integer multiple of the shortest, within a relative
  1e-12. Distances are minimum-image; tracers and given randoms outside the
  limits are wrapped into them by whole periods instead of being left out;
  `Result::tracers` holds the wrapped tracers, each displacement is the
  minimum-image vector, in [−L/2, L/2) per axis within rounding, and
  `lagrangian` is wrapped into the box. A computed mps is (V/N)^(1/3), V the
  volume of the limits and N the number of tracers given, all of which are
  kept. `Result::periodic` (`Result.periodic`)
  records the flag; `SelectionCounts` gains `periodic`, `tracersWrapped` and
  `randomsWrapped` (`tracers_wrapped`, `randoms_wrapped`), reported at
  Normal as "wrapped <a> of <N> tracers into the periodic box ...".
  `Result.from_arrays` takes `periodic=` and `box=`; `box=` alone records
  the limits of a box result that is not periodic.
- `Result::cellSide` (`Result.cell_side`): the cell side of the neighbour
  grids actually used, in coordinate units. In a periodic box it can be
  smaller than `cellSize` times the mps, since the cells divide the shortest
  side a whole number of times; the error for a grid too large to address
  covers that reduced side. Not written to files; NaN, None in Python, for a
  result not made by a reconstruction.
- Examples `periodic_box`, a reconstruction and correction of halos in
  redshift space in a periodic box, and `lightcone_cartesian`, which also
  finds each tracer's distance and direction with `pointToSky` and
  measures the mps profile alone, in C++ and Python, and
  `examples/executable/periodic_box.ini`; the data
  `periodic_halos.dat`, real- and redshift-space positions.
- Periodic correction of redshift-space distortions: `realSpaceBox`
  (`real_space_box`) of a periodic result finds the neighbours across the
  faces and wraps the corrected coordinate into the box. `neighbourAverage`
  gains an overload with the box and the flag (Python `neighbour_average`
  with `box=` and `periodic=True`). `shiftAlongAxis` still does not wrap.
- Lightcones in comoving Cartesian coordinates, the observer at the origin,
  without a cosmology: `reconstructLightcone(tracers, randoms, skyAreaDeg2
  or mask, nBins, config, distances = nullptr, cut = {})` and the two
  matching `measureMpsProfile` overloads (Python
  `reconstruct_lightcone_cartesian`, `measure_mps_profile_cartesian`). The
  mps is measured in bins of comoving distance from the origin, a profile in
  distance (`MpsProfile::distance`, `distanceMin`, `distanceMax`; Python
  `distance`, `distance_min`, `distance_max`); the mask
  is applied on the angles of `pointToSky`. The optional distance table
  plays no part in the reconstruction: it fills the redshifts of
  `tracersSky` and `lagrangianSky`, NaN without it, and is kept in the
  result. `Result::lagrangianDistance` (`lagrangian_distance`) holds the
  distance of each mean Lagrangian position.
- `DistanceCut` (Python `DistanceCut`, keyword `distance_cut`, or a pair
  `(min, max)`), the cut in comoving distance of the Cartesian overloads.
  Tracers outside it keep their row, flagged in `Result::outsideDistanceCut`
  (`outside_distance_cut`); `SelectionCounts` gains `distanceCut`,
  `tracersOutsideDistanceCut` and `randomsOutsideDistanceCut`
  (`distance_cut`, `tracers_outside_distance_cut`,
  `randoms_outside_distance_cut`).
  `Result.from_arrays` takes `outside_distance_cut=`.
- `realSpaceLightcone` (`real_space_lightcone`) corrects the result of a
  Cartesian overload: a tracer's redshift is read from `tracersSky` when it
  is not NaN, and otherwise derived from its distance with the distance
  table given. Results of the sky overloads are unchanged.
- `pointToSky` (Python `point_to_sky`): the distance, right ascension and
  declination of a Cartesian point, without a distance table; the geometric
  part of `toSky`, bit for bit.
- Files. The displacement files write `PERIODIC = 1` for a periodic box,
  with the box keywords at full precision in ASCII too, and the counts
  `NTRWRAP` and `NRNWRAP`; for a lightcone given in Cartesian coordinates,
  the column `outsideDistanceCut` (in the `Selection` group of the mean
  file and last in the complete file), the keywords `DCUTMIN`, `DCUTMAX`,
  `NTRDCUT` and `NRNDCUT`, and a profile in distance: the column `distance`
  in place of `redshift`, the keywords `DISTMIN` and `DISTMAX` in place of
  `ZMIN` and `ZMAX`, also in `writeMpsProfile`. Without a distance table no
  `DISTANCE` part is written. The readers read all of these back, and the
  files of otswap 0.3 as before.
- The executable takes four new keys, 45 in all: `periodic` for a box,
  `lightcone_input` (`sky` or `cartesian`) for a lightcone, and `d_min`,
  `d_max`, the distance cut of a Cartesian lightcone. With
  `lightcone_input = cartesian` the cosmology is optional unless
  `perform_recons` is set. `examples/executable/lightcone_cartesian.ini`
  runs the lightcone example on its catalogues converted to Cartesian
  coordinates, which `examples/data/` now holds.

### Changed

- **Breaking (CMake):** requires meshsearch 0.3.0 or a later 0.3.x, was 0.2.1 or a
  later 0.2.x; v0.3.0 is fetched when none is installed.
- **Breaking (C++, Python modules sharing the types):** the layouts of `Result`, `SelectionCounts` and `MpsProfile`
  gain the fields above. A nanobind module or a program built against the
  headers of 0.3 must be rebuilt.
- **Breaking (CMake):** the CMake package answers requests for 0.4
  (`find_package(otswap 0.4 CONFIG)`); a request for 0.3,
  `find_package(otswap 0.3 CONFIG)`, does not find 0.4.0.
- **Breaking (Python):** `MpsProfile.at` takes `value`, a redshift or a
  distance as the profile is: code that called `at(z=...)` must call
  `at(value)` or `at(value=...)`; positional calls are unaffected. In C++
  only the parameter name of `MpsProfile::at` changed. `at` raises when a
  profile holds both redshift and distance nodes.
- `Config::cellSize` (`cell_size`) defaults to 3, was 4, and the
  executable's `mesh_cellsize` with it. The cell size changes the speed
  only: results are unchanged. The files of a run with the default now
  record `CELLSIZE = 3`.
- `CorrectionStatus::LeftOut` (`LEFT_OUT`) also covers the tracers left out
  by the distance cut.
- The displacement writers refuse a periodic result whose flags are not
  consistent (not a box, no limits, a tracer flagged outside the box, or
  `selection.periodic` different from `periodic`);
  `writeMpsProfile` refuses a profile holding both redshift and distance
  nodes.
- `writeDisplacementMean` in ASCII rounds the box keywords `BOXXLO` to
  `BOXZHI` of any box that is not periodic, given by limits or taken from
  the tracers, outward to the 9 significant
  digits of the coordinates, so the tracers as written lie within the
  limits as written. FITS files, periodic boxes and the complete file keep
  the limits exactly. The values of these keywords in ASCII mean files may
  differ from those of 0.3 in the last digit.
- `uniformInt` is defined inline in `otswap/Random.h`; its values are
  unchanged.
- The error for a redshift outside the distance table names the tracer or
  random and its row in the input arrays; from `measureMpsProfile`, which
  takes no randoms, it says the redshifts of the tracers kept must lie in
  the table. The correction of a Cartesian
  lightcone result raises, naming the tracer, when its distance lies
  outside the table, and suggests a table of wider range.
- The examples `box`, `lightcone`, `full_recons_box` and
  `full_recons_lightcone` run every feature they present, with
  `cell_size = 3`: the box example gives the box limits, so its output
  differs from that of 0.3, keeps every realization in a file with IDs and
  reads it back; the lightcone example measures the mps profile alone and
  reads the mean file back. Their parameter files in `examples/executable/`
  set `periodic` or `lightcone_input`.
- Messages say "catalogue", where some said "catalog".

### Removed

- **Breaking (C++, Python):** the overloads of `reconstructLightcone` that
  took both Cartesian and sky arrays (C++), and the `tracers=` and `randoms=`
  arguments of Python `reconstruct_lightcone`, which called them. Arrays made with
  `toCartesian` and the same distance table gave the result of the sky
  overload, which takes their place; for Cartesian coordinates without sky
  ones, use the new Cartesian overloads.

### Build and documentation

- The package metadata lists classifiers: Python 3.10 to 3.14, Linux and
  macOS, astronomy and physics. The supported Python versions are
  unchanged.
- Continuous integration builds and tests the Python package in a clean
  environment with Python 3.10, 3.12 and 3.14, on Linux and macOS.
- `tools/wheels/build_deps.sh` downloads zlib from zlib.net, then its
  fossils directory, then GitHub, with the same SHA-256.
- The guide gives a numbered path of steps for each geometry (box, periodic
  box, lightcone on the sky, lightcone in Cartesian coordinates), each
  linking its example, and sections on the distance table, the inputs of the
  redshift-space correction, and files. It states that a `Result` must not
  be read from one thread while another modifies it in place.
- A new page, "Using otswap from another package", gathers the use of the
  library from CMake, the sharing of the Python types with another nanobind
  module and the loading of the same OpenMP runtime, formerly on the pages
  of names and of the Python reference.
- The Python reference lists the members of each class in the order of the
  stubs, and `BoxLimits` under "Configuration".
- The method page writes the redshift-space correction with the mean
  displacement D = −Ψ_OT and the line of sight l̂, as `otswap/RSD.h` and the
  module docstring do, and describes the first guess and the box as
  implemented. It cites the equations of Sartori et al. (2026) by their
  numbers in arXiv:2601.15378v1, and describes the neighbour average in a
  periodic box, with weights and with σ = 0.
- Each path of the guide ends with the redshift-space correction, linking
  its example, and sets the cut before the choice of the bins.

## 0.3.0 — 2026-11-05

### Added

- The header `otswap/Random.h`: `otswap::uniformInt`, `uniformReal` and
  `shuffle`, random draws over `std::mt19937` that give the same values on
  every supported platform and compiler. Their algorithms are fixed: later
  versions do not change them. `uniformReal` refuses a bound that is not
  finite, a lower bound not below the upper one, or an infinite width, and
  `shuffle` more than 2^32 − 1 elements, each with `otswap::Error`.
- `io::readDisplacementMean` (Python `io.read_displacement_mean`) reads a mean
  displacement file, of this version or of otswap 0.2, into a `Result` of one
  realization holding the mean, with the counts of valid realizations and
  the number of realizations averaged.
- `Result::nRealizationsMean` (Python `Result.n_realizations_mean`), the
  number of realizations averaged into a mean read from file; 1 otherwise.
  The displacement writers record a value above 1 as the new keyword
  NRECMEAN, with NREC = 1, and refuse 0, or a value above 1 with several
  realizations. `validRealizations` is then the count of the realizations
  averaged that were valid, always at most nRealizations x
  nRealizationsMean; `rejectMaskCrossings` on such a result sets it to 0 for
  the objects it rejects and keeps it for the others. `readDisplacementComplete`
  reads NRECMEAN back, 1 when the file has none.
- `io::columnNames` (Python `io.column_names`): the names of the columns of
  a table, those `io::read` finds.
- The mean displacement file of a lightcone holds the mps(z) profile and the
  distance table, as the complete file does: FITS extensions MPSPROF and
  DISTANCE, or ASCII twin files `.mps` and `.distances`, each with its own
  PRODUCT, and the keywords MPSPART and DISTPART. `writeDisplacementMean`
  (Python `write_displacement_mean`) takes `lightconeTables`
  (`lightcone_tables`), true by default, to leave them out.
- `otswap.openmp_root()`: the directory of the OpenMP runtime a macOS wheel
  ships, `otswap/openmp/` (`lib/libomp.dylib` with install name
  `@rpath/libomp.dylib`, `include/omp.h`, `include/ompx.h`, LLVM's
  `LICENSE.TXT`), or None. A package with OpenMP code of its own links that
  file and loads the same runtime as otswap; "Linking the same OpenMP
  runtime" in the Python reference says how.

### Changed

- `rejectMaskCrossings` (Python `reject_mask_crossings`), and the filter of
  `reconstructLightcone` with a mask, judge a displacement by the longest run
  of unobserved pixels its arc crosses one after another, no longer by the
  number of distinct unobserved pixels: `maxUnobservedPixelsCrossed`
  (`max_unobserved_pixels_crossed`, keywords CFMAXPIX and MAXPIXCR) is the
  longest run tolerated. An arc crossing several small holes, each within the
  threshold, is now kept. Results change for thresholds above 0; with 0, the
  default, they are unchanged. The selection report says "unobserved pixels in
  a row", and the comments of CFMAXPIX and MAXPIXCR say so too, now short
  enough to be written whole in FITS.
- `io::read` (Python `io.read`) reads ASCII strictly: every line that is not
  empty or a comment is a data row, and each requested field, and the first
  field of each row, must be a number as a whole (decimal, hexadecimal,
  `nan`, `inf`, `infinity`). Files that read before may now be refused where
  a field was skipped or cut silently: `NA` (the columns after it shifted),
  glued tokens such as `x12`, `1,5` with a blank delimiter (read as two
  fields), `4.5abc`, `1e`, and a header line or any line of text that does
  not start with the comment character, now an error where it was skipped.
  The message names the file, the line and the field. A comment may follow
  blanks, tabs or the delimiter: such a line, read as data before when it
  held a number, is now skipped. A line starting with `###` is the names
  line, skipped as data whatever the comment character, so columns are read
  by name with any comment character. A carriage return ending a line is
  ignored.
- Messages count as the file shows it, and say so: an ASCII line or FITS row
  from 1 ("line 7 (from 1)"), a field from 1, a row `io::write` is given in
  memory from 0. `readBiasTable` names an ASCII entry by its file line, no
  longer its data row, and `readDisplacementComplete` names a row by its
  file line or FITS row from 1, no longer from 0. The checks of coordinate
  arrays and of a distance table built from arrays say "(from 0)". The
  executable names a repeated ID by its file line or FITS row, from 1, no
  longer its data row.
- The other messages that number something say from what too: an object,
  an entry or a redshift bin of an array, from 0 ("object 3 (from 0)"), in
  the coordinate conversions, the correction of redshift-space distortions,
  the checks of a `Result` and the Python arguments; a line or an element of
  a parameter file, from 1.
- `io::write` and `io::writeParts` (Python `io.write`, `io.write_parts`)
  refuse a column description longer than 68 characters, an apostrophe
  counting twice, and a keyword comment longer than the room its FITS card
  leaves: 47 characters after a value within 20 columns, fewer after a
  longer one. FITS cut both silently before; ASCII files are held to the
  same limits. The comment of MPSREPR and the descriptions of lagrRA,
  lagrDec and lagrRed, of displX, displY and displZ in the complete file,
  and of nNeighbours, which FITS cut, are now shorter and written whole.
- Each line of a report is written as soon as its step finishes, the log
  flushed after it, no longer all of them when the call ends; the text is
  unchanged. A call that throws leaves the lines of the steps it finished.
  Python prints each line to `sys.stdout` while the call runs, and Ctrl+C
  stops a reconstruction, `measure_mps_profile` or a correction at its next
  report line with `KeyboardInterrupt`.
- `rejectMaskCrossings` (Python `reject_mask_crossings`) requires finite
  values only in the valid displacements, which it reads, so it accepts a
  NaN in a cleared entry, as a mean read from file holds for a tracer
  without valid realizations.
- `reconstructBox` and `reconstructLightcone` refuse more than 2^32 − 1
  tracers.
- The executable writes the `MPS_` file of `printMps` only when neither
  displacement file is written, those files holding the profile; in ASCII it
  reports the twin files of the mean file too.
- The executable takes `key =`, with no value, as leaving the key out, for
  the keys whose default is none: `id_column`, `random_file` (box), `z_min`,
  `z_max`, `box_min`, `box_max` and `mps`. Every other key still needs a
  value. The parameter files in `examples/executable/` set every key that
  applies to their run, defaults included, and the keys without a default as
  `key =`, instead of leaving them commented out.
- The macOS wheels hold libomp in `otswap/openmp/lib/`, where the module
  finds it through its rpath, instead of `otswap/.dylibs/`, so that
  another package can link the same file (see `openmp_root()`).
- The README and the citation files say that a version of otswap is
  included in CosmoBolognaLib, and no longer ask to cite CosmoBolognaLib.

## 0.2.1 — 2026-10-28

### Added

- `io::readKeywords` (Python `io.read_keywords`) reads the keywords of a
  table, ASCII or FITS, from the first extension or one named by EXTNAME.
- `io::writeParts`, `io::readPart`, `io::hasPart` and `io::twinFile` (Python
  `io.write_parts`, `read_part`, `has_part`, `twin_file`, with `io.TablePart`)
  write and read files of several tables: FITS extensions named by EXTNAME, or
  ASCII twin files.

### Changed

- `reconstructBox` without box limits and with randoms given: `Result::box`
  is now the bounding box of the tracers and the given randoms, not of the
  tracers alone, so every random matched to a tracer lies inside it; also when
  mps is given. A computed mps follows that box, so such runs give a slightly
  different mps, and so different displacements, than in 0.2.0. Tracers that
  span no volume themselves, all in one plane for example, are accepted when
  the given randoms do. Drawn randoms and box limits are unchanged.
- Lengths are documented in coordinate units, the unit of the Cartesian
  coordinates: that of the tracers in a box, that of the distance table in a
  lightcone (Mpc/h from the cosmology constructor). The files otswap writes
  give length columns no unit (no TUNIT, no `[Mpc/h]`) and say "in coordinate
  units" in the comments of their length keywords; screen and error messages
  say "in coordinate units". The values are unchanged.
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
