# Names side by side

The two APIs are one library. A C++ name in camelCase is the Python name in
snake_case; a C++ configuration struct is a Python class of the same name,
passed as `config=` or `correction=`; a flat
`std::vector<double>` of `3 * N` entries becomes a NumPy array of shape
(N, 3). Python angles are in the `angle_unit` given, C++ angles in radians.

## Errors

Invalid input and failures raise `otswap::Error` in C++, `otswap.Error` in
Python, whose message names the object at fault. An error of the neighbour
search is one too: its message starts with `meshsearch: `, and in C++ the
original `meshsearch::Error` is nested in it (`std::rethrow_if_nested`
recovers it). The extrapolation of b(z) beyond its table is not an error:
Python issues an `ExtrapolationWarning`, C++ reports it and counts it in
`nExtrapolated`. A correction raises for the whole call, with no partial
catalogue, when one corrected position falls outside the distance table.

| C++ | Python | Notes |
|---|---|---|
| `Error` | `Error` | C++ records the throw location: `file()`, `function()`, `line()` |
| a line in the report, the count `nExtrapolated` | `ExtrapolationWarning` | b(z) extrapolated beyond its table |
| `meshsearch::Error` | — | rethrown as `Error`, the message starting `meshsearch: `; in C++ the original is nested |
| — | `KeyboardInterrupt` | Ctrl+C, at the next report line of a call |

## Configuration

| C++ | Python | Notes |
|---|---|---|
| `Config`, its fields | `Config`, `n_realizations`, `convergence`, `seed`, `cell_size`, `reject_crossings`, `max_unobserved_pixels_crossed`, `verbosity` | keyword `config` of `reconstruct_box`, `reconstruct_lightcone` and `reconstruct_lightcone_cartesian`, and of `Result.from_arrays` |
| `CorrectionConfig`, its fields | `CorrectionConfig`, `weight_by_realizations`, `verbosity` | keyword `correction` of `real_space_lightcone` and `real_space_box` |
| `Verbosity` | `"silent"`, `"normal"`, `"detailed"` | the `verbosity` field; the `verbosity` keyword of `measure_mps_profile` and `measure_mps_profile_cartesian` |
| `Config::log`, `CorrectionConfig::log`, the `log` of `measureMpsProfile` | — | Python prints to `sys.stdout` |
| `RedshiftCut`, `min`, `max` | `RedshiftCut`, `min`, `max`; keyword `redshift_cut` | also a pair `(min, max)` |
| `DistanceCut`, `min`, `max` | `DistanceCut`, `min`, `max`; keyword `distance_cut` | also a pair `(min, max)`; of the `_cartesian` functions only |
| `BoxLimits`, `lower`, `upper` | `BoxLimits`, `lower`, `upper`; keyword `box` | also a pair `((x0, y0, z0), (x1, y1, z1))` |
| `Geometry` | `"box"`, `"lightcone"` | the `geometry` property |
| `PixelOrdering` | keyword `nest` of `Mask.from_array` | |

## Distances and coordinates

| C++ | Python | Notes |
|---|---|---|
| `DistanceTable(OmegaM, h, w0, wa, zMin, zMax, nSamples)` | `DistanceTable.flat(omega_m, h, w0, wa, z_min, z_max, n_samples)` | |
| `DistanceTable(redshift, distance, growthRate)` | `DistanceTable.from_table(redshift, distance, growth_rate)` | |
| `distanceAt`, `redshiftAt`, `growthRateAt` | `distance_at`, `redshift_at`, `growth_rate_at` | Python takes arrays |
| `hasGrowthRate()`, `minRedshift()`, `maxRedshift()` | `has_growth_rate`, `min_redshift`, `max_redshift` | properties |
| `redshift()`, `distance()`, `growthRate()` | `redshift`, `distance`, `growth_rate` | read-only views; `growth_rate` None when the table carries none |
| `toCartesian`, `toSky` | `to_cartesian`, `to_sky` | |
| `pointToSky` | `point_to_sky` | Python takes an (N, 3) array and returns distance, right ascension, declination per row |
| `skyToRadians`, `skyToDegrees` | — | Python takes `angle_unit` instead |

## Mask

| C++ | Python | Notes |
|---|---|---|
| `Mask(fitsFile)` | `Mask(fits_file)` | |
| `Mask(values, ordering)`, `Mask(values, count, ordering)` | `Mask.from_array(values, nest)` | |
| `allows(ra, dec)` | `allows(ra, dec, *, angle_unit)` | Python takes arrays |
| `skyAreaDeg2()`, `nside()` | `sky_area_deg2`, `nside` | properties |

## Reconstruction

| C++ | Python | Notes |
|---|---|---|
| `reconstructBox`, two overloads | `reconstruct_box` | `mps` omitted to compute it; randoms omitted to draw them; `limits` is `box=`, `periodic` is `periodic=` |
| `reconstructLightcone`, sky overloads | `reconstruct_lightcone` | `sky_area_deg2=` or `mask=`; `tracers_sky`, `randoms_sky`, `distances` and `angle_unit` required |
| `reconstructLightcone`, Cartesian overloads | `reconstruct_lightcone_cartesian` | `sky_area_deg2=` or `mask=`; `tracers` and `randoms`; `distances=` optional; `distance_cut=`; `angle_unit` required, for the sky arrays of the result. In C++ the optional table comes after the `Config`, as a pointer, unlike the sky overloads |
| `measureMpsProfile`, sky overloads | `measure_mps_profile` | `sky_area_deg2=` or `mask=`; `verbosity` an argument, as no `Config` is taken |
| `measureMpsProfile`, Cartesian overloads | `measure_mps_profile_cartesian` | `tracers`; `distance_cut=`; no table, no angle unit |
| `rejectMaskCrossings` | `reject_mask_crossings` | |
| `recomputeMeans` | `recompute_means` | |

## Results

| C++ | Python | Notes |
|---|---|---|
| `Result` fields, camelCase | `Result` properties, snake_case | read-only views in Python |
| `Result::nRealizationsMean` | `Result.n_realizations_mean` | read-only in Python; set by `readDisplacementMean` |
| `Result::ids`, empty without IDs | `Result.ids`, 0..n_objects-1 without IDs | settable in Python, checked on assignment; a new copy on each access |
| `Result::config` | `Result.config` | a copy; `seed` is the one used |
| `Result::distances` | `Result.distances` | a copy; None in a box |
| — | `Result.angle_unit` | the unit of `tracers_sky` and `lagrangian_sky`; `"rad"` for a result made by another module |
| — | `Result.from_arrays` | a result from arrays read back from a file |
| `Result::mps`, NaN in a lightcone | `Result.mps`, None in a lightcone | |
| `Result::box`, empty in a lightcone | `Result.box`, None in a lightcone | a copy in Python |
| `Result::periodic` | `Result.periodic` | False in a lightcone |
| `Result::cellSide` | `Result.cell_side` | the cell side of the neighbour grids used; None when not made by a reconstruction |
| `Result::outsideBox` | `Result.outside_box` | a bool view in Python, as `outside_mask` |
| `Result::outsideDistanceCut` | `Result.outside_distance_cut` | a bool view in Python, as `outside_mask` |
| `Result::lagrangianDistance` | `Result.lagrangian_distance` | None in a box; in no file: `pointToSky` (`point_to_sky`) of the Lagrangian columns gives it |
| `SelectionCounts::boxLimits` | `SelectionCounts.box_limits` | None when no limits were given |
| `MpsProfile`, `at` | `MpsProfile`, `at` | Python's `at` takes arrays; `value` a redshift or a distance, as the profile |
| `SelectionCounts`, `message()` | `SelectionCounts`, `repr()` | |

## Redshift-space correction

| C++ | Python | Notes |
|---|---|---|
| `realSpaceLightcone`, `realSpaceBox` | `real_space_lightcone`, `real_space_box` | |
| `radialProjection`, `axisProjection` | `radial_projection`, `axis_projection` | |
| `neighbourAverage`, `NeighbourAverage` | `neighbour_average`, `NeighbourAverage` | a struct; a named tuple. The periodic overload is `box=` with `periodic=True` |
| `rsdFactor`, `rsdFactorBox` | `rsd_factor`, `rsd_factor_box` | |
| `shiftRadially`, `shiftAlongAxis` | `shift_radially`, `shift_along_axis` | |
| `BiasTable` | `BiasTable` | Python checks the table when it is built |
| `RealSpaceCatalog` | `RealSpaceCatalog` | `sky` in the `angle_unit` of the result corrected in Python, radians for a catalogue made by another module; `axis`, `box_redshift`, `box_bias` None in a lightcone |
| `RealSpaceCatalog::ids` | `RealSpaceCatalog.ids` | read-only in Python; a new copy on each access, as `Result.ids` |
| `CorrectionStatus` | `CorrectionStatus` | an `IntEnum`; `MovedByNeighbours` is `MOVED_BY_NEIGHBOURS` |

## Tables

| C++ | Python | Notes |
|---|---|---|
| `io::Column`, `io::Table` | `io.Column`, `io.Table` | |
| `io::read`, two overloads | `io.read` | `integer_columns=` for the second |
| `io::write`, `WriteOptions`, `Keyword` | `io.write`, keywords `precision`, `keywords` | keywords as `(name, value, comment)` tuples |
| `io::write` with a row function | — | |
| `io::readKeywords` | `io.read_keywords`, keyword `extname` | `(name, value, comment)` tuples; comment empty |
| `io::columnNames` | `io.column_names`, keywords `extname`, `delimiter`, `comment` | a list of str |
| `io::TablePart` | `io.TablePart` | columns from their data; no row function |
| `io::writeParts`, `readPart`, `hasPart`, `twinFile` | `io.write_parts`, `read_part`, `has_part`, `twin_file` | `extname`, `twin` keyword-only |
| `io::readBiasTable` | `io.read_bias_table` | |
| `io::writeDisplacementMean`, `writeDisplacementComplete`, `writeRealSpaceCatalog`, `writeMpsProfile` | `io.write_displacement_mean`, `write_displacement_complete`, `write_real_space_catalog`, `write_mps_profile` | |
| `io::readDisplacementMean`, `readDisplacementComplete` | `io.read_displacement_mean`, `read_displacement_complete` | `angle_unit=` required |
| `lightconeTables` of `writeDisplacementMean`, `writeDisplacementComplete` | `lightcone_tables` | keyword-only |
| `io::DisplacementGroup`, `io::CatalogGroup` | strings | `TracerSky` is `"tracer_sky"` |

## OpenMP runtime

| C++ | Python | Notes |
|---|---|---|
| — | `openmp_root` | the directory of the OpenMP runtime of a macOS wheel, or None; see [Using otswap from another package](packages.md#linking-the-same-openmp-runtime) |

## Random draws

| C++ | Python | Notes |
|---|---|---|
| `uniformInt`, `uniformReal`, `shuffle` | — | C++ only, in `otswap/Random.h`; the same values on every platform and compiler |
