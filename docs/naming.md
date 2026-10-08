# Names side by side

The two APIs are one library. A C++ name in camelCase is the Python name in
snake_case; a C++ configuration struct becomes keyword arguments; a flat
`std::vector<double>` of `3 * N` entries becomes a NumPy array of shape
(N, 3). Python angles are in the `angle_unit` given, C++ angles in radians.

## Errors

| C++ | Python | Notes |
|---|---|---|
| `Error` | `Error` | C++ records the throw location: `file()`, `function()`, `line()` |
| a line in the report, the count `nExtrapolated` | `ExtrapolationWarning` | b(z) extrapolated beyond its table |
| `meshsearch::Error` | `RuntimeError` | positions that span no volume, in the neighbour search |

## Configuration

| C++ | Python | Notes |
|---|---|---|
| `Config::nRealizations`, `convergence`, `seed`, `cellSize` | keywords `n_realizations`, `convergence`, `seed`, `cell_size` | of `reconstruct_box` and `reconstruct_lightcone` |
| `Config::rejectCrossings`, `maxUnobservedPixelsCrossed` | keywords `reject_crossings`, `max_unobserved_pixels_crossed` | of `reconstruct_lightcone` |
| `Config::verbosity`, `CorrectionConfig::verbosity` | keyword `verbosity` | `Verbosity::Detailed` is `"detailed"` |
| `Config::log`, `CorrectionConfig::log` | — | Python prints to `sys.stdout` |
| `CorrectionConfig::weightByRealizations` | keyword `weight_by_realizations` | |
| `RedshiftCut` | keyword `redshift_cut`, a `(min, max)` tuple | |
| `Geometry` | `"box"`, `"lightcone"` | the `geometry` property |
| `PixelOrdering` | keyword `nest` of `Mask.from_array` | |

## Distances and coordinates

| C++ | Python | Notes |
|---|---|---|
| `DistanceTable(OmegaM, h, w0, wa, zMin, zMax, nSamples)` | `DistanceTable.flat(omega_m, h, w0, wa, z_min, z_max, n_samples)` | |
| `DistanceTable(redshift, distance, growthRate)` | `DistanceTable.from_table(redshift, distance, growth_rate)` | |
| `distanceAt`, `redshiftAt`, `growthRateAt` | `distance_at`, `redshift_at`, `growth_rate_at` | Python takes arrays |
| `hasGrowthRate()`, `minRedshift()`, `maxRedshift()` | `has_growth_rate`, `min_redshift`, `max_redshift` | properties |
| `toCartesian`, `toSky` | `to_cartesian`, `to_sky` | |
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
| `reconstructBox`, two overloads | `reconstruct_box` | `mps` omitted to compute it; randoms omitted to draw them |
| `reconstructLightcone`, four overloads | `reconstruct_lightcone` | `sky_area_deg2=` or `mask=`; `tracers=` and `randoms=` for the Cartesian overloads |
| `rejectMaskCrossings` | `reject_mask_crossings` | |
| `recomputeMeans` | `recompute_means` | |

## Results

| C++ | Python | Notes |
|---|---|---|
| `Result` fields, camelCase | `Result` properties, snake_case | read-only views in Python |
| `Result::config` | `Result.seed` | the seed used; the other parameters are those of the call |
| `Result::distances` | — | kept inside the Python result |
| — | `Result.angle_unit` | the unit of `tracers_sky` and `lagrangian_sky` |
| — | `Result.from_arrays` | a result from arrays read back from a file |
| `Result::mps`, NaN in a lightcone | `Result.mps`, None in a lightcone | |
| `MpsProfile`, `at` | `MpsProfile`, `at` | Python's `at` takes arrays |
| `SelectionCounts`, `message()` | `SelectionCounts`, `repr()` | |

## Redshift-space correction

| C++ | Python | Notes |
|---|---|---|
| `realSpaceLightcone`, `realSpaceBox` | `real_space_lightcone`, `real_space_box` | |
| `radialProjection`, `axisProjection` | `radial_projection`, `axis_projection` | |
| `neighbourAverage`, `NeighbourAverage` | `neighbour_average`, `NeighbourAverage` | a struct; a named tuple |
| `rsdFactor`, `rsdFactorBox` | `rsd_factor`, `rsd_factor_box` | |
| `shiftRadially`, `shiftAlongAxis` | `shift_radially`, `shift_along_axis` | |
| `BiasTable` | `BiasTable` | Python checks the table when it is built |
| `RealSpaceCatalog` | `RealSpaceCatalog` | `sky` in the result's `angle_unit` in Python; `axis`, `box_redshift`, `box_bias` None in a lightcone |
| `CorrectionStatus` | `CorrectionStatus` | an `IntEnum`; `MovedByNeighbours` is `MOVED_BY_NEIGHBOURS` |

## Tables

| C++ | Python | Notes |
|---|---|---|
| `io::Column`, `io::Table` | `io.Column`, `io.Table` | |
| `io::read`, two overloads | `io.read` | `integer_columns=` for the second |
| `io::write`, `WriteOptions`, `Keyword` | `io.write`, keywords `precision`, `keywords` | keywords as `(name, value, comment)` tuples |
| `io::write` with a row function | — | |
| `io::readBiasTable` | `io.read_bias_table` | |
| `io::writeDisplacements`, `writeDisplacementField`, `writeRealSpaceCatalog`, `writeMpsProfile` | `io.write_displacements`, `write_displacement_field`, `write_real_space_catalog`, `write_mps_profile` | |
| `io::DisplacementGroup`, `io::CatalogGroup` | strings | `TracerSky` is `"tracer_sky"` |
