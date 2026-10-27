# Names side by side

The two APIs are one library. A C++ name in camelCase is the Python name in
snake_case; a C++ configuration struct is a Python class of the same name,
passed as `config=` or `correction=`; a flat
`std::vector<double>` of `3 * N` entries becomes a NumPy array of shape
(N, 3). Python angles are in the `angle_unit` given, C++ angles in radians.

## Errors

| C++ | Python | Notes |
|---|---|---|
| `Error` | `Error` | C++ records the throw location: `file()`, `function()`, `line()` |
| a line in the report, the count `nExtrapolated` | `ExtrapolationWarning` | b(z) extrapolated beyond its table |
| `meshsearch::Error` | — | rethrown as `Error`, the message starting `meshsearch: `; in C++ the original is nested |

## Configuration

| C++ | Python | Notes |
|---|---|---|
| `Config`, its fields | `Config`, `n_realizations`, `convergence`, `seed`, `cell_size`, `reject_crossings`, `max_unobserved_pixels_crossed`, `verbosity` | keyword `config` of `reconstruct_box` and `reconstruct_lightcone`, and of `Result.from_arrays` |
| `CorrectionConfig`, its fields | `CorrectionConfig`, `weight_by_realizations`, `verbosity` | keyword `correction` of `real_space_lightcone` and `real_space_box` |
| `Verbosity` | `"silent"`, `"normal"`, `"detailed"` | the `verbosity` field; the `verbosity` keyword of `measure_mps_profile` |
| `Config::log`, `CorrectionConfig::log`, the `log` of `measureMpsProfile` | — | Python prints to `sys.stdout` |
| `RedshiftCut`, `min`, `max` | `RedshiftCut`, `min`, `max`; keyword `redshift_cut` | also a pair `(min, max)` |
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
| `reconstructBox`, two overloads | `reconstruct_box` | `mps` omitted to compute it; randoms omitted to draw them; `limits` is `box=` |
| `reconstructLightcone`, four overloads | `reconstruct_lightcone` | `sky_area_deg2=` or `mask=`; `tracers=` and `randoms=` for the Cartesian overloads |
| `measureMpsProfile`, two overloads | `measure_mps_profile` | `sky_area_deg2=` or `mask=`; `verbosity` an argument, as no `Config` is taken |
| `rejectMaskCrossings` | `reject_mask_crossings` | |
| `recomputeMeans` | `recompute_means` | |

## Results

| C++ | Python | Notes |
|---|---|---|
| `Result` fields, camelCase | `Result` properties, snake_case | read-only views in Python |
| `Result::ids`, empty without IDs | `Result.ids`, 0..n_objects-1 without IDs | settable in Python, checked on assignment; a new copy on each access |
| `Result::config` | `Result.config` | a copy; `seed` is the one used |
| `Result::distances` | `Result.distances` | a copy; None in a box |
| — | `Result.angle_unit` | the unit of `tracers_sky` and `lagrangian_sky`; `"rad"` for a result made by another module |
| — | `Result.from_arrays` | a result from arrays read back from a file |
| `Result::mps`, NaN in a lightcone | `Result.mps`, None in a lightcone | |
| `Result::box`, empty in a lightcone | `Result.box`, None in a lightcone | a copy in Python |
| `Result::outsideBox` | `Result.outside_box` | a bool view in Python, as `outside_mask` |
| `SelectionCounts::boxLimits` | `SelectionCounts.box_limits` | None when no limits were given |
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
| `io::readBiasTable` | `io.read_bias_table` | |
| `io::writeDisplacementMean`, `writeDisplacementComplete`, `writeRealSpaceCatalog`, `writeMpsProfile` | `io.write_displacement_mean`, `write_displacement_complete`, `write_real_space_catalog`, `write_mps_profile` | |
| `io::readDisplacementComplete` | `io.read_displacement_complete` | `angle_unit=` required |
| `io::DisplacementGroup`, `io::CatalogGroup` | strings | `TracerSky` is `"tracer_sky"` |

## Use from another nanobind module

The Python types are the C++ structs themselves, so another nanobind module
can take and return them, naming the C++ types in its own signatures:
`Config`, `CorrectionConfig`, `RedshiftCut`, `BoxLimits`, `Mask`,
`DistanceTable`, `MpsProfile`, `BiasTable`, `Result` and `RealSpaceCatalog`.
For example, `double f (const otswap::Result&)` takes an `otswap.Result`, and
a function returning an `otswap::Config` returns an `otswap.Config`. The
module binds none of these types itself.

- Both modules are built with nanobind 3.1.0 and the same nanobind ABI (the
  same compiler family and C++ standard library), and neither defines
  `NB_DOMAIN`. `import otswap` comes before the other module's functions are
  used: it registers the types. The pairs `(min, max)` and
  `((x0, y0, z0), (x1, y1, z1))` are converted in the other module's
  functions too.
- The layouts of the structs are part of the interface across 0.2.x: a field
  added in a patch release would break a module built against the older
  header.
- A `Config` or `CorrectionConfig` received from Python holds
  `log = &std::clog`; the other module sets it before using it.
- A `Result` or `RealSpaceCatalog` made by the other module is in radians:
  its `angle_unit` is `"rad"`. The unit of one made by otswap from angles in
  degrees is kept by otswap, not in the object.
- Take a `Result` as `const otswap::Result&`. One changed in place behind
  otswap's back keeps the degree arrays already read from it; to filter in
  place, call otswap's own functions.
- An `otswap::Error` thrown by otswap's Python functions reaches Python as
  `otswap.Error`. One thrown where the other module calls otswap in C++ is
  that module's to translate: otswap's translator does not recognise it
  with every standard library (with libc++, two modules that each link
  otswap do not share the type information of `otswap::Error`), and it then
  arrives as a plain `RuntimeError`. The recommended way is to catch it
  where the module calls otswap and rethrow the module's own error, its
  message prefixed by `"otswap: "`, as otswap does with the errors of
  meshsearch. For a malformed `Result` or catalogue, which only another
  module can make, the library's consistency messages name the C++ members
  (`meanDisplacement`), and the array properties raise `otswap.Error` naming
  the Python ones (`mean_displacement`).
