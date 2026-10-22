# The executable

`otswap` runs a reconstruction, and optionally the correction of its
redshift-space distortions, from a parameter file, with no code to write. It
is built and installed with the C++ library (`bin/otswap` under the install
prefix); it is not in the Python package.

```
otswap <parameter-file>
otswap --help
otswap --version
```

`--help` prints the usage, the two options and the address of this page,
`--version` prints `otswap <version>`; both go to stdout. Any other argument
that starts with `--` is refused with the usage line; a parameter file whose
name starts with `--` is given as `./--name`. Without an argument, or with
more than one, it prints the usage line to stderr and exits with 1.

The parameter file has the syntax of {cpp:class}`otswap::io::ParameterFile`,
described on the [file formats](formats.md#parameter-files) page: one
`key = value` per line, `#` for comments, double quotes around a value that
holds spaces. Two commented files, one per geometry, are in
{repo-tree}`examples/executable`, at the end of this page.

## What it does

One run is one reconstruction, of a box or of a lightcone, on the tracers
and randoms read from the files given. It writes the files asked for: the
mean displacement of each tracer, every realization, and in a lightcone the
mean particle separation profile. With `perform_recons = true` it then moves
the tracers to real space, as {cpp:func}`otswap::realSpaceBox` or
{cpp:func}`otswap::realSpaceLightcone` do, and writes the corrected
catalogue.

There is one reconstruction only, on the observed tracers, in redshift
space: no second reconstruction on the corrected catalogue, and no void
search. With `perform_recons` the reconstruction's files therefore have
`zspace_` in their names, and the corrected catalogue is the only real-space
product.

The corrected catalogue has one row per tracer, in the order of the input,
with NaN coordinates for the tracers that were not moved. It cannot be given
back to otswap as `tracer_file`, since a reconstruction refuses non-finite
coordinates.

The two-step use with BitVF, otswap's complete displacement file
(`displacement_complete_zspace_<name>`) given to a one-pass BitVF run, finds
voids in redshift space only: voids in the corrected catalogue need BitVF's
own two-pass run.

Paths in the parameter file are used as written: relative ones resolve
against the current directory, not the parameter file's.

## Exit status and messages

The exit status is 0 when every file asked for is written, and 1 on any
error. Nothing is written to stdout, except by `--help` and `--version`.

Every error is one line on stderr, `otswap: error: <message>`. The run stops
at the first error, except for keys that are unknown or do not apply to the
run: those are all reported at once. Everything that can be checked without
the data is checked before any file is read, so a broken parameter file fails
at once rather than after a long run. The messages name the parameter file,
the line and the key:

```
otswap: error: run.ini:14: key 'Nrec': 0 is out of range; it must be at least 1
otswap: error: run.ini:9: key 'bins': does not apply to a box
otswap: error: run.ini: key 'units' is missing
```

Errors found by the library, in reading the data, the mask or the bias
table, or during the reconstruction or the correction, are printed with the
library's own message. The library's report, set by `verbosity`, also goes
to stderr, in lines that start with `otswap: ` and not with
`otswap: error: `; at `normal` and `detailed` the executable adds a line
`otswap: wrote <path>` for each file it writes.

## The keys

A key that is not listed here is an error (`unknown key`), and so is a key
that does not apply to the run (`does not apply to a box`, `to a lightcone`,
or `to a run without <switch>`): an applies-to column below marks with *B*
the keys of a box, with *L* those of a lightcone, and with a switch those
that apply only when it is on: `perform_recons = true`, `random_file` set, or
`reject_mask_crossings = true`. A key that is not required and not set takes
its default, which is the library's where there is one; a key with no
library default is required in the runs it applies to.

### Run

| Key | Type | Default | Applies to | Meaning |
|---|---|---|---|---|
| `geometry` | string | required | B, L | `box` or `lightcone` |
| `verbosity` | string | `normal` | B, L | `silent`, `normal` or `detailed`: what the reconstruction and the correction report |

### Input

| Key | Type | Default | Applies to | Meaning |
|---|---|---|---|---|
| `tracer_file` | string | required | B, L | the tracers, ASCII or FITS by the extension |
| `column_tracer` | 3 strings | required | B, L | box: x, y, z, in Mpc/h; lightcone: right ascension, declination, redshift. Each an index from 0, or a column name |
| `id_column` | string | none, no IDs | B, L | an integer column of the tracer file, read exactly as the tracers' IDs; it must not be one of `column_tracer`, and no ID may appear twice |
| `random_file` | string | box: none, the randoms are drawn uniformly in the box of `box_min` and `box_max`, or without it in the tracers' bounding box; lightcone: required | B, L | the randoms, as `tracer_file` |
| `column_random` | 3 strings | required | B with `random_file`, L | as `column_tracer` |
| `units` | string | required | L | `degrees` or `radians`: the angles of tracers and randoms |
| `delimiter` | character | `" "` | B, L | the column separator of every ASCII input, the bias file included; `" "` for a space, `\t` for a tab |
| `comment` | character | `#` | B, L | the character that starts a comment line in the same files |

`delimiter` and `comment` apply even when every input is FITS, where they
change nothing, so a parameter file stays valid when an input changes format.

### Lightcone selection

| Key | Type | Default | Applies to | Meaning |
|---|---|---|---|---|
| `veto_mask_filename` | string | required | L | the HEALPix FITS map of the observed sky; it also gives the sky area |
| `z_min` | double | none | L | the lower bound of the redshift cut, applied to tracers and randoms |
| `z_max` | double | none | L | the upper bound; not below `z_min` |

### Cosmology

| Key | Type | Default | Applies to | Meaning |
|---|---|---|---|---|
| `Omega_matter` | double | required | L, B with `perform_recons` | Ω_m of a flat w0-wa cosmology, in (0, 1] |
| `little_h` | double | required | L, B with `perform_recons` | H0 / (100 km/s/Mpc), positive |
| `w0` | double | -1 | L, B with `perform_recons` | w(a) = w0 + wa (1 - a) |
| `wa` | double | 0 | L, B with `perform_recons` | |

The cosmology gives a {cpp:class}`otswap::DistanceTable` over redshifts 0
to 10. In a lightcone it turns redshifts into distances; in a box it is used
only for the growth rate at `redshift`, so a box without `perform_recons`
takes none of these keys. A lightcone beyond redshift 10 fails with the
library's message.

### Reconstruction

| Key | Type | Default | Applies to | Meaning |
|---|---|---|---|---|
| `Nrec` | unsigned | 30 | B, L | the number of realizations, at least 1 |
| `threshold` | double | 1e-3 | B, L | the convergence threshold of the swap loop, strictly between 0 and 1 |
| `seed` | unsigned | 0 | B, L | the seed of the random streams; 0 draws one, written in the displacement files |
| `mesh_cellsize` | double | 4 | B, L | the grid cell size, in mean particle separations, positive; it changes the speed only |
| `bins` | unsigned | required | L | the redshift bins of the mean particle separation profile, at least 1 |
| `reject_mask_crossings` | bool | true | L | reject a displacement whose path crosses unobserved pixels |
| `max_unobserved_pixels_crossed` | unsigned | 0 | L with `reject_mask_crossings` | the distinct unobserved pixels a path may cross; 0 rejects at the first |
| `box_min` | 3 doubles | none | B | the lower x, y, z of the box, in Mpc/h; with `box_max` or not at all |
| `box_max` | 3 doubles | none | B | the upper x, y, z, each above `box_min`'s |
| `mps` | double | none, measured | B | the mean particle separation, in Mpc/h, positive |

Without `box_min` and `box_max` the box is the tracers' bounding box. The
[guide](guide.md) describes each parameter of the reconstruction.

### RSD correction

| Key | Type | Default | Applies to | Meaning |
|---|---|---|---|---|
| `perform_recons` | bool | false | B, L | correct the tracers for redshift-space distortions after the reconstruction |
| `rsd_smoothing_scale` | double | required | B, L with `perform_recons` | the width σ of the gaussian average of the displacements, in Mpc/h, non-negative; 0 for no average |
| `rsd_weight_by_realizations` | bool | false | B, L with `perform_recons` | weight each neighbour by its number of valid realizations |
| `redshift` | double | required | B with `perform_recons` | the box's redshift, in [0, 10] |
| `bias` | double | required | B with `perform_recons` | the tracers' linear bias, positive |
| `los_axis` | unsigned | required | B with `perform_recons` | the line of sight: 0, 1 or 2 for x, y or z |
| `bias_file` | string | required | L with `perform_recons` | the bias b(z): ASCII columns 0 and 1, redshift and bias, or a FITS table with `REDSHIFT` and `BIAS` |

10 Mpc/h is a reasonable starting value for `rsd_smoothing_scale`, but the
best value depends on the sample and should be checked in each analysis.

### Output

| Key | Type | Default | Applies to | Meaning |
|---|---|---|---|---|
| `output_dir` | string | required | B, L | the folder of the output files, created with its parents if missing |
| `output_file` | string | required | B, L | the `<name>` in the file names below: not empty, no `/` |
| `output_format` | string | `dat` | B, L | `dat` (ASCII) or `fits` |
| `printDisplMean` | bool | false | B, L | write the mean displacements |
| `printDisplComplete` | bool | false | B, L | write the complete displacements |
| `printMps` | bool | false | L | write the mean particle separation profile |
| `printRecCat` | bool | false | B, L with `perform_recons` | write the corrected catalogue; it must be true with `perform_recons` |

A run must write something: without `perform_recons`, one of the `print`
keys must be true. With `perform_recons`, `printRecCat = true` is required,
the corrected catalogue being the correction's only output.

## Output files

`<dir>` is `output_dir`, `<name>` `output_file`, and `<ext>` `dat` or `fits`
as `output_format` says. `zspace_` is in the names exactly when
`perform_recons` is true.

| Key | File | Written by |
|---|---|---|
| `printDisplMean` | `<dir>/displacement_mean_[zspace_]<name>.<ext>` | {cpp:func}`otswap::io::writeDisplacementMean` |
| `printDisplComplete` | `<dir>/displacement_complete_[zspace_]<name>.<ext>` | {cpp:func}`otswap::io::writeDisplacementComplete` |
| `printMps` | `<dir>/MPS_[zspace_]<name>.<ext>` | {cpp:func}`otswap::io::writeMpsProfile` |
| `printRecCat` | `<dir>/reconstructed_catalogue_<name>.<ext>` | {cpp:func}`otswap::io::writeRealSpaceCatalog` |

In ASCII, the complete displacement file of a lightcone comes with two more
files, the profile and the distance table a reader needs:
`displacement_complete_[zspace_]<name>.mps.dat` and
`displacement_complete_[zspace_]<name>.distances.dat`. In FITS they are
extensions of the same file.

The [file formats](formats.md) page gives the two layouts, ASCII (`dat`)
and FITS (`fits`), and the columns and keywords of each file. The reconstruction's files are written as soon as it
ends, before the correction. Existing files of the same names are written
over; no other file in `<dir>` is touched.

## The example parameter files

Both are run from the repository root, and write to
`examples/executable/output/` the corrected catalogue of the corresponding
example, `full_recons_box` or `full_recons_lightcone` (see
[examples](examples.md)), byte for byte. Every key that applies is listed;
those left at their default are commented out.

::::{tab-set}

:::{tab-item} Box

```{literalinclude} ../examples/executable/box.ini
:language: ini
```
:::

:::{tab-item} Lightcone

```{literalinclude} ../examples/executable/lightcone.ini
:language: ini
```
:::

::::
