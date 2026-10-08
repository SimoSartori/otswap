# Examples

Four examples, each a C++ program and a Python script that do the same thing:
the same inputs, parameters and seed, and the same output file, byte for
byte. Each folder in {repo-tree}`examples` is a CMake project of its own,
built against an installed otswap.

## The data

The examples read the files of {repo-tree}`examples/data`, which come with a
clone of the repository; they are not in the Python package or the source
distribution. They are provided for running the examples, may be
redistributed unchanged with them, and may not be used for any other purpose:
{repo-file}`examples/data/README.md` gives the terms and the content of each
file.

## Building and running

For an example `<name>`, with otswap installed under `<prefix>`:

```
cmake -S examples/<name> -B build/examples/<name> -DCMAKE_PREFIX_PATH=<prefix>
cmake --build build/examples/<name>
build/examples/<name>/<name> [data_dir]
python examples/<name>/<name>.py [data_dir]
```

Both read the examples' data folder unless `data_dir` is given, and write to
`examples/<name>/output/`, created if missing. Where cfitsio is not in a
system path, give its prefix with `-DCFITSIO_ROOT=...` (and `-DZLIB_ROOT=...`
for a static cfitsio).

## Box

Reconstructs the displacement field of halos in a cube, in real space,
with box geometry: Cartesian coordinates and a mean particle separation
computed from the tracers. Eight realizations, each on its own subset of the
randoms.

Reads:

- {repo-file}`examples/data/box_halos_real_space.dat`
- {repo-file}`examples/data/box_randoms.dat`

Writes `output/displacement_box.dat`, with {py:func}`otswap.io.write_displacements`; its columns
and keywords are on the [file formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/box/box.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/box/box.cpp
:language: cpp
```
:::

::::

## Lightcone

Reconstructs the displacement field of a survey-like catalogue with
lightcone geometry: sky coordinates, the mean particle separation measured as
a function of redshift, a redshift cut and the survey's HEALPix mask, with
the displacements that cross unobserved pixels rejected.

Reads:

- {repo-file}`examples/data/lightcone_tracers.dat`
- {repo-file}`examples/data/lightcone_randoms.dat`
- {repo-file}`examples/data/lightcone_mask.fits`

Writes `output/displacement_lightcone.dat`, with {py:func}`otswap.io.write_displacements`; its columns
and keywords are on the [file formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/lightcone/lightcone.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/lightcone/lightcone.cpp
:language: cpp
```
:::

::::

## Redshift-space correction in a box

Runs the box reconstruction on the same halos in redshift space, distorted
along one axis, then moves each halo along that axis to its real-space
position with the factor f/(b + 3f/5).

Reads:

- {repo-file}`examples/data/box_halos_redshift_space.dat`
- {repo-file}`examples/data/box_randoms.dat`

Writes `output/reconstructed_catalogue_box.dat`, with {py:func}`otswap.io.write_real_space_catalog`; its columns
and keywords are on the [file formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/full_recons_box/full_recons_box.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/full_recons_box/full_recons_box.cpp
:language: cpp
```
:::

::::

## Redshift-space correction in a lightcone

Runs the lightcone reconstruction, without the redshift cut, then moves each
tracer along its line of sight to its real-space position, with the growth
rate from the distance table and the bias from a b(z) table.

Reads:

- {repo-file}`examples/data/lightcone_tracers.dat`
- {repo-file}`examples/data/lightcone_randoms.dat`
- {repo-file}`examples/data/lightcone_mask.fits`
- {repo-file}`examples/data/lightcone_bias.dat`

Writes `output/reconstructed_catalogue_lightcone.dat`, with {py:func}`otswap.io.write_real_space_catalog`; its columns
and keywords are on the [file formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/full_recons_lightcone/full_recons_lightcone.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/full_recons_lightcone/full_recons_lightcone.cpp
:language: cpp
```
:::

::::
