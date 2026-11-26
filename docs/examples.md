# Examples

Six examples, each a C++ program and a Python script that do the same thing:
the same inputs, parameters and seed, and the same output files, byte for
byte. Each example folder in {repo-tree}`examples` is a CMake project of its
own, built against an installed otswap.

## The data

The examples read the files of {repo-tree}`examples/data`, which come with a
clone of the repository; they are not in the Python package or the source
distribution. They are provided for running the examples, may be
redistributed unchanged with them, and may not be used for any other purpose:
{repo-file}`examples/data/README.md` gives the terms and the content of each
file. Lengths are in Mpc/h (the boxes are cubes of side 400 and 1000 Mpc/h)
and angles in degrees.

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
with box geometry: Cartesian coordinates, the box limits, and a mean
particle separation computed from the volume of the limits and the number of
tracers. Eight realizations, each on its own subset of the randoms. It then
sets IDs on the result, writes every realization to a FITS file, reads it
back, and compares the mean over half of the realizations with the mean over
all of them.

Reads:

- {repo-file}`examples/data/box_halos_real_space.dat`
- {repo-file}`examples/data/box_randoms.dat`

Writes `output/displacement_box.dat`, with {py:func}`otswap.io.write_displacement_mean`, and
`output/displacement_complete_box.fits`, with {py:func}`otswap.io.write_displacement_complete`; their
columns and keywords are on the [file formats](formats.md) page.

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
a function of redshift, a redshift cut that keeps a slice of the catalogue,
and the survey's HEALPix mask, with the displacements that cross unobserved
pixels rejected. Before the reconstruction it measures the mean particle
separation alone, from the mask and from the sky area; after it, it reads
the file written back.

Reads:

- {repo-file}`examples/data/lightcone_tracers.dat`
- {repo-file}`examples/data/lightcone_randoms.dat`
- {repo-file}`examples/data/lightcone_mask.fits`

Writes `output/mps_lightcone.dat`, with {py:func}`otswap.io.write_mps_profile`;
`output/displacement_lightcone.dat`, with {py:func}`otswap.io.write_displacement_mean`, and its
twin files `displacement_lightcone.mps.dat` and `displacement_lightcone.distances.dat`, the mps(z)
profile and the distance table; their columns and keywords are on the [file formats](formats.md) page.

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

Runs the box reconstruction, without box limits, on the same halos in
redshift space, distorted along the y axis, then moves each halo along that
axis to its real-space position with the factor f/(b + 3f/5). It reads the
corrected catalogue back by column name and compares it with the halos in
real space.

Reads:

- {repo-file}`examples/data/box_halos_redshift_space.dat`
- {repo-file}`examples/data/box_randoms.dat`
- {repo-file}`examples/data/box_halos_real_space.dat`, for the comparison

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

## Periodic box

Reconstructs the displacement field of halos in redshift space in a cube of a
simulation, with periodic boundary conditions: the box is given by its
limits, every distance is the minimum image, and nothing is cut at the
faces. No randoms are given: otswap draws them in the box. The result
reports the cell side of the neighbour grids actually used, which in a
periodic box divides the shortest side. The halos are then corrected to real
space, with neighbours across the faces and the corrected coordinate wrapped
into the box, and compared with their real-space positions, which the data
also give, along the line of sight and by the minimum image.

Reads:

- {repo-file}`examples/data/periodic_halos.dat`

Writes `output/displacement_periodic_box.dat`, with {py:func}`otswap.io.write_displacement_mean`, and
`output/reconstructed_catalogue_periodic_box.dat`, with {py:func}`otswap.io.write_real_space_catalog`; their
columns and keywords, those of the periodic box included, are on the [file formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/periodic_box/periodic_box.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/periodic_box/periodic_box.cpp
:language: cpp
```
:::

::::

## Lightcone in Cartesian coordinates

Reconstructs and corrects the lightcone of the two lightcone examples given
in comoving Cartesian coordinates, the observer at the origin: the mean
particle separation is measured in distance bins, a cut on the comoving
distance keeps a shell of the catalogue, and the mask acts on the angles of
the objects. It runs the steps of the two sky examples in one run, then the
reconstruction again without the distance table, which gives the same
displacements.

Reads:

- {repo-file}`examples/data/lightcone_tracers_cartesian.dat`
- {repo-file}`examples/data/lightcone_randoms_cartesian.dat`
- {repo-file}`examples/data/lightcone_mask.fits`
- {repo-file}`examples/data/lightcone_bias.dat`

Writes `output/displacement_lightcone_cartesian.dat`, with {py:func}`otswap.io.write_displacement_mean`, and its
twin files `.mps.dat` and `.distances.dat`, the mps(d) profile and the distance table;
`output/displacement_lightcone_cartesian_no_table.dat`, the result without the table, and its twin
`.mps.dat`; and `output/reconstructed_catalogue_lightcone_cartesian.dat`, with {py:func}`otswap.io.write_real_space_catalog`; their columns
and keywords are on the [file formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/lightcone_cartesian/lightcone_cartesian.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/lightcone_cartesian/lightcone_cartesian.cpp
:language: cpp
```
:::

::::
