# Guide

The concepts of otswap and their use, in the order a new user needs them.
Everything applies to both languages; where the names differ, the C++ name
comes first and the Python one follows in brackets, and [names side by
side](naming.md) lists them all. Each geometry has a numbered path of steps
and the [example](examples.md) that runs them.

## Tracers, randoms and the result

A reconstruction takes $N$ *tracers*, the observed objects, and a catalogue
of *randoms*, points that sample the same volume uniformly. It pairs each
tracer with one random so that the sum of the squared distances is small (the
[method](method.md) page describes how), and repeats the pairing on
`nRealizations` (`n_realizations`) disjoint sets of $N$ randoms. In each
realization a tracer's *displacement* points from the tracer to its random;
its mean over the realizations is the *mean displacement*, and the tracer
plus its mean displacement is its *mean Lagrangian position*. All of this is
returned in a {cpp:struct}`otswap::Result` ({py:class}`otswap.Result`), with
one row per input object, in input order.

The scale at which otswap works is the *mean particle separation* (mps),
$(V/N)^{1/3}$ for $N$ tracers in a volume $V$: the first guess pairs tracers
with randoms within 8 mps, and the neighbour grids have cells of
`cellSize` (`cell_size`) mps. In a box it is a single number; in a lightcone
it changes with redshift, and otswap measures it from the tracers.

There are four geometries, each described below: a box, a periodic box, a
lightcone on the sky and a lightcone in Cartesian coordinates.

## Coordinates and units

- Cartesian coordinates and every length (the mps, box limits,
  displacements, σ) are in coordinate units, as the layout below defines;
  otswap assumes no unit, and the files it writes give lengths none.
- Sky coordinates are ordered right ascension, declination, redshift.
- In C++ angles are in radians; {cpp:func}`otswap::skyToRadians` and
  {cpp:func}`otswap::skyToDegrees` convert a sky array. In Python every
  function that takes or returns angles has a required `angle_unit`, `"deg"`
  or `"rad"`.
- Right ascensions that otswap computes come back folded into [0, 2π), or
  [0, 360) degrees.
- In Python a coordinate array has shape (N, 3) and is converted to float64;
  the arrays of a result have the shapes below, with the realization first.

The layout of the C++ arrays, which the Python arrays map onto:

```{doxygengroup} coordinate_layout
```

## Box

**For a catalogue in comoving Cartesian coordinates whose density is about
constant**, cubic or not, such as a sub-volume of a simulation.

1. Read the tracers and, if there are any, the randoms ([Files](#files)).
2. Set the parameters in a {cpp:struct}`otswap::Config`
   ({py:class}`otswap.Config`); the defaults are the starting point for an
   analysis ([Choosing the parameters](#choosing-the-parameters)).
3. When the volume of the sample is known, give it as box limits,
   {cpp:struct}`otswap::BoxLimits` ({py:class}`otswap.BoxLimits`, or a pair of
   corners).
4. Reconstruct with {cpp:func}`otswap::reconstructBox`
   ({py:func}`otswap.reconstruct_box`).
5. Write the mean displacements, or every realization ([Files](#files)).

The [box example](examples.md#box) runs these steps with limits, and the
[correction in a box](examples.md#redshift-space-correction-in-a-box)
without them:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/box/box.py
:language: python
:dedent:
:start-at: "box = ([269.0, 269.0, 269.0], [669.0, 669.0, 669.0])"
:end-at: "result = otswap.reconstruct_box(tracers, randoms, box=box, config=config)"
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/box/box.cpp
:language: cpp
:dedent:
:start-at: "otswap::BoxLimits box;"
:end-at: "otswap::Result result = otswap::reconstructBox(tracers, randoms, config, box);"
```
:::

::::

**The box.** With limits, closed on every face, the tracers outside them take
no part: they keep their row, flagged in `outsideBox` (`outside_box`), with NaN
displacements, and the randoms given outside them are dropped. Without
limits, the box is the bounding box of the tracers and the given randoms, or
of the tracers alone when the randoms are drawn. Drawn randoms fill the box.
`Result::box` (`Result.box`) records the box used.

**The mps.** Unless it is given (`mps`), it is $(V/N)^{1/3}$, with $V$ the
volume of the box and $N$ the tracers inside it; `Result::mps`
(`Result.mps`) records the value used. Give it when the mean density of the
catalogue is known independently of the extent of its tracers and randoms.

**The faces.** The box is not periodic: its faces are boundaries through which
nothing flows. Modes on the scale of the box cannot be reconstructed, and the
displacements are least reliable near the faces. To use a region whole,
reconstruct a larger volume around it and discard a margin near the faces; for
the whole volume of a periodic simulation, use a periodic box.

## Periodic box

**For the whole volume of a simulation with periodic boundary conditions.**
A sub-volume cut from such a simulation is a box, not a periodic box.

1. Read the tracers. The randoms are usually drawn by otswap, uniformly in the
   box; they may also be given.
2. Give the box limits, which are the periods, and set `periodic`
   (`periodic=True`).
3. Reconstruct with {cpp:func}`otswap::reconstructBox`
   ({py:func}`otswap.reconstruct_box`).
4. Write the mean displacements, or every realization ([Files](#files)).

The [periodic box example](examples.md#periodic-box) runs these steps and
then corrects the halos for redshift-space distortions:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/periodic_box/periodic_box.py
:language: python
:dedent:
:start-at: "result = otswap.reconstruct_box(tracers, box="
:end-at: "periodic=True, config=config)"
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/periodic_box/periodic_box.cpp
:language: cpp
:dedent:
:start-at: "otswap::BoxLimits box;"
:end-at: "const otswap::Result result = otswap::reconstructBox("
```
:::

::::

Every distance is then the minimum image, and nothing is left out: tracers
and given randoms outside the limits are wrapped into the box by whole
periods and counted (`tracersWrapped`, `randomsWrapped` in `Result::selection`),
and the result holds the wrapped tracers. Each displacement is the
minimum-image vector, within half a period on each axis, and the mean
Lagrangian position is wrapped into the box. The mps, when computed, is that
of the volume of the limits and of every tracer; it rarely needs replacing.

**Each side must be a whole multiple of the shortest one**, as in a cube or
a cuboid of 500 × 500 × 1000: the neighbour grids use cubic cells that divide
every side. A box of other proportions can still be reconstructed without
`periodic`. The cells may therefore be smaller than `cellSize` mps;
`Result::cellSide` (`Result.cell_side`) gives the side used.
{cpp:func}`otswap::reconstructBox` states the rule in full.

## Lightcone on the sky

**For a survey catalogue in right ascension, declination and redshift.**

1. Build a distance table from the fiducial cosmology
   ([The distance table](#the-distance-table)).
2. Read the survey's mask ([The mask](#the-mask)), or know its sky area in
   square degrees.
3. Read the tracers and the randoms, which are required: they carry the
   survey's geometry, selection function and completeness, which otswap does
   not model ([Choosing the parameters](#choosing-the-parameters), "Randoms").
4. Choose the number of redshift bins of the mps, `nBins` (`n_bins`), with
   {cpp:func}`otswap::measureMpsProfile` ({py:func}`otswap.measure_mps_profile`).
5. To reconstruct a slice of the catalogue, set a redshift cut,
   {cpp:struct}`otswap::RedshiftCut` ({py:class}`otswap.RedshiftCut`, or a
   pair).
6. Reconstruct with {cpp:func}`otswap::reconstructLightcone`
   ({py:func}`otswap.reconstruct_lightcone`).
7. Write the mean displacements, or every realization ([Files](#files)).

The [lightcone example](examples.md#lightcone) runs these steps:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/lightcone/lightcone.py
:language: python
:dedent:
:start-at: "result = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask, n_bins=n_bins,"
:end-before: "#"
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/lightcone/lightcone.cpp
:language: cpp
:dedent:
:start-at: "const otswap::Result result = otswap::reconstructLightcone(tracersSky, randomsSky, mask, nBins,"
:end-before: "//"
```
:::

::::

The tracers and randoms outside the cut or on unobserved pixels of the mask
are left out; the tracers keep their row, flagged in `outsideRedshiftCut`
and `outsideMask` (`outside_redshift_cut`, `outside_mask`), with NaN
displacements. $N$ is the number of tracers kept. The mps is measured from
them in `nBins` uniform redshift bins over their range, and kept in the
result as a piecewise linear profile ({cpp:struct}`otswap::MpsProfile`,
{py:class}`otswap.MpsProfile`).

**Mask or sky area.** Given the sky area instead of a mask, otswap selects
nothing on the sky and rejects no displacement, so the tracers and randoms
must already be restricted to the footprint. The mask is the better choice
whenever one exists.

## Lightcone in Cartesian coordinates

**For a lightcone already in comoving Cartesian coordinates, with the
observer at the origin**, such as a mock catalogue built in a simulation.

1. Read the tracers and the randoms as Cartesian coordinates; the randoms are
   required, as on the sky.
2. Read the mask, or know the sky area. The mask is applied to the direction
   of each point, {cpp:func}`otswap::pointToSky` ({py:func}`otswap.point_to_sky`).
3. Choose the number of distance bins of the mps with the Cartesian
   overload of {cpp:func}`otswap::measureMpsProfile`
   ({py:func}`otswap.measure_mps_profile_cartesian`).
4. To reconstruct a shell of the catalogue, set a distance cut,
   {cpp:struct}`otswap::DistanceCut` ({py:class}`otswap.DistanceCut`).
5. Optionally, give a distance table, whose distances must be in the unit of
   the coordinates ([The distance table](#the-distance-table)).
6. Reconstruct with the Cartesian overload of
   {cpp:func}`otswap::reconstructLightcone`
   ({py:func}`otswap.reconstruct_lightcone_cartesian`).
7. Write the mean displacements, or every realization ([Files](#files)).

The [Cartesian lightcone example](examples.md#lightcone-in-cartesian-coordinates)
runs these steps, with the table, then again without it:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/lightcone_cartesian/lightcone_cartesian.py
:language: python
:dedent:
:start-at: "result = otswap.reconstruct_lightcone_cartesian(tracers, randoms, mask=mask, n_bins=n_bins,"
:end-before: "#"
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/lightcone_cartesian/lightcone_cartesian.cpp
:language: cpp
:dedent:
:start-at: "const otswap::Result result = otswap::reconstructLightcone(tracers, randoms, mask, nBins, config, &distances, cut);"
:end-at: "const otswap::Result result = otswap::reconstructLightcone(tracers, randoms, mask, nBins, config, &distances, cut);"
```
:::

::::

The reconstruction needs no cosmology. The mps is measured in bins of
comoving distance from the origin, a profile in distance. The distance cut
tests the distance of the given positions, tracers and randoms, before the
reconstruction: the tracers outside it keep their row, flagged in
`outsideDistanceCut` (`outside_distance_cut`), and the randoms outside it are
dropped; a mean Lagrangian position may lie outside it. The distance table
plays no part in the reconstruction: it only fills the redshifts of the
result, which are NaN without it. In C++ it is passed by pointer, after the
`Config`, since it is optional.

## The distance table

A {cpp:class}`otswap::DistanceTable` ({py:class}`otswap.DistanceTable`) holds
the comoving distance and the linear growth rate $f$ on a grid of redshifts,
and interpolates linearly between its nodes. It is built in one of two ways:

- from a flat w0-wa cosmology, $\Omega_m$ and $h$ with $w_0 = -1$ and
  $w_a = 0$ by default, over redshifts 0 to 10 by default, with distances in
  Mpc/h: `DistanceTable(OmegaM, h, ...)` ({py:meth}`otswap.DistanceTable.flat`);
- from the caller's own columns of redshift, distance and, optionally,
  growth rate: `DistanceTable(redshift, distance, growthRate)`
  ({py:meth}`otswap.DistanceTable.from_table`), for instance to give a growth
  rate of one's own.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/lightcone/lightcone.py
:language: python
:dedent:
:start-at: "distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)"
:end-at: "distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)"
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/lightcone/lightcone.cpp
:language: cpp
:dedent:
:start-at: "const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);"
:end-at: "const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);"
```
:::

::::

**The table is never extrapolated.** The redshifts of the tracers and randoms
kept, and the corrected positions of the redshift-space correction, must lie
in its range; otherwise the call raises, naming the row. A mean Lagrangian
position beyond the table gets a NaN redshift, without error.

**Its distances must be in the coordinate units.** For a lightcone on the sky
they define them (Mpc/h from a cosmology). For a lightcone in Cartesian
coordinates they must be in the unit of the coordinates: a table in Mpc/h with
coordinates in Mpc gives wrong redshifts and growth rates, and nothing can
detect it. The redshift-space correction must be given the table the result
was reconstructed with, or one with the same distances.

## Realizations and validity

The reconstruction is repeated `nRealizations` (`n_realizations`) times, each
on its own set of $N$ randoms, disjoint from the others. Each displacement has
a validity flag, `valid`, per realization and object; every flag of a tracer
that took part is set after the reconstruction, and the mask filter only ever
clears them. For each tracer, `validRealizations` (`valid_realizations`)
counts the valid ones, and `meanDisplacement` (`mean_displacement`) is their
mean, NaN when there is none. `lagrangian`, the tracer plus its mean
displacement, follows from it; the [method](method.md#realizations) page
explains what it approximates.

After changing `valid` by hand, {cpp:func}`otswap::recomputeMeans`
({py:func}`otswap.recompute_means`) recomputes the fields that follow from it,
as the reconstruction computes them.

Every output keeps one row per input object, in input order: a tracer left
out by a box, a cut or a mask keeps its row, flagged, with NaN displacements
and no valid realization.

## The mask

The mask describes the survey's footprint on the sky, so that a lightcone
reconstruction uses only what was observed and discards the displacements that
pass through what was not. It is a full-sky HEALPix map:

- from a FITS file, {cpp:class}`otswap::Mask` ({py:class}`otswap.Mask`): a
  binary table in HDU 2 whose first column holds the values, 1024 to a row,
  so that NSIDE is a multiple of 16. The ordering, RING or NESTED, is read
  from the header. A cut-sky map is refused;
- from an array in memory, `Mask(values, ordering)`
  ({py:meth}`otswap.Mask.from_array`, with `nest`). The ordering is the
  caller's to give, and a wrong one gives a wrong footprint without error.

A pixel is observed when its value is greater than 0. The value is not a
weight: a fractional map is turned into observed and unobserved pixels, and
the completeness of the survey belongs in the randoms. Given to the
lightcone reconstruction, the mask is applied twice:

- before it, tracers and randoms on unobserved pixels are left out;
- after it, a displacement whose great-circle arc crosses more than
  `maxUnobservedPixelsCrossed` (`max_unobserved_pixels_crossed`) unobserved
  pixels in a row is marked invalid; 0, the default, rejects it at the first
  one. {cpp:func}`otswap::rejectMaskCrossings`
  ({py:func}`otswap.reject_mask_crossings`) gives the rule.

The count is of pixels, not an angle, so the same threshold means a different
criterion at a different NSIDE. A second filter with a higher threshold
changes nothing, and one with a lower threshold rejects more: set the
threshold in the call that reconstructs, or set `rejectCrossings`
(`reject_crossings`) to false and filter afterwards with
{cpp:func}`otswap::rejectMaskCrossings`.

## The redshift-space correction

The displacements of a reconstruction in redshift space move the catalogue
to real space: {cpp:func}`otswap::realSpaceBox` and
{cpp:func}`otswap::realSpaceLightcone` ({py:func}`otswap.real_space_box`,
{py:func}`otswap.real_space_lightcone`) shift each tracer along its line of
sight by the line-of-sight component of its mean displacement, averaged over
its neighbours with a Gaussian of width σ and scaled by $f/(b + 3f/5)$. The
[method](method.md#redshift-space-correction) page derives the factor. The
inputs are:

- **the result**, of a reconstruction of the observed catalogue;
- **the line of sight**: a Cartesian axis in a box, 0, 1 or 2; in a
  lightcone, the direction of each tracer from the origin;
- **the growth rate**, from a distance table carrying it: in a box at the
  box's redshift, in a lightcone at each tracer's redshift. A result of a
  Cartesian lightcone reconstructed without a table can still be corrected:
  the table given to the correction supplies the redshifts;
- **the bias** of the tracers, the caller's: one number for a box; for a
  lightcone a b(z) table ({cpp:struct}`otswap::BiasTable`,
  {py:class}`otswap.BiasTable`, or read with `io::readBiasTable`,
  {py:func}`otswap.io.read_bias_table`), extrapolated linearly beyond its
  nodes with a warning;
- **σ**, in coordinate units; 0 applies no average
  ([Choosing the parameters](#choosing-the-parameters)).

A periodic result needs nothing more: the neighbours are found across the
faces, and the corrected coordinate is wrapped into the box. In a lightcone, a
corrected position outside the distance table raises for the whole call; a
table covering a wider redshift range solves it.

The [two correction examples](examples.md#redshift-space-correction-in-a-box)
correct a box and a lightcone, the [periodic box
example](examples.md#periodic-box) a periodic box and the [Cartesian
lightcone example](examples.md#lightcone-in-cartesian-coordinates) a
Cartesian lightcone. The result, an {cpp:struct}`otswap::RealSpaceCatalog`
({py:class}`otswap.RealSpaceCatalog`), gives each tracer its corrected
position and a status ({cpp:enum}`otswap::CorrectionStatus`,
{py:class}`otswap.CorrectionStatus`): moved with its own displacement; moved
with its neighbours' average, when it has no valid realization; not moved,
when no valid tracer lies within 3σ; or left out of the reconstruction.

## Files

Tables are read and written by `otswap::io` (`otswap.io`), in ASCII or FITS
by the extension, and a file written from Python is the same, byte for byte,
as the one C++ writes from the same values. The [file formats](formats.md)
page gives the layouts, and the columns and keywords of every file.

- **Inputs**: {cpp:func}`otswap::io::read` ({py:func}`otswap.io.read`) reads
  columns by index or by name; `Table::values` (`Table.values`) are the
  coordinates.
- **Displacements**: `writeDisplacementMean` (`write_displacement_mean`)
  writes the mean displacement of each tracer; `writeDisplacementComplete`
  (`write_displacement_complete`) writes every realization, without loss.
  Each has its reader, which rebuilds a `Result`.
- **The corrected catalogue**: `writeRealSpaceCatalog`
  (`write_real_space_catalog`).
- **The mps profile**: `writeMpsProfile` (`write_mps_profile`).

**IDs.** To carry the catalogue's own IDs into the files, read them exactly as
integers (the `integerColumns` overload of `io::read`, `integer_columns=` in
Python) and set them on the result before writing or correcting it:
`Result::ids` (`result.ids = ...`), one 64-bit integer per object, no value
twice. The writers then add an `id` column, and the corrected catalogue
copies them. IDs equal to 0, 1, ..., N − 1 count as none; an `id` column,
where written without IDs, holds each object's row, from 0.

## Choosing the parameters

The defaults are the starting point for an analysis; the examples use smaller
values to run quickly.

- **Randoms.** Each realization takes $N$ randoms, disjoint from those of the
  others, so the catalogue must hold at least `nRealizations` × $N$; a surplus
  is not used, the randoms of each realization being drawn at random from the
  catalogue. $N$ is the number of tracers reconstructed: those inside the box
  limits, or kept by the cut and the mask, the randoms being counted after the
  same selection. In a box, randoms may also be drawn by otswap, uniformly in
  the box; given randoms are needed when the volume is not filled uniformly.
  The default of 30 realizations needs 30 randoms per tracer; the example
  catalogues hold about 8, which is why the examples run 8 realizations.
- **`convergence`**, default 1e-3: the sweeps stop once the fraction of
  successful swaps in a sweep, swaps per tracer visited, no longer exceeds
  this threshold. In the void analysis of Sartori et al. (2026, Appendix B),
  void abundances were stable for thresholds of about 1e-2 or below, and
  stacked void profiles for 1e-3 or below; larger values are faster.
- **`nRealizations`** (`n_realizations`), default 30: that analysis used 10 to
  50, depending on the geometry, the small scales needing the most.
- **`nBins`** (`n_bins`), lightcone only: the bins are uniform over the range
  of the tracers kept, and each must hold at least 278 tracers, for a 2%
  accuracy on its mps (see {cpp:struct}`otswap::MpsProfile`). Take enough
  bins to follow the change of the density with redshift or distance. Choose
  them before reconstructing: call `measureMpsProfile`
  ({py:func}`otswap.measure_mps_profile`,
  {py:func}`otswap.measure_mps_profile_cartesian`) with a trial number and read
  `count`. A refused number comes with an upper bound, the number of tracers
  divided by 278; a sparse tail of the distribution may need fewer bins.
- **`maxUnobservedPixelsCrossed`** (`max_unobserved_pixels_crossed`),
  default 0: 0 rejects any crossing; to tolerate small holes such as
  bright-star masks, set the size of those holes in pixels at the mask's
  NSIDE.
- **`mps`**, box only: computed by default ([Box](#box)); give it when the
  mean density of the catalogue is known independently of the extent of its
  tracers and randoms. In a periodic box the computed value rarely needs
  replacing.
- **`periodic`**, box only: for the whole volume of a simulation with
  periodic boundary conditions ([Periodic box](#periodic-box)).
- **Redshift or distance cut**: to reconstruct a slice of a wider catalogue.
- **`cellSize`** (`cell_size`), default 3 mps: affects the speed only, never
  the result.
- **σ of the redshift-space correction:** for coordinates in Mpc/h, 10 Mpc/h
  is a starting value, to be calibrated on the sample.
- **`weightByRealizations`** (`weight_by_realizations`), default false:
  set it to give less weight to tracers whose displacement averages few valid
  realizations, as near the edges of a masked footprint. It changes nothing
  where every tracer has all its realizations, as in a box.

## Reproducibility and threads

```{include} ../README.md
:start-after: <!-- docs:reproducibility-start -->
:end-before: <!-- docs:reproducibility-end -->
```

The seed used is `Result::config.seed` (`Result.config.seed`). The number of
threads follows `OMP_NUM_THREADS`. The reconstructions, the mask filter, the
redshift-space correction and the table readers and writers of `otswap.io`
release the GIL while they run. A `Result` must not be read from one thread
while another modifies it in place, as {cpp:func}`otswap::rejectMaskCrossings`
and {cpp:func}`otswap::recomputeMeans` ({py:func}`otswap.reject_mask_crossings`,
{py:func}`otswap.recompute_means`) do.

## What it prints

Each call reports what it did, in lines starting with `otswap:`, at the level
`verbosity` sets: silent, normal (the default: what was selected and lost)
or detailed ({cpp:enum}`otswap::Verbosity`, {py:data}`otswap.Verbosity`, which
list the lines). The periodic box example prints, at the normal level:

```text
otswap: wrapped 0 of 50000 tracers into the periodic box [0, 1000] x [0, 1000] x [0, 1000]
otswap: reconstructBox: 8 realizations of 50000 tracers in 0.88 s; 0 without a valid displacement
otswap: realSpaceBox: corrected 50000 of 50000 tracers in 0.01 s; 0 left uncorrected
```

C++ writes to the stream set in `Config::log`, `CorrectionConfig::log` or the
`log` argument of `measureMpsProfile`, `std::clog` by default; Python to
`sys.stdout`, which a notebook shows. In Python, Ctrl+C stops a reconstruction,
`measure_mps_profile` or a correction at its next line with
`KeyboardInterrupt`.

## Errors

Invalid input and failures raise {cpp:class}`otswap::Error`
({py:class}`otswap.Error`), whose message names the object at fault;
[names side by side](naming.md#errors) describes them in both languages.
Among the conditions every reconstruction checks: at least 114 tracers
reconstructed (each needs 113 distinct neighbours), and enough randoms.
