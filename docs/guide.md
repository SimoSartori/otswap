# Guide

What applies to both APIs. The names are the C++ ones and the Python ones
side by side; [names side by side](naming.md) maps them all.

## Two geometries

**Box.** Cartesian coordinates, and a single mean particle separation (mps),
given or computed from the tracers as $(V/N)^{1/3}$, with $V$ the volume of
the bounding box of the tracers and the given randoms (of the tracers alone
when the randoms are drawn):
{cpp:func}`otswap::reconstructBox`, {py:func}`otswap.reconstruct_box`. Randoms
may be given, or drawn uniformly in the bounding box of the tracers, which is
then the box. Use it for any catalogue in comoving Cartesian
coordinates whose density is about constant, cubic or not.

The box may be given instead, as closed limits on each axis
({cpp:struct}`otswap::BoxLimits`, {py:class}`otswap.BoxLimits`, or in Python a
pair of corners): the tracers outside it are left out, the randoms given
outside it are dropped, drawn randoms fill it, and the mps is computed from its
volume and the tracers inside. Use it when the volume of the sample is known,
such as a simulation box, rather than inferred from its outermost tracers. The
result records the box used: the limits, or the bounding box of the tracers and
the given randoms.

By default the box is not periodic: nothing flows through its faces, so modes
on the scale of the box itself cannot be reconstructed, and the displacements
are least reliable near the faces. A sub-volume cut from a periodic simulation
is such a box.

The whole volume of a periodic simulation is reconstructed as a periodic box:
`periodic` (`periodic=True`) with the limits, which are then required. Each
side must be an integer multiple of the shortest, as in a cube or a
500 × 500 × 1000 cuboid. Every distance is then the minimum image, and nothing
is left out: tracers and given randoms outside the limits are wrapped into
them by whole periods and counted (`tracersWrapped`, `randomsWrapped`), and
the result holds the wrapped tracers. Each displacement is the minimum-image
vector, within [−L/2, L/2) on an axis of period L, and the mean Lagrangian
position is wrapped into the box. The mps, when computed, is that of the
volume of the limits and every tracer. The redshift-space correction of a
periodic result finds its neighbours across the faces and wraps the corrected
positions.

**Lightcone.** Sky coordinates, converted with a
{cpp:class}`otswap::DistanceTable` ({py:class}`otswap.DistanceTable`), and an
mps that changes with redshift, measured from the tracers in `nBins`
(`n_bins`) redshift bins and kept in the result
({cpp:struct}`otswap::MpsProfile`, {py:class}`otswap.MpsProfile`):
{cpp:func}`otswap::reconstructLightcone`,
{py:func}`otswap.reconstruct_lightcone`. The survey's area comes from a
{cpp:class}`otswap::Mask` or is given. Randoms are required: they carry the
survey's geometry, selection function and completeness, which otswap does not
model. A redshift cut may restrict the reconstruction to a slice. The same
profile, bit for bit, is available without the randoms and without
reconstructing, to choose `nBins` beforehand:
{cpp:func}`otswap::measureMpsProfile`, {py:func}`otswap.measure_mps_profile`.

A lightcone may also be given in comoving Cartesian coordinates, the observer
at the origin: the Cartesian overloads of `reconstructLightcone` and
`measureMpsProfile`, {py:func}`otswap.reconstruct_lightcone_cartesian` and
{py:func}`otswap.measure_mps_profile_cartesian` in Python. They need no
cosmology. The mps is measured in `nBins` bins of comoving distance from the
origin, a profile in distance, and a distance cut
({cpp:struct}`otswap::DistanceCut`, {py:class}`otswap.DistanceCut`) takes the
place of the redshift cut. The mask is applied on the right ascension and
declination of each point, {cpp:func}`otswap::pointToSky`
({py:func}`otswap.point_to_sky`). A distance table may be given: it plays no
part in the reconstruction, and only fills the redshifts of the result, NaN
without it. The redshift-space correction needs one, for the growth rate; it
takes a tracer's redshift from the result, or from its distance when the
result has none.

## Coordinates and units

- Cartesian coordinates and every length (the mps, box limits,
  displacements, σ) are in coordinate units: the unit of the tracers given in a
  box or in a lightcone given in Cartesian coordinates, otswap assuming none;
  in a lightcone given on the sky that of the distance table, Mpc/h when it
  comes from the cosmology constructor, the unit of the table given when it
  comes from a table. A table given with Cartesian coordinates must be in
  their unit. The files otswap writes give lengths no unit for that reason.
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

## Realizations and validity

The reconstruction is repeated `nRealizations` (`n_realizations`) times, each
on its own set of N points extracted from the random catalogue, disjoint from
the others, so the catalogue must hold at least `nRealizations` × N points.
Each displacement has a validity flag, `valid`, per realization and object;
every flag of a tracer that took part is set after the reconstruction, and
the mask filter only ever clears them. For each tracer, `validRealizations`
(`valid_realizations`) counts the valid ones, and `meanDisplacement`
(`mean_displacement`) is their mean, NaN when there is none. `lagrangian`,
the tracer plus its mean displacement, follows. As the number of realizations
grows, `lagrangian` tends to the centroid of the tracer's equal-volume Laguerre
cell.

After changing `valid` by hand, {cpp:func}`otswap::recomputeMeans`
({py:func}`otswap.recompute_means`) recomputes the fields that follow from it,
as the reconstruction computes them.

Tracers left out of a lightcone reconstruction, by its redshift or distance
cut or its mask, keep their row: they are flagged in `outsideRedshiftCut`,
`outsideDistanceCut` and `outsideMask` (`outside_redshift_cut`,
`outside_distance_cut`, `outside_mask`), with NaN displacements and no valid
realization. Tracers outside the limits of a box that is not periodic are
flagged in the same way, in `outsideBox` (`outside_box`). Every output keeps
one row per input object, in input order.

A mean displacement file read back with {cpp:func}`otswap::io::readDisplacementMean`
({py:func}`otswap.io.read_displacement_mean`) gives a result of one
realization holding the mean: `nRealizationsMean` (`n_realizations_mean`)
records how many realizations were averaged, and `validRealizations` how many
of them were valid, so it may exceed 1.

To carry the catalogue's own IDs into the files, set them on the result
before writing or correcting it: `Result::ids` (`result.ids = ...`), one
64-bit integer per object, no value twice. The writers then add an `id`
column to their default columns, and the real-space catalogue copies them.
Without IDs, an `id` column, where written, holds each object's row, from 0.

## The mask

A HEALPix map, from a FITS file or from an array in memory. A pixel is
observed when its value is greater than 0; the value is not a weight. Given to
the lightcone reconstruction, the mask is applied twice:

- before it, tracers and randoms on unobserved pixels are left out;
- after it, a displacement whose great-circle arc crosses more than
  `maxUnobservedPixelsCrossed` (`max_unobserved_pixels_crossed`) unobserved
  pixels in a row is marked invalid, with 0, the default, at the first one
  ({cpp:func}`otswap::rejectMaskCrossings`,
  {py:func}`otswap.reject_mask_crossings`).

The count is of pixels, not an angle, so the same threshold means a different
criterion at a different NSIDE. A second filter with a higher threshold
changes nothing: set the threshold in the call that reconstructs.

## Reproducibility and threads

```{include} ../README.md
:start-after: <!-- docs:reproducibility-start -->
:end-before: <!-- docs:reproducibility-end -->
```

The seed used is `Result::config.seed` in C++, `Result.config.seed` in Python. The
number of threads follows `OMP_NUM_THREADS`. The reconstructions, the mask
filter, the redshift-space correction and the table readers and writers of
`otswap.io` release the GIL while they run. A `Result` must not be read from
one thread while another modifies it in place, as
{cpp:func}`otswap::rejectMaskCrossings` and {cpp:func}`otswap::recomputeMeans`
({py:func}`otswap.reject_mask_crossings`, {py:func}`otswap.recompute_means`)
do.

## The parameters that matter

- **`convergence`**, default 1e-3: the sweeps stop once the fraction of
  successful swaps in a sweep, swaps per tracer visited, no longer exceeds
  this threshold. Larger values are faster and less accurate.
- **`nRealizations`** (`n_realizations`), default 30: more realizations
  bring the mean Lagrangian position closer to the centroid of the tracer's
  Laguerre cell, and need N more randoms each.
- **`cellSize`** (`cell_size`), default 3 mps: affects the speed only, never
  the result.
- **σ of the redshift-space correction:** for coordinates in Mpc/h, 10 Mpc/h
  is a reasonable starting
  value; the best value depends on the sample, and should be checked in each
  analysis.

The [method](method.md) page gives the reference for the choice of the
parameters.

## What it prints

Each call reports, as `verbosity` says ({cpp:enum}`otswap::Verbosity`,
{py:data}`otswap.Verbosity`), in lines starting with `otswap:`:

- **silent:** nothing;
- **normal:** the selection report of a lightcone, when a selection applied,
  or of a box, when limits were given, and in C++ the extrapolation of b(z),
  when it happened; then one line per call with its time and what was lost;
- **detailed:** as normal, with the lines that explain it: the mps of a box
  and its source, the mps profile of a lightcone, the breakdown of the
  corrected and of the uncorrected tracers.

C++ writes to `Config::log`, `std::clog` by default; Python to `sys.stdout`,
which a notebook shows. Each line is written as soon as its step finishes, so
a long run shows its progress, and a call that fails keeps the lines of the
steps it finished. In Python, Ctrl+C stops a reconstruction,
`measure_mps_profile` or a correction at its next line with
`KeyboardInterrupt`. The counts behind every line are in the object returned.

## Errors

Invalid input and failures raise {cpp:class}`otswap::Error`
({py:class}`otswap.Error`), whose message names the object at fault. An error
of the neighbour search is one too: its message starts with `meshsearch: `,
and in C++ the original `meshsearch::Error` is nested in it
(`std::rethrow_if_nested` recovers it). The extrapolation of b(z) beyond its
table is not an error: Python issues an
{py:class}`otswap.ExtrapolationWarning`, C++ reports it and counts it in
`nExtrapolated`.
