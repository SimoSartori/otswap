# Guide

How to use otswap, in both APIs: the two geometries, coordinates and units,
realizations and validity, the mask, the redshift-space correction,
reproducibility, and the parameters to choose. The names are the C++ ones and
the Python ones side by side; [names side by side](naming.md) maps them all.

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

By default the box is not periodic: its faces are boundaries, with the
consequences that {cpp:func}`otswap::reconstructBox` describes. A sub-volume
cut from a periodic simulation is such a box.

The whole volume of a periodic simulation is reconstructed as a periodic box:
`periodic` (`periodic=True`) with the limits, which are then required, each
side an integer multiple of the shortest, as in a cube or a 500 × 500 × 1000
cuboid ({cpp:func}`otswap::reconstructBox` gives the rule). Every distance is then the minimum image, and nothing
is left out: tracers and given randoms outside the limits are wrapped into
them by whole periods and counted (`tracersWrapped`, `randomsWrapped`), and
the result holds the wrapped tracers. Each displacement is the minimum-image
vector, within [−L/2, L/2) on an axis of period L within rounding, and the mean Lagrangian
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
{cpp:class}`otswap::Mask` or is given in square degrees (`skyAreaDeg2`,
`sky_area_deg2`). Randoms are required: they carry the
survey's geometry, selection function and completeness, which otswap does not
model. A redshift cut may restrict the reconstruction to a slice. The same
profile, bit for bit, is available without the randoms and without
reconstructing, to choose `nBins` beforehand (see [Choosing the
parameters](#choosing-the-parameters)):
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

## Realizations and validity

The reconstruction is repeated `nRealizations` (`n_realizations`) times, each
on its own set of N points extracted from the random catalogue, disjoint from
the others, so the catalogue must hold at least `nRealizations` × N points.
Each displacement has a validity flag, `valid`, per realization and object;
every flag of a tracer that took part is set after the reconstruction, and
the mask filter only ever clears them. For each tracer, `validRealizations`
(`valid_realizations`) counts the valid ones, and `meanDisplacement`
(`mean_displacement`) is their mean, NaN when there is none. `lagrangian`,
the tracer plus its mean displacement, follows; it approximates the centroid
of the tracer's Laguerre cell, as the [method](method.md) page explains.

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

The mask describes the survey's footprint on the sky, so that a lightcone
reconstruction uses only what was observed and discards the displacements that
pass through what was not. It is a HEALPix map, from a FITS file or from an
array in memory. A pixel is observed when its value is greater than 0; the
value is not a weight. Given to the lightcone reconstruction, the mask is
applied twice:

- before it, tracers and randoms on unobserved pixels are left out;
- after it, a displacement whose great-circle arc crosses more than
  `maxUnobservedPixelsCrossed` (`max_unobserved_pixels_crossed`) unobserved
  pixels in a row is marked invalid, with 0, the default, at the first one;
  {cpp:func}`otswap::rejectMaskCrossings`
  ({py:func}`otswap.reject_mask_crossings`) gives the rule.

The count is of pixels, not an angle, so the same threshold means a different
criterion at a different NSIDE. A second filter with a higher threshold
changes nothing: set the threshold in the call that reconstructs, or set
`rejectCrossings` (`reject_crossings`) to false and filter afterwards with
{cpp:func}`otswap::rejectMaskCrossings`.

## The redshift-space correction

The displacements of a reconstruction in redshift space move the catalogue to
real space: {cpp:func}`otswap::realSpaceBox` and
{cpp:func}`otswap::realSpaceLightcone` ({py:func}`otswap.real_space_box`,
{py:func}`otswap.real_space_lightcone`) shift each tracer along the line of
sight by the line-of-sight component of its displacement, smoothed over its
neighbours with a Gaussian of width σ and scaled by $f/(b + 3f/5)$; the
[method](method.md) page derives the factor. A box takes one bias and one
redshift; a lightcone takes b(z) as a table
({cpp:struct}`otswap::BiasTable`, {py:class}`otswap.BiasTable`), extrapolated
linearly beyond its nodes with a warning, and the growth rate f from the
distance table. The result gives each tracer a status
({cpp:enum}`otswap::CorrectionStatus`, {py:class}`otswap.CorrectionStatus`):
moved with its own displacement, moved with its neighbours' average when it
has no valid realization, not moved when no valid tracer lies within 3σ, or
left out of the reconstruction. `weightByRealizations`
(`weight_by_realizations`) also weights each neighbour by its count of valid
realizations.

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

## Choosing the parameters

The defaults are the starting point for an analysis; the examples use smaller
values to run quickly.

- **`convergence`**, default 1e-3: the sweeps stop once the fraction of
  successful swaps in a sweep, swaps per tracer visited, no longer exceeds
  this threshold. In the void analysis of Sartori et al. (2026, Appendix B),
  void abundances were stable for thresholds of about 1e-2 or below, and
  stacked void profiles for 1e-3 or below; larger values are faster.
- **`nRealizations`** (`n_realizations`), default 30: that analysis used 10 to
  50, depending on the geometry, the small scales needing the most. The random
  catalogue must hold `nRealizations` × N points.
- **`nBins`** (`n_bins`), lightcone only: each bin must hold at least 278
  tracers, for a 2% accuracy on its mps, and the error names the largest
  number of bins the catalogue supports. Take enough bins to follow the
  variation of the density with redshift or distance, and choose them with
  {cpp:func}`otswap::measureMpsProfile` ({py:func}`otswap.measure_mps_profile`)
  before reconstructing.
- **`maxUnobservedPixelsCrossed`** (`max_unobserved_pixels_crossed`),
  default 0: 0 rejects any crossing; to tolerate small holes such as
  bright-star masks, set the size of those holes in pixels at the mask's
  NSIDE.
- **`mps`**, box only: computed by default as $(V/N)^{1/3}$ of the box; give
  it when the catalogue's mean density is known better than its volume and
  count give it.
- **`periodic`**, box only: for the whole volume of a simulation with periodic
  boundary conditions; a sub-volume cut from one is not periodic.
- **Redshift or distance cut**
  ({cpp:struct}`otswap::RedshiftCut`, {cpp:struct}`otswap::DistanceCut`):
  to reconstruct a slice of a wider catalogue; the randoms outside it are
  dropped with the tracers.
- **`cellSize`** (`cell_size`), default 3 mps: affects the speed only, never
  the result. In a periodic box the cells divide the shortest side a whole
  number of times, so the side used can be smaller: `Result::cellSide`
  (`Result.cell_side`) reports it.
- **σ of the redshift-space correction:** for coordinates in Mpc/h, 10 Mpc/h
  is a reasonable starting value; the best value depends on the sample, and
  should be checked in each analysis.
- **`weightByRealizations`** (`weight_by_realizations`), default false:
  set it to give less weight to tracers whose displacement averages few valid
  realizations, as near the edges of a masked footprint.

## What it prints

Each call reports what it did, in lines starting with `otswap:`, at the level
`verbosity` sets: silent, normal (the default: what was selected and lost)
or detailed ({cpp:enum}`otswap::Verbosity`, {py:data}`otswap.Verbosity`, which
list the lines). C++ writes to the stream set in `Config::log`,
`CorrectionConfig::log` or the `log` argument of `measureMpsProfile`,
`std::clog` by default; Python to `sys.stdout`, which a notebook shows. In
Python, Ctrl+C stops a reconstruction, `measure_mps_profile` or a correction
at its next line with `KeyboardInterrupt`.

## Errors

Invalid input and failures raise {cpp:class}`otswap::Error`
({py:class}`otswap.Error`); [names side by side](naming.md#errors) describes
them in both languages.
