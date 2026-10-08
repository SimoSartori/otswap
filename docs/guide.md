# Guide

What applies to both APIs. The names are the C++ ones and the Python ones
side by side; [names side by side](naming.md) maps them all.

## Two geometries

**Box.** Cartesian coordinates, and a single mean particle separation (mps),
given or computed from the tracers as $(V/N)^{1/3}$, with $V$ the volume of
their bounding box: {cpp:func}`otswap::reconstructBox`,
{py:func}`otswap.reconstruct_box`. Randoms may be given, or drawn uniformly in
the tracers' bounding box. Use it for any catalogue in comoving Cartesian
coordinates whose density is about constant, cubic or not.

The box is not periodic: nothing flows through its faces, so modes on the
scale of the box itself cannot be reconstructed, and the displacements are
least reliable near the faces. A box cut from a periodic simulation is treated
as a plain sub-volume.

**Lightcone.** Sky coordinates, converted with a
{cpp:class}`otswap::DistanceTable` ({py:class}`otswap.DistanceTable`), and an
mps that changes with redshift, measured from the tracers in `nBins`
(`n_bins`) redshift bins and kept in the result
({cpp:struct}`otswap::MpsProfile`, {py:class}`otswap.MpsProfile`):
{cpp:func}`otswap::reconstructLightcone`,
{py:func}`otswap.reconstruct_lightcone`. The survey's area comes from a
{cpp:class}`otswap::Mask` or is given. Randoms are required: they carry the
survey's geometry, selection function and completeness, which otswap does not
model. A redshift cut may restrict the reconstruction to a slice.

## Coordinates and units

- Distances are in Mpc/h, the distance table's unit.
- Sky coordinates are ordered right ascension, declination, redshift.
- In C++ angles are in radians; {cpp:func}`otswap::skyToRadians` and
  {cpp:func}`otswap::skyToDegrees` convert a sky array. In Python every
  function that takes or returns angles has a required `angle_unit`, `"deg"`
  or `"rad"`.
- Right ascensions come back folded into [0, 2π), or [0, 360) degrees.
- In Python a coordinate array has shape (N, 3) and is converted to float64;
  the arrays of a result have the shapes below, with the realization first.

The layout of the C++ arrays, which the Python arrays map onto:

```{doxygengroup} coordinate_layout
```

## Realizations and validity

The reconstruction is repeated `nRealizations` (`n_realizations`) times, each
on its own set of N points extracted from the random catalogue, disjoint from
the others, so the catalogue must hold at least `nRealizations` × N points. Each displacement has a validity flag, `valid`, per
realization and object; every flag is set after the reconstruction, and the
mask filter only ever clears them. For each tracer, `validRealizations`
(`valid_realizations`) counts the valid ones, and `meanDisplacement`
(`mean_displacement`) is their mean, NaN when there is none. `lagrangian`,
the tracer plus its mean displacement, follows. As the number of realizations
grows, `lagrangian` tends to the centroid of the tracer's equal-volume Laguerre
cell.

After changing `valid` by hand, {cpp:func}`otswap::recomputeMeans`
({py:func}`otswap.recompute_means`) recomputes the fields that follow from it,
as the reconstruction computes them.

Tracers left out of a lightcone reconstruction, by its redshift cut or its
mask, keep their row: they are flagged in `outsideRedshiftCut` and
`outsideMask` (`outside_redshift_cut`, `outside_mask`), with NaN
displacements and no valid realization. Every output keeps one row per input
object, in input order.

## The mask

A HEALPix map, from a FITS file or from an array in memory. A pixel is
observed when its value is greater than 0; the value is not a weight. Given to
the lightcone reconstruction, the mask is applied twice:

- before it, tracers and randoms on unobserved pixels are left out;
- after it, a displacement whose great-circle arc crosses more than
  `maxUnobservedPixelsCrossed` (`max_unobserved_pixels_crossed`) distinct
  unobserved pixels is marked invalid, with 0, the default, at the first one
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

The seed used is `Result::config.seed` in C++, `Result.seed` in Python. The
number of threads follows `OMP_NUM_THREADS`. The reconstructions, the mask
filter, the redshift-space correction and the table readers and writers of
`otswap.io` release the GIL while they run.

## The parameters that matter

- **`convergence`**, default 1e-3: the sweeps stop once the fraction of
  successful swaps in a sweep, swaps per tracer visited, no longer exceeds
  this threshold. Larger values are faster and less accurate.
- **`nRealizations`** (`n_realizations`): more realizations bring the mean
  Lagrangian position closer to the centroid of the tracer's Laguerre cell,
  and need N more randoms each.
- **`cellSize`** (`cell_size`), default 4 mps: affects the speed only, never
  the result.
- **σ of the redshift-space correction:** 10 Mpc/h is a reasonable starting
  value; the best value depends on the sample, and should be checked in each
  analysis.

The [method](method.md) page gives the reference for the choice of the
parameters.

## What it prints

Each call reports, as `verbosity` says ({cpp:enum}`otswap::Verbosity`,
{py:data}`otswap.Verbosity`), in lines starting with `otswap:`:

- **silent:** nothing;
- **normal:** the selection report of a lightcone, when a selection applied,
  and in C++ the extrapolation of b(z), when it happened; then one line per
  call with its time and what was lost;
- **detailed:** as normal, with the lines that explain it: the mps of a box
  and its source, the mps(z) profile of a lightcone, the breakdown of the
  corrected and of the uncorrected tracers.

C++ writes to `Config::log`, `std::clog` by default; Python to `sys.stdout`,
which a notebook shows. The counts behind every line are in the object
returned.

## Errors

Invalid input and failures raise {cpp:class}`otswap::Error`
({py:class}`otswap.Error`), whose message names the object at fault. Two
exceptions:

- the neighbour search of the redshift-space correction raises
  `meshsearch::Error` (in Python a `RuntimeError`) for positions that span no
  volume;
- the extrapolation of b(z) beyond its table is not an error: Python issues an
  {py:class}`otswap.ExtrapolationWarning`, C++ reports it and counts it in
  `nExtrapolated`.
