"""Optimal transport reconstruction of the displacement field of a tracer
catalogue, by local swapping.

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates are in
Mpc/h; sky coordinates are ordered right ascension, declination, redshift,
with the angular unit given explicitly by ``angle_unit``. Inputs are
converted to float64; the reconstruction works on its own copy of them.

The number of threads follows ``OMP_NUM_THREADS``. The reconstruction
functions release the GIL while they run.

Tables are read and written by the submodule ``otswap.io``.

The redshift-space correction of tracer i is a shift along its line of
sight::

    s_i = f(z_i) / (b(z_i) + 3 f(z_i) / 5) * <Psi . r_hat>_i

with Psi the reconstructed displacement (``Result.mean_displacement``, from
the observed to the reconstructed position), r_hat the line of sight, f the
linear growth rate, b the linear bias and < > a gaussian average over the
neighbouring tracers. ``real_space_lightcone`` and ``real_space_box`` run the
whole chain; the four steps are also available on their own.
"""

from enum import IntEnum
from typing import Literal, NamedTuple, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray
from typing_extensions import TypeAlias

from . import io as io

AngleUnit: TypeAlias = Literal["deg", "rad"]

Verbosity: TypeAlias = Literal["silent", "normal", "detailed"]
"""What a reconstruction or a correction prints, to ``sys.stdout``, each line
starting with ``otswap:``:

- ``"silent"``: nothing;
- ``"normal"``: the selection report of a lightcone, when a selection
  applied, then one line per call with its time and what was lost::

      otswap: reconstructLightcone: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement
      otswap: realSpaceLightcone: corrected <C> of <N> tracers in <T> s; <U> left uncorrected

- ``"detailed"``: as ``"normal"``, with the lines that explain it: the mean
  particle separation of a box and its source, the mps(z) profile of a
  lightcone, and the breakdown of the corrected and of the uncorrected.

The counts behind every line are in the object returned. The extrapolation
of b(z) is an ``ExtrapolationWarning``, whatever the verbosity.
"""


class Error(RuntimeError):
    """Raised by every otswap function on invalid input or failure."""


class ExtrapolationWarning(UserWarning):
    """Issued, once per call, when b(z) is extrapolated beyond the redshift
    range of its table."""


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------

class Result:
    """Displacement field produced by a reconstruction.

    Displacements point from the observed (Eulerian) position of a tracer to
    its reconstructed (Lagrangian) position, which is the position of the
    random it is matched to.

    Every array attribute is read-only. The sky arrays (``tracers_sky``,
    ``lagrangian_sky``) are in the ``angle_unit`` of the call that made the
    result; in degrees they are converted once, on first access, and kept.
    Every other array is a view on memory owned by this object, and keeps it
    alive. Use ``.copy()`` to obtain a writable array. Filtering with
    ``reject_mask_crossings`` updates ``valid``, ``valid_realizations``,
    ``mean_displacement`` and ``lagrangian`` in place, so views taken earlier
    reflect it; ``lagrangian_sky`` is recomputed, and an array taken from it
    in degrees before the filter is not updated.
    """

    @staticmethod
    def from_arrays(
        displacement: ArrayLike,
        matched_random: ArrayLike,
        valid: ArrayLike,
        tracers: ArrayLike,
        *,
        tracers_sky: Optional[ArrayLike] = None,
        angle_unit: Optional[AngleUnit] = None,
        distances: Optional["DistanceTable"] = None,
        outside_redshift_cut: Optional[ArrayLike] = None,
        outside_mask: Optional[ArrayLike] = None,
        seed: int = 0,
    ) -> "Result":
        """Build a result from its arrays, as read back from a file.

        Parameters
        ----------
        displacement, matched_random : array_like, shape (n_realizations, n_objects, 3)
        valid : array_like, shape (n_realizations, n_objects)
            Each entry 0 or 1 (or bool).
        tracers : array_like, shape (n_objects, 3)
            Cartesian.
        tracers_sky : array_like, shape (n_objects, 3), optional
            For a lightcone result, with its ``angle_unit`` and the
            ``distances`` table of the reconstruction, from which
            ``lagrangian_sky`` is computed. Without it the result is a box
            result, and neither may be given.
        outside_redshift_cut, outside_mask : array_like, shape (n_objects,), optional
            0 or 1; all 0 when omitted.
        seed : int, optional
            The seed to record.

        Notes
        -----
        ``valid_realizations``, ``mean_displacement``, ``lagrangian`` and
        ``lagrangian_sky`` are computed as ``recompute_means`` computes them,
        so arrays read back exactly give the reconstruction's values, bit for
        bit. The profile, the mps and the selection counts are not
        reconstructed: ``mps`` and ``mps_profile`` are None.
        """

    @property
    def n_objects(self) -> int: ...

    @property
    def n_realizations(self) -> int: ...

    @property
    def geometry(self) -> Literal["box", "lightcone"]: ...

    @property
    def angle_unit(self) -> Optional[AngleUnit]:
        """The unit of ``tracers_sky`` and ``lagrangian_sky``: the
        ``angle_unit`` of the call that made the result. None for a box."""

    @property
    def seed(self) -> int:
        """The seed used: the one drawn when ``seed=0`` was given, so that the
        result can be reproduced."""

    @property
    def tracers(self) -> NDArray[np.float64]:
        """Cartesian position of each tracer, the start of its
        displacements, shape (n_objects, 3): the array given, or the
        conversion of ``tracers_sky``. A tracer left out by the redshift cut
        or the mask has its conversion when its redshift lies in the
        distance table, NaN otherwise."""

    @property
    def tracers_sky(self) -> Optional[NDArray[np.float64]]:
        """Sky coordinates of each tracer, shape (n_objects, 3), in
        ``angle_unit``; None for a box. In degrees, the angles come back as
        ``x * (pi/180) * (180/pi)``, which can differ from the input in the
        last bit."""

    @property
    def lagrangian(self) -> NDArray[np.float64]:
        """Mean Lagrangian position of each tracer, ``tracers +
        mean_displacement``, shape (n_objects, 3). NaN rows where the mean
        displacement is NaN."""

    @property
    def mps(self) -> Optional[float]:
        """Box: the mean particle separation used, given or computed, in
        Mpc/h. None for a lightcone, whose separation is ``mps_profile``."""

    @property
    def mps_profile(self) -> Optional["MpsProfile"]:
        """Lightcone: the mean particle separation as a function of redshift,
        measured on the tracers kept. None for a box. Keeps this result
        alive."""

    @property
    def elapsed_seconds(self) -> float:
        """Wall time of the call that made the result, in seconds. Two runs
        with the same seed give the same result but for this value."""

    @property
    def displacement(self) -> NDArray[np.float64]:
        """Displacements, shape (n_realizations, n_objects, 3)."""

    @property
    def matched_random(self) -> NDArray[np.float64]:
        """Cartesian position of the random matched to each tracer, shape
        (n_realizations, n_objects, 3)."""

    @property
    def valid(self) -> NDArray[np.bool_]:
        """Whether each displacement is valid, shape
        (n_realizations, n_objects). All true unless a filter was applied."""

    @property
    def valid_realizations(self) -> NDArray[np.uint32]:
        """Number of valid realizations per tracer, shape (n_objects,)."""

    @property
    def mean_displacement(self) -> NDArray[np.float64]:
        """Displacement averaged over the valid realizations, shape
        (n_objects, 3). NaN for a tracer with no valid realization."""

    @property
    def filtered_nside(self) -> int:
        """NSIDE of the mask the result was filtered against, 0 if none."""

    @property
    def outside_redshift_cut(self) -> NDArray[np.bool_]:
        """Whether each tracer lay outside the ``redshift_cut`` of
        ``reconstruct_lightcone`` and took no part in the reconstruction,
        shape (n_objects,). All false without a cut. A cut tracer has NaN
        displacements and matched randoms, no valid realization and a NaN
        mean displacement."""

    @property
    def outside_mask(self) -> NDArray[np.bool_]:
        """Whether each tracer fell on an unobserved pixel of the ``mask`` of
        ``reconstruct_lightcone`` and took no part in the reconstruction,
        shape (n_objects,). All false without a mask. Independent of
        ``outside_redshift_cut``: a tracer may carry both flags. A flagged
        tracer's rows are as for a cut one."""

    @property
    def selection(self) -> "SelectionCounts":
        """What ``reconstruct_lightcone`` left out, and why. Keeps this
        result alive."""

    @property
    def lagrangian_sky(self) -> Optional[NDArray[np.float64]]:
        """Sky coordinates of each tracer's mean Lagrangian position,
        ``lagrangian``, shape (n_objects, 3), in ``angle_unit``: right
        ascension in [0, 360) degrees ([0, 2 pi) radians), declination,
        redshift, as ``to_sky`` computes them with the reconstruction's
        distance table. NaN rows where ``lagrangian`` is NaN: the tracers left
        out by a selection and those with no valid realization; a position
        whose distance falls outside the table has a NaN redshift. None for a
        box. Recomputed by ``reject_mask_crossings``."""


class MpsProfile:
    """The mean particle separation of a lightcone as a function of
    redshift, measured by ``reconstruct_lightcone``: the tracers kept in
    ``n_bins`` uniform redshift bins over [``redshift_min``,
    ``redshift_max``], mps = (N / V)^(-1/3) in each, joined linearly and
    extrapolated linearly beyond the outermost bins. Read-only."""

    @property
    def redshift(self) -> NDArray[np.float64]:
        """Bin centres, shape (n_bins,)."""

    @property
    def mps(self) -> NDArray[np.float64]:
        """Mean particle separation at each centre, Mpc/h, shape (n_bins,)."""

    @property
    def count(self) -> NDArray[np.uint64]:
        """Tracers in each bin, shape (n_bins,)."""

    @property
    def redshift_min(self) -> float: ...

    @property
    def redshift_max(self) -> float: ...

    @property
    def representative(self) -> float:
        """The count-weighted median of ``mps``, which sizes the grids of the
        reconstruction."""

    def at(self, z: ArrayLike) -> NDArray[np.float64]:
        """The mps at each redshift, Mpc/h, the shape of ``z`` kept: the value
        the reconstruction used at a tracer of that redshift. Raises where the
        extrapolation is not positive."""


class SelectionCounts:
    """What a lightcone reconstruction left out, and why. Read-only.

    The tracers removed are ``tracers_outside_redshift_cut +
    tracers_outside_mask - tracers_outside_both``, and the same for the
    randoms. All counts are 0, and nothing is applied, for ``reconstruct_box``.
    ``repr`` gives the report ``reconstruct_lightcone`` prints.
    """

    @property
    def redshift_cut(self) -> Optional[tuple[float, float]]:
        """The (min, max) of the redshift cut; None when no cut was applied."""

    @property
    def mask_applied(self) -> bool:
        """Whether a mask was given."""

    @property
    def tracers(self) -> int:
        """Tracers given."""

    @property
    def tracers_outside_redshift_cut(self) -> int: ...

    @property
    def tracers_outside_mask(self) -> int: ...

    @property
    def tracers_outside_both(self) -> int:
        """Tracers counted in both of the two above."""

    @property
    def randoms(self) -> int:
        """Randoms given."""

    @property
    def randoms_outside_redshift_cut(self) -> int: ...

    @property
    def randoms_outside_mask(self) -> int: ...

    @property
    def randoms_outside_both(self) -> int:
        """Randoms counted in both of the two above."""

    @property
    def max_unobserved_pixels_crossed(self) -> Optional[int]:
        """The threshold of the latest mask filter applied to the result, by
        ``reconstruct_lightcone`` (``reject_crossings``) or by
        ``reject_mask_crossings`` afterwards; None when none was."""

    @property
    def displacements(self) -> int:
        """Valid displacements before the first mask filter."""

    @property
    def displacements_crossing_mask(self) -> int:
        """Of them, those rejected by every mask filter so far."""


# ---------------------------------------------------------------------------
# Cosmology
# ---------------------------------------------------------------------------

class DistanceTable:
    """Sampled relation between redshift, comoving distance and, optionally,
    the linear growth rate, interpolated in both directions. Immutable.

    Build it with ``DistanceTable.flat`` or ``DistanceTable.from_table``.
    """

    @staticmethod
    def flat(
        omega_m: float,
        h: float,
        w0: float = -1.0,
        wa: float = 0.0,
        z_min: float = 0.0,
        z_max: float = 10.0,
        n_samples: int = 50_000,
    ) -> "DistanceTable":
        """Tabulate a flat cosmology with dark energy w(a) = w0 + wa (1 - a).

        Comoving distances are in Mpc/h. The growth rate f = dlnD/dlna comes
        from integrating the linear growth equation. The range must cover
        every redshift to be converted, including reconstructed positions,
        which can fall slightly outside the range of the data.
        """

    @staticmethod
    def from_table(
        redshift: ArrayLike,
        distance: ArrayLike,
        growth_rate: Optional[ArrayLike] = None,
    ) -> "DistanceTable":
        """Use a table computed by the caller. Redshift and distance must be
        strictly increasing, with at least two entries; distance in Mpc/h.
        Without ``growth_rate``, ``growth_rate_at`` raises."""

    @property
    def min_redshift(self) -> float: ...

    @property
    def max_redshift(self) -> float: ...

    @property
    def has_growth_rate(self) -> bool: ...

    def distance_at(self, z: ArrayLike) -> NDArray[np.float64]:
        """Comoving distance at each redshift. Raises outside the table."""

    def redshift_at(self, distance: ArrayLike) -> NDArray[np.float64]:
        """Redshift at each comoving distance. Raises outside the table."""

    def growth_rate_at(self, z: ArrayLike) -> NDArray[np.float64]:
        """Linear growth rate at each redshift. Raises outside the table, or
        if the table has no growth rate."""


def to_cartesian(
    sky: ArrayLike,
    distances: DistanceTable,
    *,
    angle_unit: AngleUnit,
) -> NDArray[np.float64]:
    """Convert sky coordinates, shape (N, 3), to Cartesian ones in Mpc/h,
    using the same convention as ``reconstruct_lightcone``. Raises if a
    declination lies outside [-90, 90] degrees ([-pi/2, pi/2] radians); the
    message gives it in ``angle_unit``."""


def to_sky(
    cartesian: ArrayLike,
    distances: DistanceTable,
    *,
    angle_unit: AngleUnit,
) -> NDArray[np.float64]:
    """Convert Cartesian coordinates, shape (N, 3), in Mpc/h, to sky
    coordinates, shape (N, 3): right ascension in [0, 360) degrees
    ([0, 2 pi) radians), declination, redshift. The inverse of
    ``to_cartesian``, with the same observer; a position at the origin gets
    right ascension and declination 0. A row holding a NaN gives a NaN
    row; an infinite entry, or a distance outside the table, raises."""


# ---------------------------------------------------------------------------
# Angular mask
# ---------------------------------------------------------------------------

class Mask:
    """HEALPix mask, read from a FITS file or built from a full-sky map in
    memory with ``Mask.from_array``. Immutable.

    A pixel is observed when its value is greater than 0, whatever its
    magnitude: fractional values, the smallest positive ones and +inf
    included. 0, negative values, -inf, NaN and UNSEEN mark an unobserved
    pixel. Each value is tested once, in float64, when the mask is built; no
    value is kept or used as a weight. The rule is the same for both
    sources.

    NSIDE may be any integer from 1 to 2**29, HEALPix's limit, and must be a
    power of 2 in the NESTED scheme. The mask holds one byte per pixel,
    12 NSIDE**2 bytes: 12.6 MB at NSIDE 1024, 805 MB at NSIDE 8192, 3.2 GB
    at NSIDE 16384."""

    def __init__(self, fits_file: str) -> None:
        """NSIDE and ORDERING are read from the file header; RING and NESTED
        are both supported, and any numeric column type. The values are read
        a block of rows at a time, so reading needs memory for the mask's
        bytes only. A null value of an integer column is unobserved."""

    @staticmethod
    def from_array(values: ArrayLike, nest: bool = False) -> "Mask":
        """A full-sky HEALPix map given in memory, as healpy holds one.

        Parameters
        ----------
        values : array_like, shape (12 NSIDE**2,)
            One value per pixel, of any real dtype: integer, unsigned,
            floating, or bool as 0 and 1. It is read once and not kept. A
            C-contiguous float64 array is read in place; any other input is
            first converted to one, which takes 8 bytes per pixel while the
            mask is built.
        nest : bool, optional
            True for NESTED, False for RING, as healpy's ``nest``.

        Raises
        ------
        Error
            For another length, NSIDE above 2**29, NESTED with an NSIDE that
            is not a power of 2, an array that is not 1-D, values that are
            not real numbers, or a ``nest`` that is not a bool.
        """

    @property
    def nside(self) -> int:
        """NSIDE of the map, from 1 to 2**29."""

    @property
    def sky_area_deg2(self) -> float:
        """Area covered by the observed pixels, in square degrees: their
        count times the pixel area, whatever their values."""

    def allows(
        self,
        ra: ArrayLike,
        dec: ArrayLike,
        *,
        angle_unit: AngleUnit,
    ) -> NDArray[np.bool_]:
        """Whether each direction falls in an observed pixel. Raises if a
        right ascension is not finite, or a declination is not in
        [-90, 90] degrees ([-pi/2, pi/2] radians); the message gives it in
        ``angle_unit``."""


# ---------------------------------------------------------------------------
# Reconstruction
# ---------------------------------------------------------------------------

def reconstruct_box(
    tracers: ArrayLike,
    randoms: Optional[ArrayLike] = None,
    *,
    mps: Optional[float] = None,
    n_realizations: int = 1,
    convergence: float = 1e-3,
    seed: int = 0,
    cell_size: float = 4.0,
    verbosity: Verbosity = "normal",
) -> Result:
    """Reconstruct in box geometry, with a constant mean particle separation.

    Parameters
    ----------
    tracers : array_like, shape (N, 3)
        Cartesian coordinates.
    randoms : array_like, shape (M, 3), optional
        Cartesian coordinates, M >= n_realizations * N; each realization uses
        a disjoint subset. If omitted, n_realizations * N randoms are drawn
        uniformly in the bounding box of the tracers.
    mps : float, optional
        Mean particle separation, in the units of the coordinates. If
        omitted, (V / N)^(1/3), with V the volume of the tracers' bounding
        box, the box drawn randoms fill; it is ``Result.mps``. Raises if the
        tracers span no volume.
    n_realizations : int, optional
        Independent reconstructions to run.
    convergence : float, optional
        The sweeps stop once the fraction of successful swaps in a sweep,
        swaps per tracer visited, no longer exceeds this threshold.
    seed : int, optional
        Seed of the random streams; 0 draws one at random. With a fixed seed
        the result does not depend on the number of threads.
    cell_size : float, optional
        Grid cell size in units of the mps. Affects speed only.
    verbosity : Verbosity, optional
        What the call prints; see ``Verbosity``. ``"normal"`` prints one
        line, ``"detailed"`` the mps and its source before it.

    Notes
    -----
    The box is not periodic: nothing flows through its faces, so modes on the
    scale of the box itself cannot be reconstructed.
    """


def reconstruct_lightcone(
    tracers_sky: ArrayLike,
    randoms_sky: ArrayLike,
    *,
    sky_area_deg2: Optional[float] = None,
    mask: Optional[Mask] = None,
    n_bins: int,
    distances: DistanceTable,
    angle_unit: AngleUnit,
    tracers: Optional[ArrayLike] = None,
    randoms: Optional[ArrayLike] = None,
    n_realizations: int = 1,
    convergence: float = 1e-3,
    seed: int = 0,
    cell_size: float = 4.0,
    redshift_cut: Optional[tuple[float, float]] = None,
    reject_crossings: bool = True,
    max_unobserved_pixels_crossed: int = 0,
    verbosity: Verbosity = "normal",
) -> Result:
    """Reconstruct in lightcone geometry.

    The mean particle separation is measured from the tracers: n_bins
    uniform redshift bins over the observed range, mps = (N / V)^(-1/3) in
    each, with V the shell volume implied by the sky area. The bin values
    are joined by linear interpolation and extrapolated linearly beyond the
    outermost bins. Raises if a bin holds too few tracers for a 2% accuracy
    on its mps.

    Randoms are required: they carry the survey geometry, the selection
    function and the completeness.

    Parameters
    ----------
    tracers_sky : array_like, shape (N, 3)
        Sky coordinates of the tracers.
    randoms_sky : array_like, shape (M, 3)
        Sky coordinates of the randoms, M >= n_realizations * N.
    sky_area_deg2 : float, optional
        Effective survey area, when no mask is given: the published value.
        Give exactly one of ``sky_area_deg2`` and ``mask``.
    mask : Mask, optional
        The survey's mask. The sky area is then ``mask.sky_area_deg2``, and
        the mask is applied twice. Before the reconstruction, tracers and
        randoms on unobserved pixels are left out, as ``Mask.allows`` decides:
        the randoms are dropped, and the tracers keep their row, flagged in
        ``Result.outside_mask``. After it, with ``reject_crossings``,
        ``reject_mask_crossings`` is applied with
        ``max_unobserved_pixels_crossed``. The mask and the redshift cut are
        both evaluated on every object; the objects kept pass both.
    n_bins : int
        Redshift bins used to measure mps(z).
    distances : DistanceTable
        Table used for the conversion to Cartesian coordinates.
    angle_unit : AngleUnit
        Unit of right ascension and declination. A declination outside
        [-90, 90] degrees ([-pi/2, pi/2] radians) raises, the message giving
        it in this unit.
    tracers, randoms : array_like, optional
        Cartesian coordinates already computed, shape (N, 3) and (M, 3). Give
        both or neither; when given, the conversion is skipped, and their
        agreement with the sky coordinates is not checked.
    redshift_cut : tuple of two floats, optional
        (min, max). Tracers and randoms whose redshift lies outside this
        closed range are left out of the reconstruction, before anything
        else: the randoms are dropped, and the tracers keep their row,
        flagged in ``Result.outside_redshift_cut``. mps(z) is measured on the
        tracers kept, and only their redshifts need lie in the distance
        table. None, the default, cuts nothing.
    reject_crossings : bool, optional
        With a mask, filter the result with ``reject_mask_crossings``.
        Ignored without a mask.
    max_unobserved_pixels_crossed : int, optional
        The threshold of that filter. To filter with another threshold, set
        it here: a second filter with the same mask and a higher threshold
        changes nothing, since the filter only marks displacements invalid.
    verbosity : Verbosity, optional
        What the call prints; see ``Verbosity``. With ``"normal"``, when a
        mask or a cut with a finite bound is applied, what was left out: one
        line for the tracers, one for the randoms, one for the rejected
        crossings, as ``Result.selection`` records them; then the line of the
        call. In the form::

            otswap: kept <n> of <N> tracers: <a> outside the redshift cut [<zmin>, <zmax>], <b> outside the mask (<c> outside both)
            otswap: kept <m> of <M> randoms: <a> outside the redshift cut [<zmin>, <zmax>], <b> outside the mask (<c> outside both)
            otswap: rejected <r> of <D> displacements crossing more than <L> unobserved pixels
            otswap: reconstructLightcone: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement

        ``"detailed"`` adds the mps(z) profile before the last line.
    n_realizations, convergence, seed, cell_size
        As in ``reconstruct_box``.
    """


def reject_mask_crossings(
    result: Result,
    mask: Mask,
    max_unobserved_pixels_crossed: int = 0,
) -> None:
    """Mark as invalid the displacements whose path crosses more than
    ``max_unobserved_pixels_crossed`` distinct unobserved pixels of the mask;
    0 rejects at the first one.

    The path is the great-circle arc between the directions of the tracer and
    of its matched random. The pixels containing the two endpoints are not
    tested. The count depends on NSIDE. Updates ``result`` in place; it only
    ever marks displacements invalid, and raises if ``result`` was already
    filtered against a mask of a different NSIDE. Tracers flagged in
    ``outside_redshift_cut`` or ``outside_mask`` are left as they are.
    ``valid_realizations``, ``mean_displacement``, ``lagrangian`` and
    ``lagrangian_sky`` are recomputed, and ``result.selection`` records the
    filter.
    """


def recompute_means(result: Result) -> None:
    """Recompute ``valid_realizations``, ``mean_displacement``, ``lagrangian``
    and ``lagrangian_sky`` from ``displacement`` and ``valid``, as the
    reconstruction computes them: the sum over the valid realizations in
    realization order, divided by their count. Updates ``result`` in place."""


# ---------------------------------------------------------------------------
# Redshift-space correction
# ---------------------------------------------------------------------------

class BiasTable:
    """A b(z) table: the linear bias of the tracers at redshift nodes,
    interpolated linearly between them and extrapolated linearly beyond them
    from the first or last segment. Immutable.

    Read one from a file with ``otswap.io.read_bias_table``."""

    def __init__(self, redshift: ArrayLike, bias: ArrayLike) -> None:
        """At least two nodes, of equal length; redshifts finite and strictly
        increasing; bias finite and positive. Raises otswap.Error otherwise,
        naming the node."""

    @property
    def redshift(self) -> NDArray[np.float64]:
        """Redshifts of the nodes, shape (n,)."""

    @property
    def bias(self) -> NDArray[np.float64]:
        """Bias at each node, shape (n,)."""


class CorrectionStatus(IntEnum):
    """The values of ``RealSpaceCatalog.status``."""

    CORRECTED = 0
    """Had a valid realization; moved."""
    MOVED_BY_NEIGHBOURS = 1
    """Had no valid realization; moved with its neighbours' average."""
    NO_VALID_NEIGHBOUR = 2
    """No valid tracer within 3 sigma (with sigma = 0: no valid realization); not moved."""
    LEFT_OUT = 3
    """Outside the reconstruction's redshift cut or mask; not moved."""


class RealSpaceCatalog:
    """A catalogue moved to real space; row i describes input tracer i.

    Every array attribute is read-only; ``sky`` is in the ``angle_unit`` of
    the result corrected, converted once and kept, and every other array is a
    view on memory owned by this object, as for Result."""

    @property
    def n_objects(self) -> int: ...

    @property
    def geometry(self) -> Literal["box", "lightcone"]: ...

    @property
    def angle_unit(self) -> Optional[AngleUnit]:
        """The unit of ``sky``, that of the result corrected; None for a box."""

    @property
    def sky(self) -> Optional[NDArray[np.float64]]:
        """Lightcone: right ascension, in [0, 360) degrees ([0, 2 pi)
        radians), and declination of each tracer, and its corrected redshift,
        shape (n_objects, 3). The angles are those of ``Result.tracers_sky``;
        in degrees they come back as ``x * (pi/180) * (180/pi)``, which can
        differ from the input in the last bit. NaN rows for the tracers not
        moved. None for a box."""

    @property
    def cartesian(self) -> NDArray[np.float64]:
        """Corrected Cartesian position, Mpc/h, shape (n_objects, 3). Box: the
        tracer with its axis coordinate moved. Lightcone: ``r (1 + shift /
        |r|)``, the tracer moved along its line of sight. NaN rows for the
        tracers not moved."""

    @property
    def shift(self) -> NDArray[np.float64]:
        """Shift applied along the line of sight, Mpc/h, shape (n_objects,):
        the factor times the averaged projection. NaN where nothing moved."""

    @property
    def factor(self) -> NDArray[np.float64]:
        """f / (b + 3 f / 5) at each tracer's observed redshift, shape
        (n_objects,); in a box the single factor, repeated. NaN for the
        tracers left out of the reconstruction."""

    @property
    def status(self) -> NDArray[np.uint8]:
        """What the correction did with each tracer, shape (n_objects,): 0
        corrected, 1 moved with the average of its neighbours (it had no
        valid realization), 2 not moved, with no valid tracer within 3 sigma
        (with sigma = 0: no valid realization), 3 not moved, left out of the
        reconstruction by its redshift cut or mask (``Result.outside_*`` say
        which). See ``CorrectionStatus``."""

    @property
    def valid_realizations(self) -> NDArray[np.uint32]:
        """Valid realizations of each tracer, as in the result corrected,
        shape (n_objects,)."""

    @property
    def n_neighbours(self) -> NDArray[np.uint32]:
        """Tracers with a valid realization averaged for each tracer, itself
        included if valid, shape (n_objects,); 0 for a tracer left out."""

    @property
    def n_realizations_averaged(self) -> NDArray[np.uint32]:
        """Sum of the valid realizations of those tracers, shape (n_objects,)."""

    @property
    def n_extrapolated(self) -> int:
        """Tracers at whose redshift b(z) was extrapolated beyond its table;
        0 for a box."""

    @property
    def sigma(self) -> float:
        """The width of the average, Mpc/h, as given."""

    @property
    def weight_by_realizations(self) -> bool: ...

    @property
    def axis(self) -> Optional[int]:
        """Box: the line-of-sight axis; None for a lightcone."""

    @property
    def box_redshift(self) -> Optional[float]:
        """Box: the redshift of the box; None for a lightcone."""

    @property
    def box_bias(self) -> Optional[float]:
        """Box: the bias of the tracers; None for a lightcone."""

    @property
    def elapsed_seconds(self) -> float:
        """Wall time of the call that made the catalogue, in seconds."""


def radial_projection(displacement: ArrayLike, positions: ArrayLike) -> NDArray[np.float64]:
    """Component of each displacement, shape (N, 3), along the line of sight
    of its position, shape (N, 3), Cartesian, finite and away from the origin
    (for example ``Result.tracers``): ``r . d / |r|``, shape (N,), positive
    when the displacement points away from the observer. NaN displacements
    give NaN; infinite ones raise."""


def axis_projection(displacement: ArrayLike, axis: int) -> NDArray[np.float64]:
    """Component of each displacement, shape (N, 3), along the Cartesian axis
    0, 1 or 2, the line of sight of a box; shape (N,)."""


class NeighbourAverage(NamedTuple):
    """What ``neighbour_average`` returns, shape (N,) each."""

    values: NDArray[np.float64]
    """The average at each object; NaN with no valid object within 3 sigma."""

    n_neighbours: NDArray[np.uint32]
    """Valid objects averaged, the object itself included when it is valid."""

    n_realizations_averaged: NDArray[np.uint32]
    """Sum of the valid realizations of those objects."""


def neighbour_average(
    positions: ArrayLike,
    values: ArrayLike,
    valid_realizations: ArrayLike,
    *,
    sigma: float,
    weight_by_realizations: bool = False,
) -> NeighbourAverage:
    """Gaussian average of ``values`` over the neighbours of each object.

    An object is valid when its entry of ``valid_realizations`` is greater
    than zero. Each object, valid or not, receives the mean of the values of
    the valid objects within 3 sigma of its position, itself included when
    valid, weighted by exp(-d^2 / (2 sigma^2)) and, with
    ``weight_by_realizations``, by their valid realizations; the weights are
    normalised, so a constant value comes back unchanged. NaN when no valid
    object lies within 3 sigma. With sigma = 0, ``values`` is returned
    unchanged.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Cartesian coordinates, finite.
    values : array_like, shape (N,)
        Finite wherever the object is valid.
    valid_realizations : array_like, shape (N,)
        Non-negative integers, such as ``Result.valid_realizations``.
    sigma : float
        Width of the gaussian, in the unit of the positions; finite and
        non-negative. 10 Mpc/h is a reasonable starting value; the best value
        depends on the sample, and should be checked in each analysis.

    Returns
    -------
    NeighbourAverage
        The averages and, per object, the number of valid objects averaged
        and the sum of their valid realizations.

    Notes
    -----
    The result does not depend on the internal neighbour grid within
    rounding, and is the same bit for bit for the same inputs, whatever the
    number of threads. With sigma > 0, positions that span no volume (a
    single object, or all in one plane or on one line) raise RuntimeError
    from the neighbour search, not otswap.Error.
    """


def rsd_factor(
    redshift: ArrayLike,
    distances: DistanceTable,
    *,
    bias: BiasTable,
) -> NDArray[np.float64]:
    """The factor f / (b + 3 f / 5) at each redshift, shape (N,).

    f is ``distances.growth_rate_at(z)``, b(z) the ``bias`` table at z.
    Beyond the table's nodes b(z) is extrapolated linearly from the end
    segments, and one ``ExtrapolationWarning`` is issued for the call. Raises
    if the extrapolated bias is not positive, if a redshift lies outside the
    distance table, or if the table has no growth rate."""


def rsd_factor_box(redshift: float, distances: DistanceTable, *, bias: float) -> float:
    """The factor f / (b + 3 f / 5) of a box at a single redshift, with a
    constant, positive bias."""


def shift_radially(positions: ArrayLike, shift: ArrayLike) -> NDArray[np.float64]:
    """Move each position, shape (N, 3), by its shift, shape (N,), along its
    own line of sight: ``r (1 + shift / |r|)``. A positive shift moves away
    from the observer. A NaN shift gives a NaN row; infinite ones raise."""


def shift_along_axis(positions: ArrayLike, shift: ArrayLike, axis: int) -> NDArray[np.float64]:
    """Move each position, shape (N, 3), by its shift, shape (N,), along the
    Cartesian axis 0, 1 or 2, without periodic wrapping. A NaN shift gives NaN
    along the axis; infinite ones raise."""


def real_space_lightcone(
    result: Result,
    *,
    distances: DistanceTable,
    bias: BiasTable,
    sigma: float,
    weight_by_realizations: bool = False,
    verbosity: Verbosity = "normal",
) -> RealSpaceCatalog:
    """Move a lightcone catalogue, the tracers of ``result``, from redshift
    space to real space.

    For each tracer: the projection of ``result.mean_displacement`` on its line
    of sight, at ``result.tracers``, the gaussian average of width ``sigma``
    over the tracers with a valid realization, times f / (b + 3 f / 5) at its
    observed redshift; the tracer keeps its right ascension and declination
    and moves to the redshift of comoving distance d(z) + shift. A tracer
    without a valid realization still moves with the average of its
    neighbours; only one with no valid tracer within 3 sigma is not moved.

    Parameters
    ----------
    result : Result
        A lightcone result, which carries its tracers; ``sky`` comes back in
        its ``angle_unit``.
    distances : DistanceTable
        Must carry the growth rate, and cover the corrected distances as well
        as the observed ones.
    bias : BiasTable
        The b(z) table, as in ``rsd_factor``; warns ``ExtrapolationWarning``
        once if b(z) is extrapolated.
    sigma : float
        Width of the average, in Mpc/h; 0 for none. 10 Mpc/h is a reasonable
        starting value; the best value depends on the sample, and should be
        checked in each analysis.
    weight_by_realizations : bool, optional
        Weight each neighbour by its number of valid realizations as well.
    verbosity : Verbosity, optional
        What the call prints; see ``Verbosity``. ``"normal"`` prints the line
        of the call; ``"detailed"`` before it how many of the corrected moved
        with their neighbours' average, having no valid realization, and how
        many of the uncorrected were left out of the reconstruction or had no
        valid tracer within 3 sigma.

    Raises
    ------
    Error
        For a box result, or if a corrected comoving distance is not positive
        or lies outside the distance table; the message names the tracer.

    Notes
    -----
    Tracers outside the reconstruction's redshift cut or mask are not
    corrected, take no part in any average, and have status 3.
    """


def real_space_box(
    result: Result,
    *,
    axis: int,
    redshift: float,
    distances: DistanceTable,
    bias: float,
    sigma: float,
    weight_by_realizations: bool = False,
    verbosity: Verbosity = "normal",
) -> RealSpaceCatalog:
    """Move a box catalogue, the tracers of ``result``, from redshift space to
    real space, with the line of sight along ``axis``: the chain of
    ``real_space_lightcone`` with the projection on the axis, a single factor
    at ``redshift`` with the constant ``bias``, and the shift along the axis,
    without periodic wrapping. Raises for a lightcone result.

    The box is not periodic: nothing flows through its faces, so modes on the
    scale of the box itself cannot be reconstructed.
    """
