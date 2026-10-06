"""Optimal transport reconstruction of the displacement field of a tracer
catalogue, by local swapping.

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates are in
Mpc/h; sky coordinates are ordered right ascension, declination, redshift,
with the angular unit given explicitly by ``angle_unit``. Inputs are
converted to float64; the reconstruction works on its own copy of them.

The number of threads follows ``OMP_NUM_THREADS``. The reconstruction
functions release the GIL while they run.
"""

from typing import Literal, Optional, overload

import numpy as np
from numpy.typing import ArrayLike, NDArray

AngleUnit = Literal["deg", "rad"]


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

    Every array attribute is a read-only view on memory owned by this object,
    and keeps it alive. Use ``.copy()`` to obtain a writable array. Filtering
    with ``reject_mask_crossings`` updates ``valid``, ``valid_realizations``
    and ``mean_displacement`` in place, so views taken earlier reflect it.
    """

    @property
    def n_objects(self) -> int: ...

    @property
    def n_realizations(self) -> int: ...

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
    using the same convention as ``reconstruct_lightcone``."""


# ---------------------------------------------------------------------------
# Angular mask
# ---------------------------------------------------------------------------

class Mask:
    """HEALPix mask read from a FITS file. Immutable.

    A pixel is observed when its value is greater than 0, fractional values
    included; 0, negative values, NaN and UNSEEN mark an unobserved pixel.
    No value is used as a weight."""

    def __init__(self, fits_file: str) -> None:
        """NSIDE and ORDERING are read from the file header; RING and NESTED
        are both supported. Any pixel value is accepted."""

    @property
    def nside(self) -> int: ...

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
        [-90, 90] degrees."""


# ---------------------------------------------------------------------------
# Reconstruction
# ---------------------------------------------------------------------------

def reconstruct_box(
    tracers: ArrayLike,
    randoms: Optional[ArrayLike] = None,
    *,
    mps: float,
    n_realizations: int = 1,
    convergence: float = 1e-3,
    seed: int = 0,
    cell_size: float = 4.0,
) -> Result:
    """Reconstruct in box geometry, with a constant mean particle separation.

    Parameters
    ----------
    tracers : (N, 3) Cartesian coordinates.
    randoms : (M, 3) Cartesian coordinates, M >= n_realizations * N; each
        realization uses a disjoint subset. If omitted, n_realizations * N
        randoms are drawn uniformly in the bounding box of the tracers.
    mps : mean particle separation, in the units of the coordinates.
    n_realizations : independent reconstructions to run.
    convergence : the loop stops after a sweep that changes at most this
        fraction of the pairs.
    seed : seed of the random streams; 0 draws one at random. With a fixed
        seed the result does not depend on the number of threads.
    cell_size : grid cell size in units of the mps. Affects speed only.

    The box is not periodic: nothing flows through its faces, so modes on the
    scale of the box itself cannot be reconstructed.
    """


def reconstruct_lightcone(
    tracers_sky: ArrayLike,
    randoms_sky: ArrayLike,
    *,
    sky_area_deg2: float,
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
) -> Result:
    """Reconstruct in lightcone geometry.

    The mean particle separation is measured from the tracers: n_bins
    uniform redshift bins over the observed range, mps = (N / V)^(-1/3) in
    each, with V the shell volume implied by ``sky_area_deg2``. The bin
    values are joined by linear interpolation and extrapolated linearly
    beyond the outermost bins. Raises if a bin holds too few tracers for a
    2% accuracy on its mps.

    Randoms are required: they carry the survey geometry, the selection
    function and the completeness.

    Parameters
    ----------
    tracers_sky : (N, 3) sky coordinates of the tracers.
    randoms_sky : (M, 3) sky coordinates of the randoms,
        M >= n_realizations * N.
    sky_area_deg2 : effective survey area. ``Mask.sky_area_deg2`` gives it
        when a mask is available.
    n_bins : redshift bins used to measure mps(z).
    distances : table used for the conversion to Cartesian coordinates.
    angle_unit : unit of right ascension and declination.
    tracers, randoms : Cartesian coordinates already computed. Give both or
        neither; when given, the conversion is skipped, and their agreement
        with the sky coordinates is not checked.
    redshift_cut : (min, max). Tracers and randoms whose redshift lies outside
        this closed range are left out of the reconstruction, before anything
        else: the randoms are dropped, and the tracers keep their row, flagged
        in ``Result.outside_redshift_cut``. mps(z) is measured on the tracers
        kept, and only their redshifts need lie in the distance table. None,
        the default, cuts nothing.

    The remaining parameters are as in ``reconstruct_box``.
    """


def reject_mask_crossings(
    result: Result,
    mask: Mask,
    max_forbidden_pixels: int = 0,
) -> None:
    """Mark as invalid the displacements whose path crosses more than
    ``max_forbidden_pixels`` distinct unobserved pixels of the mask; 0
    rejects at the first one.

    The path is the great-circle arc between the directions of the tracer and
    of its matched random. The pixels containing the two endpoints are not
    tested. The count depends on NSIDE. Updates ``result`` in place; it only
    ever marks displacements invalid, and raises if ``result`` was already
    filtered against a mask of a different NSIDE.
    """


# ---------------------------------------------------------------------------
# Redshift-space correction
#
# The correction of tracer i is a shift along its line of sight,
#
#     s_i = f(z_i) / (b(z_i) + 3 f(z_i) / 5) * <Psi . r_hat>_i ,
#
# with Psi the reconstructed displacement (Result.mean_displacement, from the
# observed to the reconstructed position), r_hat the line of sight, f the
# linear growth rate, b the linear bias and < > a gaussian average over the
# neighbouring tracers. real_space_lightcone and real_space_box run the whole
# chain; the four steps are also available on their own.
# ---------------------------------------------------------------------------

class RealSpaceCatalog:
    """A catalogue moved to real space; row i describes input tracer i.

    Every array attribute is a read-only view on memory owned by this object,
    as for Result."""

    @property
    def n_objects(self) -> int: ...

    @property
    def positions(self) -> NDArray[np.float64]:
        """Shape (n_objects, 3). Lightcone: right ascension and declination as
        given, in the angle_unit of the call, and the corrected redshift. Box:
        the Cartesian position, corrected along the axis. NaN rows for the
        tracers in ``uncorrected``."""

    @property
    def n_neighbours(self) -> NDArray[np.uint32]:
        """Tracers with a valid realization averaged for each tracer, itself
        included if valid, shape (n_objects,)."""

    @property
    def n_realizations_averaged(self) -> NDArray[np.uint32]:
        """Sum of the valid realizations of those tracers, shape (n_objects,)."""

    @property
    def uncorrected(self) -> NDArray[np.int64]:
        """Indices, increasing, of the tracers left without a correction: those
        with no tracer with a valid realization within 3 sigma, themselves
        included. With sigma = 0, those without a valid realization."""


def line_of_sight_projection(
    displacement: ArrayLike,
    *,
    positions: Optional[ArrayLike] = None,
    axis: Optional[int] = None,
) -> NDArray[np.float64]:
    """Component of each displacement, shape (N, 3), along its line of sight,
    shape (N,): positive when it points away from the observer.

    Give exactly one of ``positions``, Cartesian (N, 3), finite and away from
    the origin, for a radial line of sight, and ``axis`` (0, 1 or 2) for a
    box. NaN displacements give NaN; infinite ones raise."""


@overload
def neighbour_average(
    positions: ArrayLike,
    values: ArrayLike,
    valid_realizations: ArrayLike,
    *,
    sigma: float,
    weight_by_realizations: bool = False,
    diagnostics: Literal[False] = False,
) -> NDArray[np.float64]: ...
@overload
def neighbour_average(
    positions: ArrayLike,
    values: ArrayLike,
    valid_realizations: ArrayLike,
    *,
    sigma: float,
    weight_by_realizations: bool = False,
    diagnostics: Literal[True],
) -> tuple[NDArray[np.float64], NDArray[np.uint32], NDArray[np.uint32]]:
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
    positions : (N, 3) Cartesian coordinates, finite.
    values : (N,); finite wherever the object is valid.
    valid_realizations : (N,) non-negative integers, such as
        ``Result.valid_realizations``.
    sigma : width of the gaussian, in the unit of the positions; finite and
        non-negative. 10 Mpc/h is a reasonable starting value; the best value
        depends on the sample, and should be checked in each analysis.
    diagnostics : also return, per object, the number of valid objects
        averaged and the sum of their valid realizations.

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
    bias_redshift: ArrayLike,
    bias: ArrayLike,
) -> NDArray[np.float64]:
    """The factor f / (b + 3 f / 5) at each redshift, shape (N,).

    f is ``distances.growth_rate_at(z)``. b(z) is interpolated linearly
    between the nodes (``bias_redshift``, ``bias``): at least two, redshifts
    strictly increasing, bias finite and positive. Beyond the nodes it is
    extrapolated linearly from the end segments, and one
    ``ExtrapolationWarning`` is issued for the call. Raises if the
    extrapolated bias is not positive, if a redshift lies outside the
    distance table, or if the table has no growth rate."""


def rsd_factor_box(redshift: float, distances: DistanceTable, *, bias: float) -> float:
    """The factor f / (b + 3 f / 5) of a box at a single redshift, with a
    constant, positive bias."""


def shift_along_line_of_sight(
    positions: ArrayLike,
    shift: ArrayLike,
    *,
    axis: Optional[int] = None,
) -> NDArray[np.float64]:
    """Move each position, shape (N, 3), by its shift, shape (N,), along its
    own line of sight, r (1 + shift / |r|), or, with ``axis``, along that
    Cartesian axis without periodic wrapping. A positive shift moves away from
    the observer. A NaN shift gives NaN; infinite ones raise."""


def real_space_lightcone(
    result: Result,
    tracers_sky: ArrayLike,
    *,
    distances: DistanceTable,
    bias_redshift: ArrayLike,
    bias: ArrayLike,
    sigma: float,
    angle_unit: AngleUnit,
    weight_by_realizations: bool = False,
) -> RealSpaceCatalog:
    """Move a lightcone catalogue from redshift space to real space.

    For each tracer: the projection of ``result.mean_displacement`` on its line
    of sight, the gaussian average of width ``sigma`` over the tracers with a
    valid realization, times f / (b + 3 f / 5) at its observed redshift; the
    tracer keeps its right ascension and declination and moves to the redshift
    of comoving distance d(z) + shift. A tracer without a valid realization
    still moves with the average of its neighbours; only one with no valid
    tracer within 3 sigma is left uncorrected.

    Parameters
    ----------
    result : the reconstruction of these tracers. Only their number is
        checked: pass the arrays the reconstruction was run on, in the same
        order.
    tracers_sky : (N, 3) observed sky coordinates.
    distances : must carry the growth rate, and cover the corrected
        distances as well as the observed ones.
    bias_redshift, bias : the b(z) table, as in ``rsd_factor``; warns
        ``ExtrapolationWarning`` once if b(z) is extrapolated.
    sigma : width of the average, in Mpc/h; 0 for none. 10 Mpc/h is a
        reasonable starting value; the best value depends on the sample, and
        should be checked in each analysis.
    weight_by_realizations : weight each neighbour by its number of valid
        realizations as well.

    Tracers outside the reconstruction's redshift cut are not corrected,
    take no part in any average, and are listed in ``uncorrected``.

    Raises if a corrected comoving distance is not positive or lies outside
    the distance table; the message names the tracer.
    """


def real_space_box(
    result: Result,
    tracers: ArrayLike,
    *,
    axis: int,
    redshift: float,
    distances: DistanceTable,
    bias: float,
    sigma: float,
    weight_by_realizations: bool = False,
) -> RealSpaceCatalog:
    """Move a box catalogue from redshift space to real space, with the line of
    sight along ``axis``: the chain of ``real_space_lightcone`` with the
    projection on the axis, a single factor at ``redshift`` with the constant
    ``bias``, and the shift along the axis, without periodic wrapping.

    The box is not periodic: nothing flows through its faces, so modes on the
    scale of the box itself cannot be reconstructed.
    """
