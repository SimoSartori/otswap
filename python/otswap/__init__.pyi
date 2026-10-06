"""Optimal transport reconstruction of the displacement field of a tracer
catalogue, by local swapping.

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates are in
Mpc/h; sky coordinates are ordered right ascension, declination, redshift,
with the angular unit given explicitly by ``angle_unit``. Inputs are
converted to float64; the reconstruction works on its own copy of them.

The number of threads follows ``OMP_NUM_THREADS``. The reconstruction
functions release the GIL while they run.
"""

from typing import Literal, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

AngleUnit = Literal["deg", "rad"]


class Error(RuntimeError):
    """Raised by every otswap function on invalid input or failure."""


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
