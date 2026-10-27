"""Optimal transport reconstruction of the displacement field of a tracer
catalogue, by local swapping.

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates, and
every length, are in coordinate units: the unit of the tracers given in a box,
that of the distance table in a lightcone (Mpc/h from
``DistanceTable.flat``); sky coordinates are ordered right ascension, declination, redshift,
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
from math import inf
from typing import ClassVar, Literal, NamedTuple, Optional, Union

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
  applied, and of a box, when ``box`` was given; then one line per call with
  its time and what was lost::

      otswap: reconstructBox: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement
      otswap: reconstructLightcone: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement
      otswap: measureMpsProfile: <B> bins over [<zmin>, <zmax>]
      otswap: realSpaceLightcone: corrected <C> of <N> tracers in <T> s; <U> left uncorrected
      otswap: realSpaceBox: corrected <C> of <N> tracers in <T> s; <U> left uncorrected

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
# Configuration
# ---------------------------------------------------------------------------

class Config:
    """Free parameters of the reconstruction, as ``otswap::Config``: the
    ``config`` of ``reconstruct_box`` and ``reconstruct_lightcone``.

    Every field is a keyword of the constructor, with the default of the C++
    struct, and a read/write attribute. A value is converted to the field's
    C++ type when it is stored, and raises ``Error`` if it cannot be; its
    rules (``n_realizations`` at least 1, ...) are checked by the call that
    reads it, which refuses it with the library's message. Each call works
    on its own copy, so the object passed is never changed. ``repr`` gives
    the constructor call; two configs compare equal when every field is
    equal. Unhashable, being mutable.
    """

    def __init__(
        self,
        *,
        n_realizations: int = 30,
        convergence: float = 1e-3,
        seed: int = 0,
        cell_size: float = 4.0,
        reject_crossings: bool = True,
        max_unobserved_pixels_crossed: int = 0,
        verbosity: Verbosity = "normal",
    ) -> None: ...

    n_realizations: int
    """Independent reconstructions to run."""

    convergence: float
    """The sweeps stop once the fraction of successful swaps in a sweep,
    swaps per tracer visited, no longer exceeds this threshold."""

    seed: int
    """Seed of the random streams; 0 draws one at random. With a fixed seed
    the result does not depend on the number of threads."""

    cell_size: float
    """Grid cell size in units of the mps. Affects speed only."""

    reject_crossings: bool
    """With a mask given to ``reconstruct_lightcone``, filter the result with
    ``reject_mask_crossings``. Read by ``reconstruct_lightcone`` with a mask
    only."""

    max_unobserved_pixels_crossed: int
    """The threshold of that filter. To filter with another threshold, set
    it here: a second filter with the same mask and a higher threshold
    changes nothing, since the filter only marks displacements invalid. Read
    by ``reconstruct_lightcone`` with a mask only."""

    verbosity: Verbosity
    """What the call prints; see ``Verbosity`` and each function."""

    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]


class CorrectionConfig:
    """Options of ``real_space_lightcone`` and ``real_space_box``, as
    ``otswap::CorrectionConfig``: their ``correction``. Keywords,
    attributes, ``repr`` and equality as ``Config``.
    """

    def __init__(
        self,
        *,
        weight_by_realizations: bool = False,
        verbosity: Verbosity = "normal",
    ) -> None: ...

    weight_by_realizations: bool
    """Weight each neighbour by its number of valid realizations as well."""

    verbosity: Verbosity
    """What the call prints; see ``Verbosity`` and each function."""

    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]


class RedshiftCut:
    """A closed redshift range for ``reconstruct_lightcone`` and
    ``measure_mps_profile``, as ``otswap::RedshiftCut``. The default,
    (-inf, +inf), cuts nothing. Wherever a RedshiftCut is taken, a
    ``(min, max)`` pair is accepted too.

    ``min`` and ``max`` are read/write; the call that reads the cut checks
    them (neither NaN, ``min`` not above ``max``). ``repr`` and equality as
    ``Config``.
    """

    def __init__(self, min: float = -inf, max: float = inf) -> None: ...

    min: float
    """The lowest redshift kept."""

    max: float
    """The highest redshift kept."""

    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------

class Result:
    """Displacement field produced by a reconstruction.

    Displacements point from the observed (Eulerian) position of a tracer to
    its reconstructed (Lagrangian) position, which is the position of the
    random it is matched to.

    Every array attribute is read-only, and only ``ids`` can be assigned.
    The sky arrays (``tracers_sky``, ``lagrangian_sky``) are in the
    ``angle_unit`` of the call that made the result; in degrees they are
    converted once, on first access, and kept. ``ids`` is a new copy on each
    access. Every other array is a view on memory owned by this object, and
    keeps it alive. Use ``.copy()`` to obtain a writable array. Filtering with
    ``reject_mask_crossings`` updates ``valid``, ``valid_realizations``,
    ``mean_displacement`` and ``lagrangian`` in place, so views taken earlier
    reflect it; ``lagrangian_sky`` is recomputed, and an array taken from it
    in degrees before the filter is not updated.

    The type is ``otswap::Result`` itself, so another nanobind module can take
    and return it (see "Use from another nanobind module" on the "Names side
    by side" page of the documentation). The angle unit is kept by the module,
    not in the object: a result made by another module is in radians. Each
    array getter checks the size of the C++ array first, and raises
    ``otswap.Error`` for a malformed result, which only another module can
    make. No attribute can be added to a result.
    """

    @staticmethod
    def from_arrays(
        displacement: ArrayLike,
        valid: ArrayLike,
        tracers: ArrayLike,
        *,
        tracers_sky: Optional[ArrayLike] = None,
        angle_unit: Optional[AngleUnit] = None,
        distances: Optional["DistanceTable"] = None,
        outside_redshift_cut: Optional[ArrayLike] = None,
        outside_mask: Optional[ArrayLike] = None,
        outside_box: Optional[ArrayLike] = None,
        ids: Optional[ArrayLike] = None,
        config: Optional[Config] = None,
    ) -> "Result":
        """Build a result from its arrays, as read back from a file.

        Parameters
        ----------
        displacement : array_like, shape (n_realizations, n_objects, 3)
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
        outside_box : array_like, shape (n_objects,), optional
            0 or 1; all 0 when omitted. For a box result only: raises with
            ``tracers_sky``.
        ids : array_like, shape (n_objects,), optional
            The IDs, checked and stored as by assigning ``Result.ids``; none
            when omitted.
        config : Config, optional
            The configuration to record, as ``Result.config``; its
            ``n_realizations`` must equal the arrays'. When omitted, the
            defaults with the arrays' ``n_realizations``.

        Notes
        -----
        ``valid_realizations``, ``mean_displacement``, ``lagrangian`` and
        ``lagrangian_sky`` are computed as ``recompute_means`` computes them,
        so arrays read back exactly give the reconstruction's values, bit for
        bit. The profile, the mps, the box and the selection counts are not
        reconstructed: ``mps``, ``mps_profile`` and ``box`` are None.
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
        ``angle_unit`` of the call that made the result; ``"rad"`` for a
        result made by another module. None for a box."""

    @property
    def config(self) -> Config:
        """The parameters of the run: a new copy on each access, so changing
        it changes nothing in the result. ``seed`` is the one used: the one
        drawn when 0 was given, so that the result can be reproduced."""

    @property
    def distances(self) -> Optional["DistanceTable"]:
        """The distance table of a lightcone result, a new copy on each
        access; None for a box."""

    @property
    def ids(self) -> NDArray[np.int64]:
        """The ID of each object, shape (n_objects,), read-only; when none
        were set, 0..n_objects-1, the row of each object in the input.

        Assignable, before writing or correcting the result: one 64-bit
        integer per object, negative ones allowed, no value twice; a float
        array is accepted when every entry is an integer. None, an empty
        array or exactly 0..n_objects-1 set no IDs. Any other value raises
        ``Error``, naming the first repeated value and its two entries, and
        leaves the IDs unchanged.

        The writers write them as the id column, among their default
        columns when IDs were set; ``real_space_box`` and
        ``real_space_lightcone`` copy them into the ``RealSpaceCatalog``. No
        computation reads them, and the reconstruction sets none.

        Each access returns a new copy, so store it once (``ids =
        result.ids``) before indexing it in a loop."""

    @ids.setter
    def ids(self, value: Optional[ArrayLike]) -> None: ...

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
        coordinate units. None for a lightcone, whose separation is ``mps_profile``."""

    @property
    def box(self) -> Optional["BoxLimits"]:
        """Box: the box of the reconstruction, the ``box`` given to
        ``reconstruct_box`` or, without it, the bounding box of the tracers
        and the given randoms (of the tracers alone when the randoms were
        drawn), also when ``mps`` was given; every random matched to a tracer
        lies inside it; a copy. None for a lightcone."""

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
        displacements, no valid realization and a NaN mean displacement."""

    @property
    def outside_mask(self) -> NDArray[np.bool_]:
        """Whether each tracer fell on an unobserved pixel of the ``mask`` of
        ``reconstruct_lightcone`` and took no part in the reconstruction,
        shape (n_objects,). All false without a mask. Independent of
        ``outside_redshift_cut``: a tracer may carry both flags. A flagged
        tracer's rows are as for a cut one."""

    @property
    def outside_box(self) -> NDArray[np.bool_]:
        """Whether each tracer lay outside the ``box`` given to
        ``reconstruct_box`` and took no part in the reconstruction, shape
        (n_objects,). All false without a box and for a lightcone. A flagged
        tracer's rows are as for a cut one."""

    @property
    def selection(self) -> "SelectionCounts":
        """What the reconstruction left out, and why. Keeps this result
        alive."""

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
    redshift, measured by ``reconstruct_lightcone``, or alone by
    ``measure_mps_profile``: the tracers kept in
    ``n_bins`` uniform redshift bins over [``redshift_min``,
    ``redshift_max``], mps = (N / V)^(-1/3) in each, joined linearly and
    extrapolated linearly beyond the outermost bins. Read-only."""

    @property
    def redshift(self) -> NDArray[np.float64]:
        """Bin centres, shape (n_bins,)."""

    @property
    def mps(self) -> NDArray[np.float64]:
        """Mean particle separation at each centre, in coordinate units, shape
        (n_bins,)."""

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
        """The mps at each redshift, in coordinate units, the shape of ``z``
        kept: the value
        the reconstruction used at a tracer of that redshift. Raises if the
        profile is empty, or where the extrapolation is not positive."""


class SelectionCounts:
    """What a reconstruction left out, and why. Read-only.

    For ``reconstruct_lightcone``, the tracers removed are
    ``tracers_outside_redshift_cut + tracers_outside_mask -
    tracers_outside_both``, and the same for the randoms. For
    ``reconstruct_box``, the tracers removed are ``tracers_outside_box``, and
    the randoms removed ``randoms_outside_box``; these, ``tracers`` and
    ``randoms`` are set only when a ``box`` was given, ``randoms`` 0 when the
    randoms were drawn. Every other count is 0, and nothing else is applied,
    for ``reconstruct_box``; the box counts are 0 for
    ``reconstruct_lightcone``. ``repr`` gives the report the reconstruction
    prints.
    """

    @property
    def redshift_cut(self) -> Optional[RedshiftCut]:
        """The redshift cut applied, a new ``RedshiftCut`` on each access;
        None when no cut was applied."""

    @property
    def mask_applied(self) -> bool:
        """Whether a mask was given."""

    @property
    def box_limits(self) -> Optional["BoxLimits"]:
        """The ``box`` given to ``reconstruct_box``; None when none was."""

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
    def tracers_outside_box(self) -> int:
        """Tracers outside the box limits."""

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
    def randoms_outside_box(self) -> int:
        """Randoms given outside the box limits."""

    @property
    def max_unobserved_pixels_crossed(self) -> Optional[int]:
        """The threshold of the latest mask filter applied to the result, by
        ``reconstruct_lightcone`` (``config.reject_crossings``) or by
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

        Raises
        ------
        Error
            If any parameter is out of range: ``omega_m`` must lie in (0, 1],
            ``h`` be positive, ``w0`` and ``wa`` finite, ``z_min``
            non-negative, ``z_max`` above ``z_min``, and ``n_samples`` at
            least two.
        """

    @staticmethod
    def from_table(
        redshift: ArrayLike,
        distance: ArrayLike,
        growth_rate: Optional[ArrayLike] = None,
    ) -> "DistanceTable":
        """Use a table computed by the caller. Redshift and distance must be
        strictly increasing, with at least two entries; distance in the unit
        that Cartesian coordinates and every length will then have.
        ``growth_rate``, when given, has the same length, every entry finite
        and positive. Without ``growth_rate``, ``growth_rate_at`` raises.

        Raises
        ------
        Error
            If either of the first two arrays is not strictly increasing, if
            the lengths differ, if fewer than two points are given, if an
            entry is not finite, or if ``growth_rate`` is neither None nor of
            the same length, or holds a value that is not finite and
            positive.
        """

    @property
    def min_redshift(self) -> float: ...

    @property
    def max_redshift(self) -> float: ...

    @property
    def has_growth_rate(self) -> bool: ...

    @property
    def redshift(self) -> NDArray[np.float64]:
        """The sampled redshifts, shape (n_samples,), strictly increasing; a
        read-only view. ``DistanceTable.from_table(t.redshift, t.distance,
        t.growth_rate)`` gives the same values, bit for bit, at every
        lookup."""

    @property
    def distance(self) -> NDArray[np.float64]:
        """The comoving distance at each sampled redshift, Mpc/h from ``flat``,
        the unit of the table given to ``from_table``; a read-only view."""

    @property
    def growth_rate(self) -> Optional[NDArray[np.float64]]:
        """The growth rate at each sampled redshift, a read-only view; None
        when the table carries none."""

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
    """Convert sky coordinates, shape (N, 3), to Cartesian ones in the
    unit of the distance table,
    using the same convention as ``reconstruct_lightcone``. Raises if the
    array is not of shape (N, 3), if an entry is not finite, if a declination
    lies outside [-90, 90] degrees ([-pi/2, pi/2] radians), or if a redshift
    falls outside the distance table; the message names the object and gives
    the angle in ``angle_unit``."""


def to_sky(
    cartesian: ArrayLike,
    distances: DistanceTable,
    *,
    angle_unit: AngleUnit,
) -> NDArray[np.float64]:
    """Convert Cartesian coordinates, shape (N, 3), in the distance table's
    unit, to sky
    coordinates, shape (N, 3): right ascension in [0, 360) degrees
    ([0, 2 pi) radians), declination, redshift. The inverse of
    ``to_cartesian``, with the same observer; a position at the origin gets
    right ascension and declination 0. A row holding a NaN gives a NaN
    row. Raises if the array is not of shape (N, 3), if an entry is
    infinite, or if a distance falls outside the table; the message names the
    object and, for the distance, suggests a table covering a wider redshift
    range."""


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
        bytes only. A null value of an integer column is unobserved. Raises
        if the file cannot be opened or read, if it is not a full-sky HEALPix
        map in that layout, if NSIDE exceeds 2**29, or if the ordering is
        NESTED and NSIDE is not a power of 2."""

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

class BoxLimits:
    """The closed limits of a box for ``reconstruct_box``: the points with
    ``lower[k] <= x[k] <= upper[k]`` on each axis k, in coordinate units.
    Read-only.

    Wherever a BoxLimits is taken, a pair of corners
    ``((x0, y0, z0), (x1, y1, z1))`` is accepted too. Two box limits compare
    equal when both corners are equal; unhashable.

    Parameters
    ----------
    lower, upper : array_like, shape (3,)
        x, y, z of the lower and of the upper corner. Only their shape is
        checked here; ``reconstruct_box`` checks that each bound is finite and
        that ``lower[k] < upper[k]``.
    """

    def __init__(self, lower: ArrayLike, upper: ArrayLike) -> None: ...

    @property
    def lower(self) -> tuple[float, float, float]:
        """x, y, z of the lower corner."""

    @property
    def upper(self) -> tuple[float, float, float]:
        """x, y, z of the upper corner."""

    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]


def reconstruct_box(
    tracers: ArrayLike,
    randoms: Optional[ArrayLike] = None,
    *,
    mps: Optional[float] = None,
    box: Optional[Union[BoxLimits, tuple[ArrayLike, ArrayLike]]] = None,
    config: Config = Config(),
) -> Result:
    """Reconstruct in box geometry, with a constant mean particle separation.

    Parameters
    ----------
    tracers : array_like, shape (N, 3)
        Cartesian coordinates.
    randoms : array_like, shape (M, 3), optional
        Cartesian coordinates, M >= R * N with R ``config.n_realizations``,
        both counted inside ``box`` when it is given; each realization uses a
        disjoint subset. If omitted, R * N randoms are drawn uniformly in
        ``box``, or without it in the bounding box of the tracers.
    mps : float, optional
        Mean particle separation, in the units of the coordinates. If
        omitted, (V / N)^(1/3), with N the tracers inside the box and V its
        volume, the box being ``Result.box``: ``box`` or, without it, the
        bounding box of the tracers and the given randoms, or of the tracers
        alone when the randoms are drawn (the box drawn randoms fill); it is
        ``Result.mps``. Without ``box``, raises if the tracers, with the given
        randoms, span no volume.
    box : BoxLimits or ((x0, y0, z0), (x1, y1, z1)), optional
        The box of the reconstruction, closed: a point on a face is inside.
        The tracers outside it take no part in the reconstruction and are
        flagged in ``Result.outside_box``, the randoms given outside it are
        dropped, and the neighbour grids span it. ``Result.box`` holds it.
    config : Config, optional
        ``n_realizations``, ``convergence``, ``seed``, ``cell_size`` and
        ``verbosity``; ``reject_crossings`` and
        ``max_unobserved_pixels_crossed`` are not read. With
        ``n_realizations`` R, the randoms given must number at least R N.
        ``verbosity``: see ``Verbosity``; ``"normal"`` prints one line,
        ``"detailed"`` the mps and its source before it. With ``box``,
        ``"normal"`` prints first what the box left out::

            otswap: kept <n> of <N> tracers: <a> outside the box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]
            otswap: kept <m> of <M> randoms: <b> outside the box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]

        the second line only when randoms are given.

    Raises
    ------
    Error
        If an array is not of shape (N, 3), if the tracer array is empty, if
        the randoms are too few, if any coordinate is not finite, or if
        ``mps`` is not positive; without ``mps`` and ``box``, also if the
        tracers span no volume (all of them in one plane, on one line or at
        one point).
    Error
        If a field of ``config`` is out of range, the message describing the
        quantity and its rule.
    Error
        If a bound of ``box`` is not finite or a lower bound is not below the
        upper one, or if ``box`` keeps too few tracers or randoms; the
        message names the box.
    Error
        From the neighbour search, with a message starting
        ``"meshsearch: "``, if its grid cannot be built: a cell of
        ``config.cell_size`` mps so small against the extent of the tracers and
        randoms that the grid would need more cells than it can count.

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
    redshift_cut: Optional[Union[RedshiftCut, tuple[float, float]]] = None,
    config: Config = Config(),
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
        Sky coordinates of the randoms. Their number is checked after the
        redshift cut and the mask: the randoms kept must number at least
        ``config.n_realizations`` times the tracers kept, which is
        M >= R * N, with R ``config.n_realizations``, when nothing is left
        out.
    sky_area_deg2 : float, optional
        Effective survey area, when no mask is given: the published value.
        Give exactly one of ``sky_area_deg2`` and ``mask``.
    mask : Mask, optional
        The survey's mask. The sky area is then ``mask.sky_area_deg2``, and
        the mask is applied twice. Before the reconstruction, tracers and
        randoms on unobserved pixels are left out, as ``Mask.allows`` decides:
        the randoms are dropped, and the tracers keep their row, flagged in
        ``Result.outside_mask``. After it, with ``config.reject_crossings``,
        ``reject_mask_crossings`` is applied with
        ``config.max_unobserved_pixels_crossed``. The mask and the redshift
        cut are both evaluated on every object; the objects kept pass both.
    n_bins : int
        Redshift bins used to measure mps(z).
    distances : DistanceTable
        The redshift to distance relation, with which the sky coordinates
        are converted and the shell volumes computed; it must cover the
        redshifts of the objects kept.
    angle_unit : AngleUnit
        Unit of right ascension and declination. A declination outside
        [-90, 90] degrees ([-pi/2, pi/2] radians) raises, the message giving
        it in this unit.
    tracers, randoms : array_like, optional
        Cartesian coordinates already computed, shape (N, 3) and (M, 3). Give
        both or neither; when given, the conversion is skipped, and their
        agreement with the sky coordinates is not checked.
    redshift_cut : RedshiftCut or (min, max), optional
        Tracers and randoms whose redshift lies outside this closed range
        are left out of the reconstruction, before anything else: the
        randoms are dropped, and the tracers keep their row, flagged in
        ``Result.outside_redshift_cut``. mps(z) is measured on the tracers
        kept, and only their redshifts need lie in the distance table. None,
        the default, and ``RedshiftCut()`` cut nothing.
    config : Config, optional
        Every field; ``reject_crossings`` and
        ``max_unobserved_pixels_crossed`` are read only with a mask.
        ``verbosity``: see ``Verbosity``. With ``"normal"``, when a mask or a
        cut with a finite bound is applied, what was left out: one line for
        the tracers, one for the randoms, and one for the rejected crossings
        when the mask filter was applied, as ``Result.selection`` records
        them; then the line of the call. In the form::

            otswap: kept <n> of <N> tracers: <a> outside the redshift cut [<zmin>, <zmax>], <b> outside the mask (<c> outside both)
            otswap: kept <m> of <M> randoms: <a> outside the redshift cut [<zmin>, <zmax>], <b> outside the mask (<c> outside both)
            otswap: rejected <r> of <D> displacements crossing more than <L> unobserved pixels
            otswap: reconstructLightcone: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement

        ``"detailed"`` adds the mps(z) profile before the last line.

    Raises
    ------
    Error
        If an array is not of shape (N, 3), if ``tracers`` is given without
        ``randoms``, or if they describe a different number of objects than
        the sky arrays; if a declination is out of range, or a redshift falls
        outside the distance table; if the randoms are too few; if both or
        neither of ``sky_area_deg2`` and ``mask`` are given.
    Error
        If a bound of ``redshift_cut`` is NaN or its minimum exceeds its
        maximum, or if too few tracers or randoms remain inside the cut and
        the mask; the message gives the counts kept and dropped.
    Error
        If a bin holds too few tracers for its density to be meaningful, the
        message reporting how many bins the catalogue supports, or if the
        profile extrapolates to an mps that is not positive at a tracer's
        redshift.
    Error
        If a field of ``config`` is out of range, as in ``reconstruct_box``.
    Error
        From the neighbour search, as in ``reconstruct_box``, the cell being
        ``config.cell_size`` mean particle separations.
    """


def measure_mps_profile(
    tracers_sky: ArrayLike,
    *,
    sky_area_deg2: Optional[float] = None,
    mask: Optional[Mask] = None,
    n_bins: int,
    distances: DistanceTable,
    redshift_cut: Optional[Union[RedshiftCut, tuple[float, float]]] = None,
    angle_unit: AngleUnit,
    verbosity: Verbosity = "normal",
) -> MpsProfile:
    """Measure the mps(z) profile of a lightcone without reconstructing.

    The steps of ``reconstruct_lightcone`` that concern the tracers, in the
    same order and with the same code: the checks of the sky coordinates
    and of the cut, the redshift cut and the mask, the minimum number of
    tracers kept, and the profile with its per-bin threshold. The profile
    equals, bit for bit, ``Result.mps_profile`` of a
    ``reconstruct_lightcone`` call with the same tracers, sky area or mask,
    ``n_bins``, ``distances`` and ``redshift_cut``, whatever its randoms and
    other arguments. Use it to choose ``n_bins``, or to inspect the
    profile, before a reconstruction.

    Parameters
    ----------
    tracers_sky : array_like, shape (N, 3)
        Sky coordinates of the tracers.
    sky_area_deg2 : float, optional
        Effective survey area, when no mask is given, as in
        ``reconstruct_lightcone``. Give exactly one of ``sky_area_deg2`` and
        ``mask``.
    mask : Mask, optional
        The survey's mask. The sky area is then ``mask.sky_area_deg2``, and
        the tracers on unobserved pixels are left out, as
        ``reconstruct_lightcone`` leaves them out.
    n_bins : int
        Redshift bins used to measure mps(z).
    distances : DistanceTable
        The redshift to distance relation; it must cover the redshifts of
        the tracers kept.
    redshift_cut : RedshiftCut or (min, max), optional
        As in ``reconstruct_lightcone``: the tracers whose redshift lies
        outside this closed range are left out. None, the default, cuts
        nothing.
    angle_unit : AngleUnit
        Unit of right ascension and declination. A declination outside
        [-90, 90] degrees ([-pi/2, pi/2] radians) raises, the message giving
        it in this unit.
    verbosity : Verbosity, optional
        What the call prints; see ``Verbosity``. With ``"normal"``, when a
        mask or a cut with a finite bound is applied, the tracer line of the
        selection report of ``reconstruct_lightcone``; then the line of the
        call. In the form::

            otswap: kept <n> of <N> tracers: <a> outside the redshift cut [<zmin>, <zmax>], <b> outside the mask (<c> outside both)
            otswap: measureMpsProfile: <B> bins over [<zmin>, <zmax>]

        ``"detailed"`` adds the mps(z) profile before the last line.

    Returns
    -------
    MpsProfile
        The profile.

    Raises
    ------
    Error
        As ``reconstruct_lightcone``, with the same messages, if the tracer
        array is malformed, if a declination is out of range, if a bound of
        the cut is NaN or its minimum exceeds its maximum, if a cut or the
        mask keeps too few tracers, if a bin holds too few tracers, or if a
        redshift kept falls outside the distance table; also if both or
        neither of ``sky_area_deg2`` and ``mask`` are given.
    """


def reject_mask_crossings(
    result: Result,
    mask: Mask,
    max_unobserved_pixels_crossed: int = 0,
) -> None:
    """Mark as invalid the displacements whose path crosses more than
    ``max_unobserved_pixels_crossed`` distinct unobserved pixels of the mask;
    0 rejects at the first one.

    The path is the great-circle arc between the directions of the tracer
    (``result.tracers``) and of the tracer plus its displacement. The pixels
    containing the two endpoints are not tested. The count depends on NSIDE.
    Updates ``result`` in place; it only ever marks displacements invalid, and
    raises if ``result`` was already filtered against a mask of a different
    NSIDE, or if ``result`` is a box result, which has no sky to filter
    against. It also raises if ``result`` is malformed: a NaN outside the
    rows of the flagged tracers, a flagged tracer with a valid entry, a
    ``tracers`` array without N rows, a flag, ``tracers_sky``,
    ``lagrangian`` or ``lagrangian_sky`` array whose size is neither 0 nor
    the one of its layout, or an entry of ``valid`` neither 0 nor 1. Tracers
    flagged in ``outside_redshift_cut`` or ``outside_mask`` are left as they
    are. ``valid_realizations``, ``mean_displacement``, ``lagrangian`` and
    ``lagrangian_sky`` are recomputed, and ``result.selection`` records the
    filter.
    """


def recompute_means(result: Result) -> None:
    """Recompute ``valid_realizations``, ``mean_displacement``, ``lagrangian``
    and ``lagrangian_sky`` from ``displacement`` and ``valid``, as the
    reconstruction computes them: the sum over the valid realizations in
    realization order, divided by their count. Updates ``result`` in place.
    Raises if ``displacement`` or ``valid`` has the wrong size, an entry of
    ``valid`` is neither 0 nor 1, or a flag, ``tracers`` or ``tracers_sky``
    array has a size that is neither 0 nor the one of its layout."""


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
    """Outside the reconstruction's redshift cut, mask or box limits; not moved."""


class RealSpaceCatalog:
    """A catalogue moved to real space; row i describes input tracer i.

    Every array attribute is read-only; ``sky`` is in the ``angle_unit`` of
    the result corrected, converted once and kept, ``ids`` is a new copy on
    each access, and every other array is a view on memory owned by this
    object, as for Result. As Result, the type is
    ``otswap::RealSpaceCatalog`` itself: a catalogue made by another module
    is in radians, and its array getters check the sizes first."""

    @property
    def n_objects(self) -> int: ...

    @property
    def geometry(self) -> Literal["box", "lightcone"]: ...

    @property
    def ids(self) -> NDArray[np.int64]:
        """The IDs of the objects, copied from the ``Result`` corrected,
        shape (n_objects,), read-only; 0..n_objects-1 when it had none. Each
        access returns a new copy, so store it once (``ids = catalog.ids``)
        before indexing it in a loop."""

    @property
    def angle_unit(self) -> Optional[AngleUnit]:
        """The unit of ``sky``, that of the result corrected; ``"rad"`` for a
        catalogue made by another module. None for a box."""

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
        """Corrected Cartesian position, in coordinate units, shape
        (n_objects, 3). Box: the
        tracer with its axis coordinate moved. Lightcone: ``r (1 + shift /
        |r|)``, the tracer moved along its line of sight. NaN rows for the
        tracers not moved."""

    @property
    def shift(self) -> NDArray[np.float64]:
        """Shift applied along the line of sight, in coordinate units, shape
        (n_objects,):
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
        reconstruction by its redshift cut, mask or box limits
        (``Result.outside_*`` say which). See ``CorrectionStatus``."""

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
        """The width of the average, in coordinate units, as given."""

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
    give NaN. Raises if the shapes differ, if a position is not finite or is
    at the origin, or if a displacement is infinite; the message names the
    object."""


def axis_projection(displacement: ArrayLike, axis: int) -> NDArray[np.float64]:
    """Component of each displacement, shape (N, 3), along the Cartesian axis
    0, 1 or 2, the line of sight of a box; shape (N,). NaN displacements give
    NaN. Raises if the axis is not 0, 1 or 2, or if a displacement is
    infinite."""


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
        non-negative. For coordinates in Mpc/h, 10 Mpc/h is a reasonable
        starting value; the best value depends on the sample, and should be
        checked in each analysis.

    Returns
    -------
    NeighbourAverage
        The averages and, per object, the number of valid objects averaged
        and the sum of their valid realizations.

    Raises
    ------
    Error
        If the shapes differ, if a position is not finite, if the value of a
        valid object is not finite, or if sigma is negative or not finite;
        the message names the object.
    Error
        From the neighbour search, with a message starting
        ``"meshsearch: "``, if sigma > 0 and the positions span no volume (a
        single object, or all in one plane or on one line).

    Notes
    -----
    The result does not depend on the internal neighbour grid within
    rounding, and is the same bit for bit for the same inputs, whatever the
    number of threads.
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
    if the table has no growth rate, if a redshift is not finite or lies
    outside the distance table, or if the extrapolated bias at a redshift is
    not positive; the message names the object."""


def rsd_factor_box(redshift: float, distances: DistanceTable, *, bias: float) -> float:
    """The factor f / (b + 3 f / 5) of a box at a single redshift, with a
    constant, positive bias. Raises if the redshift is not finite or lies
    outside the distance table, if the table has no growth rate, or if the
    bias is not positive."""


def shift_radially(positions: ArrayLike, shift: ArrayLike) -> NDArray[np.float64]:
    """Move each position, shape (N, 3), by its shift, shape (N,), along its
    own line of sight: ``r (1 + shift / |r|)``. A positive shift moves away
    from the observer. A NaN shift gives a NaN row. Raises if the shapes
    differ, if a position is not finite or is at the origin, or if a shift is
    infinite; the message names the object."""


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
    correction: CorrectionConfig = CorrectionConfig(),
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
        Width of the average, in coordinate units; 0 for none. For coordinates
        in Mpc/h, 10 Mpc/h is a reasonable starting value; the best value
        depends on the sample, and should be checked in each analysis.
    correction : CorrectionConfig, optional
        ``weight_by_realizations``: weight each neighbour by its number of
        valid realizations as well. ``verbosity``: what the call prints; see
        ``Verbosity``. ``"normal"`` prints the line
        of the call; ``"detailed"`` before it how many of the corrected moved
        with their neighbours' average, having no valid realization, and how
        many of the uncorrected were left out of the reconstruction or had no
        valid tracer within 3 sigma.

    Raises
    ------
    Error
        For a box result, or one that is malformed or lacks its ``tracers`` or
        ``tracers_sky``; if its ``ids`` are neither empty nor one per object
        or hold a value twice; if a tracer kept is not finite or is at the
        origin; if the distance table has no growth rate; if ``bias``
        extrapolates to a bias that is not positive; or if a corrected
        comoving distance is not positive or lies outside the distance table.
        The message names the tracer, and in the last case suggests a table
        covering a wider redshift range.
    Error
        From the neighbour search, as in ``neighbour_average``, for the
        tracers corrected.

    Notes
    -----
    Tracers outside the reconstruction's redshift cut or mask are not
    corrected, take no part in any average, and have status 3. The catalogue
    carries the result's ``ids``.

    On a result read by ``otswap.io.read_displacement_complete`` the
    correction is that of the reconstruction itself, bit for bit, but for the
    right ascension and declination of ``sky``, which match within the
    rounding of the file's conversion to degrees and back.
    """


def real_space_box(
    result: Result,
    *,
    axis: int,
    redshift: float,
    distances: DistanceTable,
    bias: float,
    sigma: float,
    correction: CorrectionConfig = CorrectionConfig(),
) -> RealSpaceCatalog:
    """Move a box catalogue, the tracers of ``result``, from redshift space to
    real space, with the line of sight along ``axis``: the chain of
    ``real_space_lightcone`` with the projection on the axis, a single factor
    at ``redshift`` with the constant ``bias``, and the shift along the axis,
    without periodic wrapping.

    The box is not periodic: nothing flows through its faces, so modes on the
    scale of the box itself cannot be reconstructed. Tracers outside the
    reconstruction's box limits are not corrected, take no part in any
    average, and have status 3. The catalogue carries the result's ``ids``.

    Raises
    ------
    Error
        For a lightcone result, or one that is malformed or lacks its
        ``tracers``; if its ``ids`` are neither empty nor one per object or
        hold a value twice; if the axis is not 0, 1 or 2; or as
        ``rsd_factor_box``.
    Error
        From the neighbour search, as in ``neighbour_average``, for the
        tracers corrected.
    """
