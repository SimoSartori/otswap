"""Reading and writing ASCII and FITS tables.

The format of a file follows its extension: ``.fits``, ``.fit`` and
``.fits.gz`` are FITS, anything else is ASCII. The reading and writing are
the C++ library's, so a file written here is the same, byte for byte, as the
one the C++ library writes from the same values.

ASCII layout::

    ## PRODUCT = displacements / what the file holds
    # tracX   comoving position of the tracer [Mpc/h]
    # ...
    #
    ###   tracX   tracY   tracZ ...
    269.123 ...

one ``## NAME = value / comment`` line per keyword; one ``# name
description [unit]`` line per column with a description or a unit; a line
holding ``#`` alone after them, when there are any; the ``###`` line with the
column names; then one line per row, values separated by single spaces:
``"D"`` with ``precision`` significant digits, NaN as ``nan``, ``"J"`` and
``"K"`` as integers. Every header line starts with ``#``, so ``read``,
``numpy.loadtxt`` and most table readers skip them.

FITS layout: the table in HDU 2, a binary table of 1D, 1J and 1K columns;
descriptions in TCOMMn, units in TUNITn, keywords in the header of HDU 2. The
displacement writers name that table DISPLACEMENTS (EXTNAME); the complete
displacement file adds one HDU per further table, each named by its EXTNAME,
which an ASCII file writes as twin files of its own.

The writers of a ``Result`` and of a ``RealSpaceCatalog`` write one row per
input object, in input order, tracers left out of the reconstruction
included, with NaN where nothing was computed, so row i of a file describes
object i of the input. The id column, where written, holds the object's ID
(``Result.ids``, ``RealSpaceCatalog.ids``), or i when no IDs were set. Angles
are in degrees, right ascension in [0, 360) but in the complete displacement
file, which keeps it as given. Each column has a description and, where it
has one, a unit. Every header carries PRODUCT, what the file holds, and
OTSWAPV, the library version; those of a ``Result`` or a
``RealSpaceCatalog`` also carry GEOMETRY, box or lightcone, and NOBJECTS;
then those of the product, listed with each writer.
"""

import os
from typing import Literal, Optional, Sequence, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray
from typing_extensions import TypeAlias

from . import AngleUnit, BiasTable, MpsProfile, RealSpaceCatalog, Result

Path = Union[str, os.PathLike]


class Column:
    """One column of a table to write. Read-only."""

    def __init__(
        self,
        name: str,
        data: ArrayLike,
        *,
        type: Literal["D", "J", "K"] = "D",
        description: str = "",
        unit: str = "",
    ) -> None:
        """``data`` is one-dimensional. ``"D"`` (float64) and ``"J"`` (32-bit
        integer) columns convert it to float64; ``"K"`` (64-bit integer)
        columns to int64 exactly, from integer or boolean values, or from
        floating ones when each is an integer in range. The values are copied
        once, 8 bytes each. ``description`` is the ASCII header comment and
        the FITS TCOMMn, ``unit`` the ASCII ``[unit]`` and the FITS TUNITn.
        A ``"J"`` value that is not an integer in [-2**31, 2**31) raises when
        the table is written."""

    @property
    def name(self) -> str: ...

    @property
    def type(self) -> Literal["D", "J", "K"]: ...

    @property
    def description(self) -> str: ...

    @property
    def unit(self) -> str: ...

    @property
    def data(self) -> Union[NDArray[np.float64], NDArray[np.int64]]:
        """The values: int64 for a ``"K"`` column, float64 otherwise."""


class Table:
    """A table read from file. Read-only; its arrays are views that keep it
    alive."""

    @property
    def n_rows(self) -> int: ...

    @property
    def n_columns(self) -> int:
        """Columns read as float64."""

    @property
    def n_integer_columns(self) -> int:
        """Columns read as int64."""

    @property
    def values(self) -> NDArray[np.float64]:
        """The ``columns`` of ``read``, shape (n_rows, n_columns)."""

    @property
    def integers(self) -> NDArray[np.int64]:
        """The ``integer_columns`` of ``read``, shape
        (n_rows, n_integer_columns)."""


def read(
    file: Path,
    columns: Sequence[Union[str, int]],
    *,
    integer_columns: Sequence[Union[str, int]] = (),
    delimiter: str = " ",
    comment: str = "#",
) -> Table:
    """Read selected columns of a table: ``columns`` as float64,
    ``integer_columns`` exactly as int64, each in the order given.

    In FITS the columns are named, case insensitively, and must be scalar;
    an integer column must have an integer type. In ASCII an int, or a str
    of digits, is a 0-based index; any other str is a name, looked up in the
    last line starting with ``###`` before the first data row, where
    ``write`` puts the names. Lines that are empty or start with ``comment``
    are skipped; fields are separated by ``delimiter``, blanks or tabs. An
    ASCII integer field is an optional ``-`` followed by digits.

    NaN and infinities are returned as stored, in both formats. Raises if the
    file is missing or unreadable, a column is absent (the message lists the
    names found) or requested twice, a requested index is too large, a row is
    too short, an integer field is not an integer, lies outside the 64-bit
    range or comes from a FITS column of floating type, or a FITS integer
    column holds a null value; for an integer field the message names the
    row and the column.
    """


def write(
    file: Path,
    columns: Sequence[Column],
    *,
    precision: int = 9,
    keywords: Optional[Sequence[tuple[str, str, str]]] = None,
) -> None:
    """Write a table. An existing file is overwritten.

    Parameters
    ----------
    precision : int, optional
        Significant digits of a ``"D"`` value in ASCII, from 1 to 17; 17 reads
        back as the same float64. FITS stores the values themselves.
    keywords : sequence of (str, str, str), optional
        ``(name, value, comment)`` tuples, written in this order; a name is 1
        to 8 characters of A-Z, 0-9, ``-`` and ``_``, a value at most 68
        characters, and neither value nor comment holds a line break. In FITS
        a value that reads entirely as an integer or a finite number is
        written as one, anything else as a string.

    Raises
    ------
    Error
        If no column is given, the columns differ in length, or a ``"J"``
        value is not an integer in [-2**31, 2**31) (the message names the
        column and the row); if a keyword name is invalid, a keyword value
        exceeds 68 characters, or a keyword holds a line break; if the
        precision is outside [1, 17]; or if the file cannot be written.
    """


def read_bias_table(file: Path, *, delimiter: str = " ", comment: str = "#") -> BiasTable:
    """Read a b(z) table: FITS from the columns named REDSHIFT and BIAS, ASCII
    from columns 0 and 1. Checked as ``BiasTable`` checks its arrays. Raises
    if the file cannot be read (as ``read``), if it holds fewer than two
    rows, if a redshift is not finite or not greater than the one before, or
    if a bias is not finite or not positive; the message names the file and
    the row."""


DisplacementGroup: TypeAlias = Literal["id", "tracer_sky", "lagrangian_sky", "tracer",
                                       "lagrangian", "displacement", "valid_realizations",
                                       "selection", "outside_box"]
CatalogGroup: TypeAlias = Literal["id", "sky", "cartesian", "valid_realizations",
                                  "neighbours", "shift", "status"]


def write_displacement_mean(
    file: Path,
    result: Result,
    *,
    groups: Optional[Sequence[DisplacementGroup]] = None,
) -> None:
    """Write the mean displacement field of a result, one row per tracer.

    The groups and their columns:

    - ``"id"``: id (64-bit integer), ``Result.ids``: the row of the tracer in
      the input when no IDs were set;
    - ``"tracer_sky"``: tracRA, tracDec (degrees), tracRed, the observed
      tracer; lightcone only;
    - ``"lagrangian_sky"``: lagrRA, lagrDec, lagrRed, ``Result.lagrangian_sky``;
      lightcone only;
    - ``"tracer"``: tracX, tracY, tracZ (Mpc/h), ``Result.tracers``;
    - ``"lagrangian"``: lagrX, lagrY, lagrZ, ``Result.lagrangian``;
    - ``"displacement"``: displX, displY, displZ, ``Result.mean_displacement``;
    - ``"valid_realizations"``: nValidRec;
    - ``"selection"``: outsideRedshiftCut, outsideMask (0 or 1); lightcone only;
    - ``"outside_box"``: outsideBox (0 or 1), ``Result.outside_box``; box only.

    Parameters
    ----------
    groups : sequence of DisplacementGroup, optional
        In output order; None for the default, CosmoBolognaLib's columns
        first: box tracer, lagrangian, displacement, valid_realizations, and
        outside_box when a ``box`` was given to ``reconstruct_box``;
        lightcone tracer_sky, lagrangian_sky, tracer, lagrangian,
        displacement, valid_realizations, selection; id first when IDs were
        set.

    Raises
    ------
    Error
        If a group is repeated or does not apply to the geometry
        (``"outside_box"`` to a lightcone, or a lightcone group to a box), or
        if a field a group reads does not have the size of its layout;
        whatever the groups, if ``ids`` is neither empty nor one per object
        or holds a value twice, if the result's realization count differs
        from its configuration's, or if a box result's box differs from the
        limits recorded in its selection; or as ``write``.

    Notes
    -----
    Floating values are written with 9 significant digits in ASCII. Keywords:
    those of the reconstruction, listed with ``write_displacement_complete``.
    In FITS the table is named DISPLACEMENTS.
    """


def write_displacement_complete(
    file: Path,
    result: Result,
    *,
    lightcone_tables: bool = True,
) -> None:
    """Write every realization of a result, losslessly, with what a reader
    needs to rebuild it: one row per realization and tracer, ordered
    [realization][object]; and, for a lightcone, its mps(z) profile and its
    distance table.

    Columns: realization; id (``Result.ids`` as in
    ``write_displacement_mean``); for a lightcone, tracRA, tracDec (degrees,
    the right ascension as given, not folded into [0, 360)) and tracRed;
    tracX..tracZ (``Result.tracers``); displX..displZ
    (``Result.displacement``); valid; outsideRedshiftCut, outsideMask,
    outsideBox (0 or 1, in both geometries). The columns of the tracer, id and
    flags are repeated in every realization.

    Parameters
    ----------
    lightcone_tables : bool, optional
        False to leave out the profile and the table of a lightcone result;
        no effect on a box result, which has neither.

    Raises
    ------
    Error
        If a field read does not have the size of its layout
        (``tracers_sky`` included, for a lightcone), if ``ids`` is neither
        empty nor one per object or holds a value twice, if the result's
        realization count differs from its configuration's, if a box
        result's box differs from the limits recorded in its selection, if
        the profile written does not have one mps and one count per redshift,
        or as ``write``.

    Notes
    -----
    Floating values are written with 17 significant digits in ASCII, and as
    they are in FITS, so every value reads back identical, NaN included, in
    both formats; ``Result.from_arrays`` on the arrays read back gives the
    reconstruction's means bit for bit. ``read_displacement_complete`` reads
    the file back into a ``Result``.

    Keywords, those of ``write_displacement_mean`` too: NREC, SEED (the seed
    used), CONVERG, CELLSIZE, REJCROSS (0 or 1), CFMAXPIX, VERBOSE (silent,
    normal or detailed): the configuration; MPS for a box; BOXXLO, BOXXHI,
    BOXYLO, BOXYHI, BOXZLO, BOXZHI, ``Result.box`` when it is set, and
    BOXLIM = 1 when it holds the limits given to ``reconstruct_box``; ZCUTMIN
    and ZCUTMAX for the finite bounds of a redshift cut, MASK = 1 when a mask
    selected the objects, FILTNSID and MAXPIXCR when a mask filter was
    applied; NTRAC, NTRCUT, NTRMASK, NTRBOTH, NTRBOX, NRAND, NRNCUT, NRNMASK,
    NRNBOTH, NRNBOX, NDISPL and NDISPCR, the selection counts, each when it
    is not 0. Then, in this file only, MPSPART = 1 and DISTPART = 1 when the
    profile and the table are written.

    The profile (redshift, MPS, nTracers, with the keywords ZMIN, ZMAX and
    MPSREPR of ``write_mps_profile``) and the table (redshift, distance and,
    when the table has one, growthRate, with the keywords ZMIN and ZMAX) are
    written with 17 digits in ASCII, each with its own PRODUCT, complete
    displacements mps profile and complete displacements distance table. In
    FITS they are the HDUs MPSPROF and DISTANCE, after the main table
    DISPLACEMENTS. In ASCII they are twin files next to the file, the name
    with .mps and .distances inserted before the text after its last '.'
    (run.dat gives run.mps.dat and run.distances.dat), or appended when the
    name has no '.'. A file renamed needs its twins renamed by the same rule.
    Existing files of those names are overwritten: the twin run.mps.dat
    replaces a ``write_mps_profile`` file of that name, and is replaced by it
    when that is written later. Each is written when ``lightcone_tables`` is
    true and the result is a lightcone that holds it; a result made by
    ``Result.from_arrays`` holds a table but no profile.
    """


def read_displacement_complete(file: Path, *, angle_unit: AngleUnit) -> Result:
    """Read a file written by ``write_displacement_complete`` back into a
    ``Result``, its sky arrays in ``angle_unit``.

    Read as written: the geometry, ``ids``, the configuration, ``tracers``,
    ``tracers_sky``, ``displacement``, ``valid``, the three flags, the
    selection counts, the box, the mps(z) profile and the distance table.
    Recomputed, as ``recompute_means`` does: the valid
    realizations, ``mean_displacement``, ``lagrangian`` and
    ``lagrangian_sky``. ``ids`` reads as having no IDs when the file's IDs
    are 0..n_objects-1. The columns of the tracer, id and flags are taken
    from realization 0.

    Parameters
    ----------
    file : str or os.PathLike
        A complete displacement file; its extension sets the format.
    angle_unit : {"deg", "rad"}
        The unit of the sky arrays of the result; required, and without
        effect for a box file.

    Returns
    -------
    Result
        The result written, but for the differences below.

    Raises
    ------
    Error
        If the file cannot be read; if its PRODUCT is not that of
        ``write_displacement_complete``; if a keyword or a column its layout
        requires is missing or malformed; if the rows are not NOBJECTS x
        NREC, ordered in blocks of realization 0, 1, and so on; if a flag or
        valid entry is neither 0 nor 1, or an ID is repeated; if the file
        records the tracer counts of the selection (NTRAC) and one disagrees
        with the flags; if a box file declares a profile or a table; if a
        declared profile or table is missing, the message naming the file and
        the part, or is not the part written, by its PRODUCT, or is
        malformed. Every message names the file.

    Notes
    -----
    Every value reads back identical but ``tracers_sky``, whose right
    ascension and declination are equal within the rounding of the
    conversion to degrees and back; NaN reads back as NaN, with its sign and
    payload in FITS only. ``real_space_lightcone`` on the result read gives
    the catalogue of the reconstruction bit for bit, but for its right
    ascension and declination, equal within that rounding.

    The profile and the table are read when the file declares them (MPSPART,
    DISTPART): in FITS from the HDUs MPSPROF and DISTANCE, in ASCII from the
    twin files named by the rule of ``write_displacement_complete``, from the
    name of the file read. They are left empty when it does not: for a file
    written with ``lightcone_tables=False``, and for a file of otswap 0.1,
    whose result has no ``tracers_sky`` either, and whose index column is
    read as the IDs. A box file of otswap 0.1 gets the bounding box of its
    tracers as its box, the box its reconstruction used. The columns
    lagrX..lagrZ of a file of otswap 0.1 are not read.
    """


def write_real_space_catalog(
    file: Path,
    catalog: RealSpaceCatalog,
    *,
    groups: Optional[Sequence[CatalogGroup]] = None,
) -> None:
    """Write a catalogue moved to real space, one row per tracer.

    The groups and their columns:

    - ``"id"``: id (64-bit integer), ``RealSpaceCatalog.ids``;
    - ``"sky"``: tracRA, tracDec (degrees), tracRed, ``RealSpaceCatalog.sky``;
      lightcone only;
    - ``"cartesian"``: tracX, tracY, tracZ (Mpc/h), ``RealSpaceCatalog.cartesian``;
    - ``"valid_realizations"``: nValidRec;
    - ``"neighbours"``: nNeighbours, nRealizationsAveraged;
    - ``"shift"``: shift (Mpc/h), rsdFactor;
    - ``"status"``: status, as ``RealSpaceCatalog.status``.

    Parameters
    ----------
    groups : sequence of CatalogGroup, optional
        In output order; None for the default, CosmoBolognaLib's columns
        first: lightcone sky, cartesian, valid_realizations, neighbours,
        status; box the same without sky; id first when IDs were set.

    Raises
    ------
    Error
        If a group is repeated, if ``"sky"`` is given for a box, if a field a
        group reads does not have the size of its layout, if ``ids`` is
        neither empty nor one per object or holds a value twice, or as
        ``write``.

    Notes
    -----
    Floating values are written with 9 significant digits in ASCII.
    Keywords: SIGMA, WEIGHTED; NEXTRAP for a lightcone; AXIS, ZBOX and BIAS
    for a box.
    """


def write_mps_profile(file: Path, profile: MpsProfile) -> None:
    """Write a mean particle separation profile: redshift (the bin centre),
    MPS (Mpc/h), nTracers, one row per bin. Keywords: ZMIN, ZMAX, MPSREPR
    (the representative value). ``MpsProfile.at`` gives the profile
    anywhere. Raises if the profile is empty or malformed, or as ``write``."""
