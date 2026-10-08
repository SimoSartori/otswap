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
descriptions in TCOMMn, units in TUNITn, keywords in the header of HDU 2.

The writers of a ``Result``, a ``RealSpaceCatalog`` and an ``MpsProfile``
write one row per input object, in input order, tracers left out of the
reconstruction included, with NaN where nothing was computed, so a row's
index is the object's index in the input. Angles are in degrees, right
ascension in [0, 360). Each column has a description and a unit, and the
header carries keywords: PRODUCT, what the file holds; OTSWAPV, the library
version; GEOMETRY, box or lightcone; NOBJECTS; and those of the product,
listed with each writer. They are the C++ library's writers, so the C++ and
Python files are the same, byte for byte.
"""

import os
from typing import Literal, Optional, Sequence, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray
from typing_extensions import TypeAlias

from . import BiasTable, MpsProfile, RealSpaceCatalog, Result

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

    NaN and infinities are returned as stored, in both formats. Raises if a
    column is absent (the message lists the names found) or requested twice,
    a row is too short, an integer field is not an integer or lies outside
    the 64-bit range, or a FITS integer column holds a null value.
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
        characters. In FITS a value that reads entirely as an integer or a
        finite number is written as one, anything else as a string.

    Raises
    ------
    Error
        If the columns differ in length, or a ``"J"`` value is not an integer
        in [-2**31, 2**31); the message names the column and the row.
    """


def read_bias_table(file: Path, *, delimiter: str = " ", comment: str = "#") -> BiasTable:
    """Read a b(z) table: FITS from the columns named REDSHIFT and BIAS, ASCII
    from columns 0 and 1. Checked as ``BiasTable`` checks its arrays; the
    message names the file and the row."""


DisplacementGroup: TypeAlias = Literal["index", "tracer_sky", "lagrangian_sky", "tracer",
                                       "lagrangian", "displacement", "valid_realizations",
                                       "selection"]
CatalogGroup: TypeAlias = Literal["index", "sky", "cartesian", "valid_realizations",
                                  "neighbours", "shift", "status"]


def write_displacements(
    file: Path,
    result: Result,
    *,
    groups: Optional[Sequence[DisplacementGroup]] = None,
) -> None:
    """Write the mean displacement field of a result, one row per tracer.

    The groups and their columns:

    - ``"index"``: index (64-bit integer), the row of the tracer in the input;
    - ``"tracer_sky"``: tracRA, tracDec (degrees), tracRed, the observed
      tracer; lightcone only;
    - ``"lagrangian_sky"``: lagrRA, lagrDec, lagrRed, ``Result.lagrangian_sky``;
      lightcone only;
    - ``"tracer"``: tracX, tracY, tracZ (Mpc/h), ``Result.tracers``;
    - ``"lagrangian"``: lagrX, lagrY, lagrZ, ``Result.lagrangian``;
    - ``"displacement"``: displX, displY, displZ, ``Result.mean_displacement``;
    - ``"valid_realizations"``: nValidRec;
    - ``"selection"``: outsideRedshiftCut, outsideMask (0 or 1); lightcone only.

    Parameters
    ----------
    groups : sequence of DisplacementGroup, optional
        In output order; None for the default, CosmoBolognaLib's columns
        first: box tracer, lagrangian, displacement, valid_realizations;
        lightcone tracer_sky, lagrangian_sky, tracer, lagrangian,
        displacement, valid_realizations, selection.

    Notes
    -----
    Floating values are written with 9 significant digits in ASCII. Keywords:
    NREC, SEED (the seed used), CONVERG; MPS for a box; ZCUTMIN and ZCUTMAX
    for the finite bounds of a redshift cut, MASK = 1 when a mask selected
    the objects, FILTNSID and MAXPIXCR when a mask filter was applied.
    Raises for a repeated group, or one that does not apply to the geometry.
    """


def write_displacement_field(file: Path, result: Result) -> None:
    """Write every realization of a result, losslessly: one row per
    realization and tracer, ordered [realization][object].

    Columns: realization, index, tracX..tracZ (``Result.tracers``, repeated),
    lagrX..lagrZ (``Result.matched_random``), displX..displZ
    (``Result.displacement``), valid, outsideRedshiftCut, outsideMask.
    Keywords as ``write_displacements``.

    Floating values are written with 17 significant digits in ASCII, and as
    they are in FITS, so every value reads back identical, NaN included, in
    both formats; ``Result.from_arrays`` on the arrays read back gives the
    reconstruction's means bit for bit. displX..displZ equal lagr - trac bit
    for bit.
    """


def write_real_space_catalog(
    file: Path,
    catalog: RealSpaceCatalog,
    *,
    groups: Optional[Sequence[CatalogGroup]] = None,
) -> None:
    """Write a catalogue moved to real space, one row per tracer.

    The groups and their columns:

    - ``"index"``: index (64-bit integer);
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
        status; box the same without sky.

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
    anywhere."""
