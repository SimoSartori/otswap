# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Reading and writing ASCII and FITS tables.

The format of a file follows its extension: ``.fits``, ``.fit`` and
``.fits.gz`` are FITS, anything else is ASCII. The reading and writing are
the C++ library's, so a file written here is the same, byte for byte, as the
one the C++ library writes from the same values.
"""

from . import _otswap

Column = _otswap.io.Column
Table = _otswap.io.Table
read = _otswap.io.read
write = _otswap.io.write
read_keywords = _otswap.io.read_keywords
read_bias_table = _otswap.io.read_bias_table
read_displacement_complete = _otswap.io.read_displacement_complete
write_displacement_mean = _otswap.io.write_displacement_mean
write_displacement_complete = _otswap.io.write_displacement_complete
write_real_space_catalog = _otswap.io.write_real_space_catalog
write_mps_profile = _otswap.io.write_mps_profile
read_part = _otswap.io.read_part
has_part = _otswap.io.has_part
twin_file = _otswap.io.twin_file


class TablePart:
    """One table of a file holding several, for write_parts. The precision
    and the keywords are checked by write_parts, as write checks them."""

    __slots__ = ("_columns", "_extname", "_twin", "_precision", "_keywords")

    def __init__(self, columns, *, extname="", twin="", precision=9, keywords=None):
        try:
            self._columns = tuple(columns)
        except TypeError:
            raise _otswap.Error(f"columns must be a sequence of otswap.io.Column; got {columns!r}") from None
        for c in self._columns:
            if not isinstance(c, Column):
                raise _otswap.Error(f"columns must be a sequence of otswap.io.Column; got {columns!r}")
        for name, value in (("extname", extname), ("twin", twin)):
            if not isinstance(value, str):
                raise _otswap.Error(f"{name} must be a str; got {value!r}")
        self._extname = extname
        self._twin = twin
        self._precision = precision
        try:
            self._keywords = () if keywords is None else tuple(keywords)
        except TypeError:
            raise _otswap.Error("keywords must be a sequence of (name, value, comment) tuples; "
                                f"got {keywords!r}") from None

    @property
    def columns(self):
        return self._columns

    @property
    def extname(self):
        return self._extname

    @property
    def twin(self):
        return self._twin

    @property
    def precision(self):
        return self._precision

    @property
    def keywords(self):
        return self._keywords

    def __repr__(self):
        names = ", ".join(c.name for c in self._columns)
        return f"TablePart([{names}], extname={self._extname!r}, twin={self._twin!r})"


def write_parts(file, parts):
    """Write a file holding several tables, one per TablePart."""
    try:
        given = list(parts)
    except TypeError:
        given = None
    if given is None or not all(isinstance(p, TablePart) for p in given):
        raise _otswap.Error(f"parts must be a sequence of otswap.io.TablePart; got {parts!r}")
    _otswap.io._write_parts(file, [(p.columns, p.extname, p.twin, p.precision, p.keywords) for p in given])


__all__ = ["Column", "Table", "TablePart", "has_part", "read", "read_bias_table", "read_displacement_complete",
           "read_keywords", "read_part", "twin_file", "write", "write_displacement_complete",
           "write_displacement_mean", "write_mps_profile", "write_parts", "write_real_space_catalog"]
