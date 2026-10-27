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

__all__ = ["Column", "Table", "read", "read_bias_table", "read_displacement_complete", "read_keywords",
           "write", "write_displacement_complete", "write_displacement_mean", "write_mps_profile",
           "write_real_space_catalog"]
