# File formats

otswap reads and writes tables in two formats, chosen by the file's extension,
whatever its case: `.fits`, `.fit` and `.fits.gz` are FITS, anything else is
ASCII. The C++ functions are in the namespace `otswap::io`, the Python ones in
the module `otswap.io`; both run the same code, so a file written from Python
is the same, byte for byte, as the one the C++ library writes from the same
values. The names of the files the `otswap` executable writes are on the
[executable](executable.md) page, and the syntax of its parameter files at
the [end of this page](#parameter-files).

## A table of one's own

A table is written from its columns and read back by column name; an ID
column is written and read exactly, as 64-bit integers, and the keywords come
back as strings:

```python
import numpy as np
import otswap

rng = np.random.default_rng(1)
x, y, z = rng.random((3, 5))
ids = np.arange(100, 105)

columns = [otswap.io.Column("id", ids, type="K", description="catalogue ID"),
           otswap.io.Column("x", x), otswap.io.Column("y", y), otswap.io.Column("z", z)]
otswap.io.write("points.dat", columns, keywords=[("SIDE", "1", "box side")])

table = otswap.io.read("points.dat", ["x", "y", "z"], integer_columns=["id"])
positions, read_ids = table.values, table.integers[:, 0]
keywords = {name: value for name, value, _ in otswap.io.read_keywords("points.dat")}
```

The C++ functions have the same names in camelCase ([names side by
side](naming.md#tables)). With `points.fits` the same lines write and read a
FITS table; `read_keywords` then also returns the structural cards of the
table (`XTENSION`, `TTYPE1`, ...) before the caller's.

## Reading

{cpp:func}`otswap::io::read` and {py:func}`otswap.io.read` read selected
columns of a table, as double, and others exactly as 64-bit integers.

- **FITS:** the table in HDU 2; columns by name, case insensitively, scalar;
  an integer column must have an integer type.
- **ASCII:** a column given by digits alone is a 0-based index; any other
  name is looked up in the last line starting with `###` before the first data
  row, where otswap writes the names. A header such as `# RA Dec z`, written
  by another program, gives no names: read its columns by index. Lines that are empty, hold only
  separators, start with the comment character (`#` by default) after any
  separators, or start with `###` (a names line, whatever the comment
  character) are skipped; every other line is a data row, so a header line
  must start with the comment character. Fields are separated by the delimiter, blanks or
  tabs. Each requested field, and the first field of every data row, must be
  a number as a whole: decimal, hexadecimal (`0x10`), `nan`, `inf` or
  `infinity`, in any case and with an optional sign; `NA`, `1,5` with a blank
  delimiter, `4.5abc` or `1e` is an error. A carriage return ending a line is
  ignored. An integer field is an optional `-` followed by digits.
- **Messages** name an ASCII line as the file shows it, counted from 1 with
  comments and blank lines, a field counted from 1, and a FITS row counted
  from 1; a row that `write` is given in memory is counted from 0. These
  messages state which, as "(from 1)" or "(from 0)".
- NaN and infinities are returned as stored, in both formats, so a file
  written with NaN rows reads back with them.

{cpp:func}`otswap::io::readKeywords` and {py:func}`otswap.io.read_keywords`
read the keywords of a table, and {cpp:func}`otswap::io::columnNames` and
{py:func}`otswap.io.column_names` the names of its columns, those `read`
finds: in FITS the `TTYPEn` of the table, in ASCII the fields of the last
`###` line before the first data row.

A file can hold several tables, written by {cpp:func}`otswap::io::writeParts`
or {py:func}`otswap.io.write_parts`: in FITS one extension per table, each
named by its EXTNAME; in ASCII the first table in the file and each other one
in a twin file of its own ({cpp:func}`otswap::io::twinFile`). Every table,
the first included, is read by {cpp:func}`otswap::io::readPart` or
{py:func}`otswap.io.read_part` from its EXTNAME and twin name, and its
keywords by `readKeywords` with the EXTNAME, or on the twin file that
`twinFile` names. `readPart` reads any FITS table by its EXTNAME, whoever
wrote it.

The other inputs:

- **A mask:** a full-sky HEALPix map, from a FITS file ({cpp:class}`otswap::Mask`,
  {py:class}`otswap.Mask`), in RING or NESTED order, or from an array in memory
  ({py:meth}`otswap.Mask.from_array`). A pixel is observed when its value is
  greater than 0.
- **A bias table:** b(z) at redshift nodes, read by
  {cpp:func}`otswap::io::readBiasTable` and {py:func}`otswap.io.read_bias_table`:
  FITS from the columns `REDSHIFT` and `BIAS`, ASCII from the first two
  columns; the bias is dimensionless.
- **A distance table:** computed from a flat cosmology, or built from the
  caller's own columns ({py:meth}`otswap.DistanceTable.from_table`), for
  example read with `read`. Its distances set the coordinate units of a lightcone given
  on the sky (Mpc/h from the cosmology); with Cartesian coordinates they must
  be in the unit of the coordinates.

## The two layouts

The layouts of every table otswap writes, with {cpp:func}`otswap::io::write`,
{py:func}`otswap.io.write` or one of the writers below:

```{doxygengroup} file_formats
```

## The writers

Each product has its writer, the same in both languages. The groups of
columns and the keywords are listed with each.

| C++ | Python | Writes |
|---|---|---|
| `io::writeDisplacementMean` | `io.write_displacement_mean` | the mean displacement of each tracer |
| `io::writeDisplacementComplete` | `io.write_displacement_complete` | every realization, losslessly, with what a reader needs |
| `io::writeRealSpaceCatalog` | `io.write_real_space_catalog` | a catalogue moved to real space |
| `io::writeMpsProfile` | `io.write_mps_profile` | a mean particle separation profile |

The two displacement files have their readers, which rebuild a `Result`:

| C++ | Python | Reads |
|---|---|---|
| `io::readDisplacementMean` | `io.read_displacement_mean` | a mean displacement file of otswap 0.2 or later, as a `Result` of one realization holding the mean |
| `io::readDisplacementComplete` | `io.read_displacement_complete` | a complete displacement file of otswap 0.1 or later |

The complete file is lossless in both formats; the precision of the mean file
is given under {cpp:func}`otswap::io::writeDisplacementMean` below. A mean
file of otswap 0.1 is refused; a complete file of otswap 0.1 is read, its
`index` column as the IDs. The corrected catalogue and the profile have no
reader of their own: {cpp:func}`otswap::io::read` reads their columns by
name.

Angles in the files are in degrees. The Python readers take a required
`angle_unit`, the unit of the sky arrays of the `Result` they return.

A periodic box and a lightcone given in Cartesian coordinates are written in
the same files, with what tells them apart: the keyword `PERIODIC` and the
counts of wrapped objects for the former; for the latter, the distance cut,
the column `outsideDistanceCut`, written for every Cartesian result, and a
profile in distance. A Cartesian result reconstructed without a distance table
has NaN in its redshift columns and no distance table part. Each writer below
lists what it writes.

In Python the groups are strings, the C++ names in snake case
(`DisplacementGroup::TracerSky` is `"tracer_sky"`).

### Displacements

```{doxygenfunction} otswap::io::writeDisplacementMean
```

```{doxygenenum} otswap::io::DisplacementGroup
```

```{doxygenfunction} otswap::io::readDisplacementMean
```

### Complete displacements

```{doxygenfunction} otswap::io::writeDisplacementComplete
```

```{doxygenfunction} otswap::io::readDisplacementComplete
```

### Real-space catalogue

```{doxygenfunction} otswap::io::writeRealSpaceCatalog
```

```{doxygenenum} otswap::io::CatalogGroup
```

### Mean particle separation profile

```{doxygenfunction} otswap::io::writeMpsProfile
```

## Parameter files

{cpp:class}`otswap::io::ParameterFile`, in `otswap/ParameterFile.h`, reads a
text file of parameters, one `key = value` per line. It is C++ only. For
example:

```text
# a reconstruction
n_realizations = 8
columns  = 0 1 2          # three elements
output   = "../my output/run 1.fits"
delimiter = \t
```

- **Keys** hold letters, digits and `_`, do not start with a digit, and are
  case-sensitive. Each key is set at most once in a file.
- **Comments** start with `#` and run to the end of the line, on a line of
  their own or after a value. Blank lines are skipped.
- **Values** are lists of elements separated by spaces or tabs. A scalar
  parameter needs exactly one element; a list parameter takes them all.
  `key =` sets an empty list.
- **Double quotes** group an element that holds spaces, tabs or `#`, as in
  `output` above, and do not change its meaning: `"5"` is 5 and `"true"` is
  true. `""` is one empty element. There is no escape, so an element cannot
  hold a double quote, and a quote must not be glued to other text. Single
  quotes, braces and commas are ordinary characters: `{0, 1, 2}` is the three
  elements `{0,`, `1,` and `2}`.
- Everything after the first `=` is the value, further `=` signs included.
- **Types:** strings; `true` or `false`; integers written in decimal, with an
  optional `-` and no `+`, so `1e3` and `5.0` are not integers; floating-point
  numbers, finite, with an optional exponent (`2.194e-9`, `.5`, `1.`, `5`),
  with `.` as the decimal separator whatever the locale; single characters, `\t` for a tab;
  and lists of each.
- Windows line endings and a UTF-8 byte-order mark at the start of the file
  are accepted.

Every error names the file and, where there is one, the line and the key, as
`run.ini:7: key 'n_realizations': '8.0' is not an integer`. A syntax error
anywhere in the file is reported when the file is read, before any parameter
is; a value of the wrong type when that parameter is read. A program can
also list the keys it never read, with their lines (`unread()`), to refuse
unknown ones.
