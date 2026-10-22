# File formats

otswap reads and writes tables in two formats, chosen by the file's extension:
`.fits`, `.fit` and `.fits.gz` are FITS, anything else is ASCII. The C++
functions are in the namespace `otswap::io`, the Python ones in the module
`otswap.io`; both run the same code, so a file written from Python is the same,
byte for byte, as the one the C++ library writes from the same values. The
names of the files the `otswap` executable writes are on the
[executable](executable.md) page.

## Reading

{cpp:func}`otswap::io::read` and {py:func}`otswap.io.read` read selected
columns of a table, as double, and others exactly as 64-bit integers.

- **FITS:** the table in HDU 2; columns by name, case insensitively, scalar;
  an integer column must have an integer type.
- **ASCII:** a column given by digits alone is a 0-based index; any other
  name is looked up in the last line starting with `###` before the first data
  row, where otswap writes the names. Lines that are empty or start with the
  comment character (`#` by default) are skipped; fields are separated by the
  delimiter, blanks or tabs. An integer field is an optional `-` followed by
  digits.
- NaN and infinities are returned as stored, in both formats, so a file
  written with NaN rows reads back with them.

The other inputs:

- **A mask:** a full-sky HEALPix map, from a FITS file ({cpp:class}`otswap::Mask`,
  {py:class}`otswap.Mask`), in RING or NESTED order, or from an array in memory
  ({py:meth}`otswap.Mask.from_array`). A pixel is observed when its value is
  greater than 0.
- **A bias table:** b(z) at redshift nodes, read by
  {cpp:func}`otswap::io::readBiasTable` and {py:func}`otswap.io.read_bias_table`:
  FITS from the columns `REDSHIFT` and `BIAS`, ASCII from the first two
  columns.
- **A distance table:** computed from a flat cosmology, or built from the
  caller's own columns ({py:meth}`otswap.DistanceTable.from_table`), for
  example read with `read`.

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

The complete displacement file has its reader, which rebuilds the `Result`:

| C++ | Python | Reads |
|---|---|---|
| `io::readDisplacementComplete` | `io.read_displacement_complete` | a complete displacement file, of this version or of otswap 0.1 |

In Python the groups are strings, the C++ names in snake case
(`DisplacementGroup::TracerSky` is `"tracer_sky"`).

### Displacements

```{doxygenfunction} otswap::io::writeDisplacementMean
```

```{doxygenenum} otswap::io::DisplacementGroup
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
  numbers, finite, with an optional exponent (`2.194e-9`, `.5`, `1.`), always
  with a decimal point whatever the locale; single characters, `\t` for a tab;
  and lists of each.
- Windows line endings and a UTF-8 byte-order mark at the start of the file
  are accepted.

Every error names the file and, where there is one, the line and the key, as
`run.ini:7: key 'n_realizations': '8.0' is not an integer`. A syntax error
anywhere in the file is reported when the file is read, before any parameter
is; a value of the wrong type when that parameter is read. A program can
also list the keys it never read, with their lines (`unread()`), to refuse
unknown ones.
