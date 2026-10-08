# File formats

otswap reads and writes tables in two formats, chosen by the file's extension:
`.fits`, `.fit` and `.fits.gz` are FITS, anything else is ASCII. The C++
functions are in the namespace `otswap::io`, the Python ones in the module
`otswap.io`; both run the same code, so a file written from Python is the same,
byte for byte, as the one the C++ library writes from the same values.

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
| `io::writeDisplacements` | `io.write_displacements` | the mean displacement of each tracer |
| `io::writeDisplacementField` | `io.write_displacement_field` | every realization, losslessly |
| `io::writeRealSpaceCatalog` | `io.write_real_space_catalog` | a catalogue moved to real space |
| `io::writeMpsProfile` | `io.write_mps_profile` | a mean particle separation profile |

In Python the groups are strings, the C++ names in snake case
(`DisplacementGroup::TracerSky` is `"tracer_sky"`).

### Displacements

```{doxygenfunction} otswap::io::writeDisplacements
```

```{doxygenenum} otswap::io::DisplacementGroup
```

### Displacement field

```{doxygenfunction} otswap::io::writeDisplacementField
```

### Real-space catalogue

```{doxygenfunction} otswap::io::writeRealSpaceCatalog
```

```{doxygenenum} otswap::io::CatalogGroup
```

### Mean particle separation profile

```{doxygenfunction} otswap::io::writeMpsProfile
```
