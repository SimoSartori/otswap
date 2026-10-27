// SPDX-License-Identifier: GPL-2.0-or-later
/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * This program is free software; you can redistribute it and/or    *
 * modify it under the terms of the GNU General Public License as   *
 * published by the Free Software Foundation; either version 2 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU General Public License for more details.                     *
 *                                                                  *
 * You should have received a copy of the GNU General Public        *
 * License along with this program; if not, write to the Free       *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,      *
 * Boston, MA 02110-1301 USA.                                       *
 ********************************************************************/

/**
 *  @file include/otswap/io.h
 *
 *  @brief Reading and writing ASCII and FITS tables.
 *
 *  The format of a file follows its extension: .fits, .fit and .fits.gz
 *  are FITS, anything else is ASCII. Both layouts, described in the group
 *  file_formats below, are stable, and a file written here reads back
 *  through io::read in either format.
 *
 *  Included by otswap/OT.h; either header may be included first.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

/**
 *  @defgroup file_formats File formats
 *
 *  @brief The layouts of the tables io::write, io::writeParts and the
 *  writers produce.
 *
 *  ASCII:
 *
 *      ## PRODUCT = displacements / what the file holds
 *      # tracX   comoving position of the tracer [Mpc/h]
 *      # ...
 *      #
 *      ###   tracX   tracY   tracZ ...
 *      269.123 ...
 *
 *  - one "## NAME = value / comment" line per keyword, when keywords are
 *    given;
 *  - one "# name   description [unit]" line per column with a description
 *    or a unit, the names padded to a common width, and a line holding "#"
 *    alone after the keyword and description lines, when there are any;
 *  - the "###" line, with the column names;
 *  - one line per row, values separated by single spaces: 'D' with
 *    WriteOptions::precision significant digits, a NaN as nan whatever its
 *    sign, 'J' and 'K' as integers.
 *
 *  Every header line starts with '#', so io::read, numpy.loadtxt and most
 *  table readers skip them.
 *
 *  FITS: a new file, its table in HDU 2 as a binary table of 1D, 1J and 1K
 *  columns; descriptions in TCOMMn, units in TUNITn, keywords in the header
 *  of HDU 2. The displacement writers name that table DISPLACEMENTS
 *  (EXTNAME); the complete displacement file adds one HDU per further
 *  table, each named by its EXTNAME, which an ASCII file writes as twin
 *  files of its own. io::writeParts writes such files from the caller's
 *  tables.
 *
 *  The writers of a Result and of a RealSpaceCatalog write one row per
 *  input object, in input order, tracers left out of the reconstruction
 *  included, with NaN where nothing was computed, so row i of a file
 *  describes object i of the input. The id column, where written, holds
 *  the object's ID (Result::ids, RealSpaceCatalog::ids), or i when no IDs
 *  were set. Angles are in degrees, right ascension in [0, 360) but in the
 *  complete displacement file, which keeps it as given. Each column has a
 *  description and, where it has one, a unit. Every header carries PRODUCT,
 *  what the file holds, and OTSWAPV, the library version; those of a Result
 *  or a RealSpaceCatalog also carry GEOMETRY, box or lightcone, and
 *  NOBJECTS; then those of the product, listed with each writer.
 */

#ifndef OTSWAP_IO_H
#define OTSWAP_IO_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "otswap/OT.h"

namespace otswap {

  struct BiasTable;
  struct RealSpaceCatalog;

  namespace io {

    // ========================================================================
    // Tables
    // ========================================================================

    /**
     *  @brief One column of a table to write.
     *
     *  'D' (double) and 'J' (32-bit integer) columns hold their values in
     *  data; 'K' (64-bit integer) columns in integers, exactly. The other
     *  vector is ignored.
     */
    struct Column {
      std::string name;
      char        type = 'D';                ///< 'D', 'J' or 'K'
      std::string description;               ///< ASCII header comment; FITS TCOMMn
      std::string unit;                      ///< ASCII "[unit]" after the description; FITS TUNITn
      std::vector<double>       data;        ///< the values of a 'D' or 'J' column
      std::vector<std::int64_t> integers;    ///< the values of a 'K' column
    };

    /// A table read from file.
    struct Table {
      std::size_t nRows = 0;
      std::size_t nColumns = 0;              ///< the columns read as double
      std::vector<double> values;            ///< flat, row-major [row][column]
      std::size_t nIntegerColumns = 0;       ///< the columns read as 64-bit integers
      std::vector<std::int64_t> integers;    ///< flat, row-major [row][integer column]
    };

    /// A header keyword: a FITS keyword in the header of a table's HDU, an
    /// ASCII "## NAME = value / comment" line. The name is at most 8 characters
    /// of A-Z, 0-9, '-' and '_'. The value is at most 68 characters, and
    /// neither value nor comment holds a line break. In FITS a value that
    /// reads entirely as an integer or a finite number is written as one,
    /// anything else as a string.
    struct Keyword {
      std::string name;
      std::string value;
      std::string comment;
    };

    /// Options of io::write.
    struct WriteOptions {
      /// Significant digits of a 'D' value in ASCII, from 1 to 17; 17
      /// reads back as the same double. FITS stores the doubles themselves.
      int precision = 9;

      /// Keywords written in the header, in this order.
      std::vector<Keyword> keywords;
    };

    /**
     *  @brief Read selected columns of a table, as double.
     *
     *  In FITS the columns are named, case insensitively, and must be
     *  scalar. In ASCII a column given by digits alone is a 0-based index;
     *  any other string is a name, looked up in the last line starting with
     *  "###" before the first data row, which is where io::write puts the
     *  names. Lines that are empty or start with @p comment are skipped, and
     *  fields are separated by @p delimiter, blanks or tabs.
     *
     *  NaN and infinities are returned as stored, in both formats: a file
     *  io::write wrote with NaN rows reads back with them. A null value of
     *  a FITS integer column (TNULL) is an error.
     *
     *  @exception Error if the file is missing or unreadable, a column is
     *  absent (the message lists the names found), a column is requested
     *  twice, a requested index is too large, a row is too short, or a null
     *  integer is met.
     */
    Table read (const std::string& file,
                const std::vector<std::string>& columns,
                char delimiter = ' ', char comment = '#');

    /**
     *  @brief Read selected columns as double, and others exactly as 64-bit
     *  integers.
     *
     *  @p columns are returned in Table::values, @p integerColumns in
     *  Table::integers, each in the order requested. In ASCII an integer
     *  field is an optional '-' followed by digits; in FITS an integer
     *  column has an integer type (B, I, J or K).
     *
     *  @exception Error as above, and if an integer field holds anything
     *  else (a decimal point, an exponent, a '+', nan), lies outside the
     *  64-bit range, is null, or comes from a FITS column of floating type;
     *  the message names the row and the column.
     */
    Table read (const std::string& file,
                const std::vector<std::string>& columns,
                const std::vector<std::string>& integerColumns,
                char delimiter = ' ', char comment = '#');

    /**
     *  @brief Write a table. An existing file is overwritten.
     *
     *  @exception Error if no column is given, a column type is not 'D',
     *  'J' or 'K', the columns differ in length (data for 'D' and 'J',
     *  integers for 'K'), a 'J' value is NaN, not an integer, or outside
     *  [-2^31, 2^31) (the message names the column and the row), a keyword
     *  name is invalid, a keyword value exceeds 68 characters, a keyword
     *  holds a line break, the precision is outside [1, 17], or the file
     *  cannot be written.
     */
    void write (const std::string& file, const std::vector<Column>& columns,
                const WriteOptions& options = {});

    /**
     *  @brief Write a table whose rows are produced on demand.
     *
     *  The data and integers of each column are ignored; fillRow is called
     *  once per row, in order, with a buffer of one entry per column.
     *  Nothing is held in memory, so this is the variant for catalogue-sized
     *  output. 'K' columns need the overload below.
     *
     *  @exception Error as above, and if a column is of type 'K'.
     */
    void write (const std::string& file, const std::vector<Column>& columns,
                std::size_t nRows,
                const std::function<void(std::size_t, std::vector<double>&)>& fillRow,
                const WriteOptions& options = {});

    /**
     *  @brief As above, with 'K' columns: fillRow receives two buffers of
     *  one entry per column, and each column reads its own entry of the
     *  buffer of its type, values for 'D' and 'J', integers for 'K'; the
     *  entry of the other buffer is ignored.
     */
    void write (const std::string& file, const std::vector<Column>& columns,
                std::size_t nRows,
                const std::function<void(std::size_t, std::vector<double>&,
                                         std::vector<std::int64_t>&)>& fillRow,
                const WriteOptions& options = {});

    /**
     *  @brief Read the keywords of a table, in file order.
     *
     *  FITS (a file name ending in .fits, .fit or .fits.gz): the header of
     *  HDU 2 or, with @p extname, of the HDU of that EXTNAME. Every card is
     *  returned but COMMENT, HISTORY and those without a name, the
     *  structural ones included (XTENSION, NAXISn, TTYPEn, TFORMn, EXTNAME
     *  and the like). A string value is returned without its quotes and
     *  outer blanks, any other value as its card text: a number as written,
     *  with an exponent in E, a logical as T or F.
     *
     *  ASCII: the "## NAME = value / comment" lines of the header, before
     *  the "###" line or the first data row; @p extname is ignored. The
     *  value is the text after " = " up to the first " / ", as written, so
     *  a value that itself holds " / " is cut there.
     *
     *  Comments are not read: Keyword::comment is empty. A keyword written
     *  by io::write reads back with the same name and, but for the
     *  exponent's case in FITS and outer blanks of a FITS string, the same
     *  value; numbers written with 17 significant digits read back as the
     *  same double.
     *
     *  @exception Error if the file cannot be opened, is not a valid FITS
     *  file, has no HDU 2, or has no HDU of that EXTNAME (the message names
     *  the file).
     */
    std::vector<Keyword> readKeywords (const std::string& file, const std::string& extname = "");

    /**
     *  @brief One table of a file holding several, for writeParts.
     *
     *  In FITS each part is a binary table of the same file, the first in
     *  HDU 2 and the others in the next HDUs, in order, each with its
     *  EXTNAME. In ASCII the first part is written to the file itself and
     *  each other one to its twin file, twinFile(file, twin); each file is
     *  a table as io::write writes it, readable by io::read.
     *
     *  The rows come from fillRow when it is set, as with the row overload
     *  of io::write: nRows rows, the data and integers of the columns
     *  ignored. Without it they come from the columns, which must then all
     *  hold the same number of entries, as with io::write; nRows is
     *  ignored.
     */
    struct TablePart {
      /// The EXTNAME of the FITS table: 1 to 68 characters of A-Z, a-z,
      /// 0-9, '-' and '_'. Optional for the first part; required for the
      /// others, and unique in the file, ignoring case.
      std::string extname;

      /// The name of the ASCII twin file: 1 to 68 characters of A-Z, a-z,
      /// 0-9, '-' and '_'. Ignored for the first part; required for the
      /// others, and unique in the file, ignoring case.
      std::string twin;

      std::vector<Column> columns;
      std::size_t nRows = 0;
      std::function<void(std::size_t, std::vector<double>&, std::vector<std::int64_t>&)> fillRow;
      WriteOptions options;
    };

    /**
     *  @brief Write a file holding several tables: FITS, one binary table
     *  per part, from HDU 2; ASCII, the first part in the file and each
     *  other one in its twin file. Existing files are overwritten; twin
     *  files of an earlier write that this one does not write are left
     *  where they are, so hasPart still finds them.
     *
     *  @exception Error if no part is given, if an extname or a twin is
     *  missing, malformed or repeated as described in TablePart, if a part
     *  holds a keyword named EXTNAME, or as io::write for any part. Every
     *  part is checked before any file is opened, but for the values of
     *  'J' columns, checked as they are written.
     */
    void writeParts (const std::string& file, const std::vector<TablePart>& parts);

    /**
     *  @brief Read selected columns of one part of a file writeParts wrote,
     *  the first included, as io::read does, given the part's extname and
     *  twin. FITS: the table of that EXTNAME, compared ignoring case, or
     *  with an empty extname the first table (HDU 2); twin is ignored.
     *  ASCII: the table in twinFile(file, twin), or with an empty twin in
     *  the file itself, read with the default delimiter and comment
     *  character; extname is ignored.
     *
     *  @exception Error as io::read, and if the file has no table of that
     *  EXTNAME, or the twin file is missing.
     */
    Table readPart (const std::string& file, const std::string& extname, const std::string& twin,
                    const std::vector<std::string>& columns,
                    const std::vector<std::string>& integerColumns = {});

    /**
     *  @brief Whether a file holds a part, the first included, given the
     *  part's extname and twin. FITS: whether it holds a table of that
     *  EXTNAME, compared ignoring case, or with an empty extname a first
     *  table (HDU 2); twin is ignored. ASCII: whether twinFile(file, twin)
     *  or, with an empty twin, the file itself can be opened; extname is
     *  ignored, and with a twin the file itself is not checked.
     *
     *  @exception Error if a FITS file cannot be opened.
     */
    bool hasPart (const std::string& file, const std::string& extname, const std::string& twin);

    /**
     *  @brief The ASCII twin file of a part: \<stem\>.\<twin\>.\<ext\>, where
     *  ext is the text after the last '.' of the file's name and stem what
     *  precedes it, or \<file\>.\<twin\> when the name holds no '.'. For
     *  example out/field.dat with twin mps gives out/field.mps.dat.
     */
    std::string twinFile (const std::string& file, const std::string& twin);

    // ========================================================================
    // Bias table
    // ========================================================================

    /**
     *  @brief Read a b(z) table.
     *
     *  The format follows the extension: FITS from the columns named
     *  REDSHIFT and BIAS (case insensitively); ASCII from columns 0 and 1,
     *  with the given delimiter and comment character.
     *
     *  The table is checked as rsdFactor checks a BiasTable, so a table read
     *  without error is one rsdFactor accepts.
     *
     *  @exception Error if the file cannot be read (as io::read), if it
     *  holds fewer than two rows, if a redshift is not finite or not greater
     *  than the one before, or if a bias is not finite or not positive; the
     *  message names the file and the row.
     */
    BiasTable readBiasTable (const std::string& file,
                             char delimiter = ' ', char comment = '#');

    // ========================================================================
    // Writers
    // ========================================================================

    /// The column groups of writeDisplacementMean.
    enum class DisplacementGroup {
      Id,                 ///< id (K): Result::ids, or the row of the tracer in the input, from 0, when it is empty
      TracerSky,          ///< tracRA, tracDec (deg), tracRed: the observed tracer; lightcone only
      LagrangianSky,      ///< lagrRA, lagrDec (deg), lagrRed: Result::lagrangianSky; lightcone only
      Tracer,             ///< tracX, tracY, tracZ (Mpc/h): Result::tracers
      Lagrangian,         ///< lagrX, lagrY, lagrZ (Mpc/h): Result::lagrangian
      Displacement,       ///< displX, displY, displZ (Mpc/h): Result::meanDisplacement
      ValidRealizations,  ///< nValidRec (J): Result::validRealizations
      Selection,          ///< outsideRedshiftCut, outsideMask (J, 0 or 1); lightcone only
      OutsideBox          ///< outsideBox (J, 0 or 1): Result::outsideBox; box only
    };

    /**
     *  @brief Write the mean displacement field of a result, one row per
     *  tracer.
     *
     *  The default groups, with the columns of CosmoBolognaLib's displacement
     *  file first, in its order, and its random columns named lagr:
     *  - box: Tracer, Lagrangian, Displacement, ValidRealizations (10
     *    columns), and OutsideBox (11) when limits were given to
     *    reconstructBox;
     *  - lightcone: TracerSky, LagrangianSky, Tracer, Lagrangian,
     *    Displacement, ValidRealizations, Selection (18 columns);
     *  - Id first (one more column) when Result::ids is not empty.
     *
     *  'D' values are written with 9 significant digits in ASCII. Keywords:
     *  those of the reconstruction, listed with writeDisplacementComplete.
     *  In FITS the table is named DISPLACEMENTS.
     *
     *  @param file the file to write; its extension sets the format.
     *  @param result a box or lightcone result.
     *  @param groups the groups to write, in this order; empty for the
     *  default of the result's geometry.
     *
     *  @exception Error if a group is repeated or does not apply to the
     *  geometry (OutsideBox to a lightcone, or a lightcone group to a box),
     *  if a field a group reads does not have the size of its
     *  layout, whatever the groups: if ids is neither empty nor one per
     *  object or holds a value twice, if result.nRealizations and
     *  result.config.nRealizations differ, or if selection.boxLimits is set
     *  and box is empty or differs from it; or as io::write.
     */
    void writeDisplacementMean (const std::string& file, const Result& result,
                                const std::vector<DisplacementGroup>& groups = {});

    /**
     *  @brief Write every realization of a result, losslessly, with what a
     *  reader needs to rebuild it: one row per realization and tracer,
     *  ordered [realization][object], the order of Result::displacement;
     *  and, for a lightcone, its mps(z) profile and its distance table.
     *
     *  Columns: realization (J, from 0); id (K, the object's ID as in
     *  writeDisplacementMean); for a lightcone, tracRA, tracDec (deg) and
     *  tracRed (Result::tracersSky, the right ascension as given, not folded
     *  into [0, 360)); tracX..tracZ (Result::tracers); displX..displZ
     *  (Result::displacement); valid (J); outsideRedshiftCut, outsideMask,
     *  outsideBox (J, 0 or 1, in both geometries). The columns of the
     *  tracer, id and flags are repeated in every realization.
     *
     *  'D' values are written with 17 significant digits in ASCII, and as
     *  they are in FITS, so every double reads back identical, NaN rows as
     *  NaN, in both formats. A Result filled from such a file, with
     *  recomputeMeans, has the reconstruction's meanDisplacement,
     *  validRealizations and lagrangian bit for bit, and its lagrangianSky
     *  given the distance table. readDisplacementComplete reads the file
     *  back into a Result.
     *
     *  Keywords, those of writeDisplacementMean too: NREC, SEED (the seed
     *  used), CONVERG, CELLSIZE, REJCROSS (0 or 1), CFMAXPIX, VERBOSE
     *  (silent, normal or detailed): the configuration; MPS for a box;
     *  BOXXLO, BOXXHI, BOXYLO, BOXYHI, BOXZLO, BOXZHI, Result::box when it
     *  is set, and BOXLIM = 1 when it holds the limits given to
     *  reconstructBox; ZCUTMIN and ZCUTMAX for the finite bounds of a
     *  redshift cut, MASK = 1 when a mask selected the objects, FILTNSID and
     *  MAXPIXCR when a mask filter was applied; NTRAC, NTRCUT, NTRMASK,
     *  NTRBOTH, NTRBOX, NRAND, NRNCUT, NRNMASK, NRNBOTH, NRNBOX, NDISPL and
     *  NDISPCR, the counts of Result::selection, each when it is not 0.
     *  Doubles are written as the shortest decimal that reads back as the
     *  same value. Then, in this file only, MPSPART = 1 and DISTPART = 1
     *  when the profile and the table are written.
     *
     *  The profile, columns redshift, MPS (Mpc/h) and nTracers (J), with the
     *  keywords ZMIN, ZMAX and MPSREPR of writeMpsProfile, and the table,
     *  columns redshift, distance (Mpc/h) and, when the table has one,
     *  growthRate, with the keywords ZMIN and ZMAX, its ends, are written
     *  with 17 digits in ASCII, each with its own PRODUCT, complete
     *  displacements mps profile and complete displacements distance table.
     *  In FITS they are the HDUs MPSPROF and DISTANCE, after the main table
     *  DISPLACEMENTS. In ASCII they are twin files next to the file, the
     *  name with .mps and .distances inserted before the text after its last
     *  '.' (run.dat gives run.mps.dat and run.distances.dat), or appended
     *  when the name has no '.'. A file renamed needs its twins renamed by
     *  the same rule. Existing files of those names are overwritten: the
     *  twin run.mps.dat replaces a writeMpsProfile file of that name, and is
     *  replaced by it when that is written later. Each is written when
     *  lightconeTables is true, the result is a lightcone and holds it: the
     *  profile when mpsProfile is not empty, the table when distances is not
     *  null.
     *
     *  @param file the file to write; its extension sets the format.
     *  @param result a box or lightcone result.
     *  @param lightconeTables false to leave out the profile and the table
     *  of a lightcone result; no effect on a box result, which has neither.
     *
     *  @exception Error if a field read does not have the size of its
     *  layout (tracersSky included, for a lightcone), if ids is neither empty
     *  nor one per object or holds a value twice, if result.nRealizations
     *  and result.config.nRealizations differ, if selection.boxLimits is set
     *  and box is empty or differs from it, if the profile written does not
     *  have one mps and one count per redshift, or as io::write.
     */
    void writeDisplacementComplete (const std::string& file, const Result& result,
                                    bool lightconeTables = true);

    /**
     *  @brief Read a file written by writeDisplacementComplete back into a
     *  Result.
     *
     *  Read as written: the geometry, ids, config (but log, which is
     *  &std::clog), tracers, tracersSky (degrees to radians), displacement,
     *  valid, the three flags, filteredNside, selection, mps,
     *  box, mpsProfile, and the distance table, as a new table equal to the
     *  one written. Recomputed with recomputeMeans: validRealizations,
     *  meanDisplacement, lagrangian and lagrangianSky. elapsedSeconds is 0.
     *  ids is empty when the file's IDs are 0..nObjects-1. The columns of the
     *  tracer, id and flags are taken from realization 0.
     *
     *  Every double reads back identical, and so does every field above but
     *  tracersSky, whose right ascension and declination are equal within
     *  the rounding of the conversion to degrees and back; NaN reads back as
     *  NaN, with its sign and payload in FITS only. realSpaceLightcone on the
     *  Result read gives the catalogue of the reconstruction bit for bit,
     *  but for its right ascension and declination, equal within that
     *  rounding.
     *
     *  The profile and the table are read when the main table declares them
     *  (MPSPART, DISTPART): in FITS from the HDUs MPSPROF and DISTANCE, in
     *  ASCII from the twin files named by the rule of
     *  writeDisplacementComplete, from the name of the file read. They are
     *  left empty when it does not: for a file written with lightconeTables
     *  false, and for a file of otswap 0.1, whose Result has no tracersSky
     *  either, and whose index column is read as the IDs. A box file of
     *  otswap 0.1 gets the bounding box of its tracers as Result::box, the
     *  box its reconstruction used. The columns lagrX..lagrZ of a file of
     *  otswap 0.1 are not read.
     *
     *  @param file a complete displacement file; its extension sets the
     *  format.
     *
     *  @return the result, as writeDisplacementComplete received it but for
     *  the differences above.
     *
     *  @exception Error if the file cannot be read (as io::read); if its
     *  PRODUCT is not that of writeDisplacementComplete; if a keyword or a
     *  column its layout requires is missing or malformed; if the rows are
     *  not NOBJECTS x NREC, ordered in blocks of realization 0, 1, and so
     *  on; if a flag or valid entry is neither 0 nor 1, or an ID is
     *  repeated; if the file records the tracer counts of the selection
     *  (NTRAC) and one disagrees with the flags; if a box file declares a
     *  profile or a table; if a declared profile or table is missing, the
     *  message naming the file and the part, or is not the part written, by
     *  its PRODUCT, or is malformed. Every message names the file.
     */
    Result readDisplacementComplete (const std::string& file);

    /// The column groups of writeRealSpaceCatalog.
    enum class CatalogGroup {
      Id,                 ///< id (K): RealSpaceCatalog::ids, or the row of the tracer in the input, from 0, when it is empty
      Sky,                ///< tracRA, tracDec (deg), tracRed: RealSpaceCatalog::sky; lightcone only
      Cartesian,          ///< tracX, tracY, tracZ (Mpc/h): RealSpaceCatalog::cartesian
      ValidRealizations,  ///< nValidRec (J)
      Neighbours,         ///< nNeighbours, nRealizationsAveraged (J)
      Shift,              ///< shift (Mpc/h), rsdFactor: RealSpaceCatalog::shift, factor
      Status              ///< status (J): RealSpaceCatalog::status
    };

    /**
     *  @brief Write a catalogue moved to real space, one row per tracer.
     *
     *  The default groups, CosmoBolognaLib's reconstructed catalogue first:
     *  - lightcone: Sky, Cartesian, ValidRealizations, Neighbours, Status
     *    (10 columns);
     *  - box: Cartesian, ValidRealizations, Neighbours, Status (7 columns);
     *  - Id first (one more column) when RealSpaceCatalog::ids is not empty.
     *
     *  'D' values are written with 9 significant digits in ASCII. Keywords:
     *  SIGMA, WEIGHTED (0 or 1); NEXTRAP for a lightcone; AXIS, ZBOX and BIAS
     *  for a box.
     *
     *  @param file the file to write; its extension sets the format.
     *  @param catalog a box or lightcone catalogue.
     *  @param groups the groups to write, in this order; empty for the
     *  default of the catalogue's geometry.
     *
     *  @exception Error if a group is repeated, if Sky is given for a box,
     *  if a field a group reads does not have the size of its layout, if
     *  ids is neither empty nor one per object or holds a value twice, or as
     *  io::write.
     */
    void writeRealSpaceCatalog (const std::string& file, const RealSpaceCatalog& catalog,
                                const std::vector<CatalogGroup>& groups = {});

    /**
     *  @brief Write a mean particle separation profile: redshift (the bin
     *  centre), MPS (Mpc/h), nTracers (J), one row per bin. Keywords: ZMIN,
     *  ZMAX, MPSREPR (the representative value). MpsProfile::at gives the
     *  profile anywhere.
     *
     *  @exception Error if the profile is empty or malformed, or as
     *  io::write.
     */
    void writeMpsProfile (const std::string& file, const MpsProfile& profile);

  }

}

#endif
