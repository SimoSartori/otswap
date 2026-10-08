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
 *  are FITS, anything else is ASCII. Both layouts are stable, and a file
 *  written here reads back through io::read in either format.
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
 *  of HDU 2.
 *
 *  The writers of a Result, a RealSpaceCatalog and an MpsProfile write one
 *  row per input object, in input order, tracers left out of the
 *  reconstruction included, with NaN where nothing was computed, so a row's
 *  index is the object's index in the input. Angles are in degrees, right
 *  ascension in [0, 360). Each column has a description and a unit, and the
 *  header carries keywords: PRODUCT, what the file holds; OTSWAPV, the
 *  library version; GEOMETRY, box or lightcone; NOBJECTS; and those of the
 *  product, listed with each writer.
 *
 *  Included by otswap/OT.h; either header may be included first.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
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

    /// A header keyword: a FITS keyword in the header of HDU 2, an ASCII
    /// "## NAME = value / comment" line. The name is at most 8 characters
    /// of A-Z, 0-9, '-' and '_'. In FITS a value that reads entirely as an
    /// integer or a finite number is written as one, anything else as a
    /// string.
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
     *  name is invalid, the precision is outside [1, 17], or the file
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

    /// The column groups of writeDisplacements.
    enum class DisplacementGroup {
      Index,              ///< index (K): the row of the tracer in the input, from 0
      TracerSky,          ///< tracRA, tracDec (deg), tracRed: the observed tracer; lightcone only
      LagrangianSky,      ///< lagrRA, lagrDec (deg), lagrRed: Result::lagrangianSky; lightcone only
      Tracer,             ///< tracX, tracY, tracZ (Mpc/h): Result::tracers
      Lagrangian,         ///< lagrX, lagrY, lagrZ (Mpc/h): Result::lagrangian
      Displacement,       ///< displX, displY, displZ (Mpc/h): Result::meanDisplacement
      ValidRealizations,  ///< nValidRec (J): Result::validRealizations
      Selection           ///< outsideRedshiftCut, outsideMask (J, 0 or 1); lightcone only
    };

    /**
     *  @brief Write the mean displacement field of a result, one row per
     *  tracer.
     *
     *  The default groups, with the columns of CosmoBolognaLib's displacement
     *  file first, in its order, and its random columns named lagr:
     *  - box: Tracer, Lagrangian, Displacement, ValidRealizations (10
     *    columns);
     *  - lightcone: TracerSky, LagrangianSky, Tracer, Lagrangian,
     *    Displacement, ValidRealizations, Selection (18 columns).
     *
     *  'D' values are written with 9 significant digits in ASCII. Keywords:
     *  NREC, SEED (the seed used), CONVERG; MPS for a box; ZCUTMIN and
     *  ZCUTMAX for the finite bounds of a redshift cut, MASK = 1 when a mask
     *  selected the objects, FILTNSID and MAXPIXCR when a mask filter was
     *  applied.
     *
     *  @param groups the groups to write, in this order; empty for the
     *  default of the result's geometry.
     *
     *  @exception Error if a group is repeated or does not apply to the
     *  geometry, if a field a group reads does not have the size of its
     *  layout, or as io::write.
     */
    void writeDisplacements (const std::string& file, const Result& result,
                             const std::vector<DisplacementGroup>& groups = {});

    /**
     *  @brief Write every realization of a result, losslessly: one row per
     *  realization and tracer, ordered [realization][object], the order of
     *  Result::displacement.
     *
     *  Columns: realization (J, from 0), index (K), tracX..tracZ
     *  (Result::tracers, repeated in every realization), lagrX..lagrZ
     *  (Result::matchedRandom), displX..displZ (Result::displacement), valid
     *  (J), outsideRedshiftCut, outsideMask (J, repeated). Keywords as
     *  writeDisplacements.
     *
     *  'D' values are written with 17 significant digits in ASCII, and as
     *  they are in FITS, so every double reads back identical, NaN rows as
     *  NaN, in both formats. A Result filled from such a file, with
     *  recomputeMeans, has the reconstruction's meanDisplacement,
     *  validRealizations and lagrangian bit for bit, and its lagrangianSky
     *  given the distance table. displacement also equals lagr - trac bit for
     *  bit, as the reconstruction computes it.
     *
     *  @exception Error if a field read does not have the size of its
     *  layout, or as io::write.
     */
    void writeDisplacementField (const std::string& file, const Result& result);

    /// The column groups of writeRealSpaceCatalog.
    enum class CatalogGroup {
      Index,              ///< index (K): the row of the tracer in the input, from 0
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
     *  - box: Cartesian, ValidRealizations, Neighbours, Status (7 columns).
     *
     *  'D' values are written with 9 significant digits in ASCII. Keywords:
     *  SIGMA, WEIGHTED (0 or 1); NEXTRAP for a lightcone; AXIS, ZBOX and BIAS
     *  for a box.
     *
     *  @param groups the groups to write, in this order; empty for the
     *  default of the catalogue's geometry.
     *
     *  @exception Error as writeDisplacements.
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
