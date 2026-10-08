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
 *  @file src/Table.cpp
 *
 *  @brief Reading and writing ASCII and FITS tables.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

#include <fitsio.h>

#include "internal.h"

namespace {

  using Fill = std::function<void(std::size_t, std::vector<double>&, std::vector<std::int64_t>&)>;

  // Largest 0-based column index an ASCII request may name. A larger one
  // is far more likely to be a mistyped column list than a real table.
  constexpr std::size_t kMaxAsciiColumnIndex = 4095;

  /// Most significant digits written to ASCII: 17 read back as the same
  /// double.
  constexpr int kMaxPrecision = 17;

  /// Longest FITS string keyword value, the room an 80-character card
  /// leaves.
  constexpr std::size_t kMaxKeywordValue = 68;

  bool has_extension (const std::string& file, const std::string& ext)
  {
    if (file.size() < ext.size()) return false;
    std::string tail = file.substr(file.size() - ext.size());
    std::transform(tail.begin(), tail.end(), tail.begin(),
                   [] (unsigned char c) { return (char)std::tolower(c); });
    return tail == ext;
  }

  bool is_fits (const std::string& file)
  {
    return has_extension(file, ".fits") || has_extension(file, ".fit") ||
           has_extension(file, ".fits.gz");
  }

  bool is_index (const std::string& column)
  {
    return !column.empty() &&
           std::all_of(column.begin(), column.end(),
                       [] (unsigned char c) { return std::isdigit(c) != 0; });
  }

  bool separator (const char c, const char delimiter)
  {
    return c == delimiter || c == ' ' || c == '\t';
  }

  /// True when the line holds a field strtod reads as a number, as a data
  /// row does.
  bool has_number (const std::string& line)
  {
    const char* p = line.c_str();
    char* end = nullptr;
    while (*p != '\0') {
      std::strtod(p, &end);
      if (end != p) return true;
      ++p;
    }
    return false;
  }

  /// The column names of an ASCII file: the fields of the last line
  /// starting with "###" before the first data row. found is false when
  /// there is no such line.
  std::vector<std::string> ascii_names (const std::string& file, const char delimiter,
                                        const char comment, bool& found)
  {
    std::ifstream input(file);
    if (!input)
      throw otswap::Error("cannot open the input file: " + file);

    std::vector<std::string> names;
    found = false;
    std::string line;
    while (std::getline(input, line)) {
      if (line.compare(0, 3, "###") == 0) {
        found = true;
        names.clear();
        std::size_t i = 3;
        while (i < line.size()) {
          while (i < line.size() && separator(line[i], delimiter)) ++i;
          const std::size_t start = i;
          while (i < line.size() && !separator(line[i], delimiter)) ++i;
          if (i > start) names.push_back(line.substr(start, i - start));
        }
        continue;
      }
      if (line.empty() || line[0] == comment) continue;
      if (has_number(line)) break;
    }
    return names;
  }

  otswap::io::Table read_ascii (const std::string& file,
                                const std::vector<std::string>& columns,
                                const std::vector<std::string>& integerColumns,
                                const char delimiter, const char comment)
  {
    std::vector<std::string> names;
    bool namesFound = false;
    const bool byName =
      std::any_of(columns.begin(), columns.end(), [] (const std::string& c) { return !is_index(c); }) ||
      std::any_of(integerColumns.begin(), integerColumns.end(),
                  [] (const std::string& c) { return !is_index(c); });
    if (byName) names = ascii_names(file, delimiter, comment, namesFound);

    auto resolve = [&] (const std::string& c) -> std::size_t {
      if (is_index(c)) {
        std::size_t value = 0;
        const auto conversion = std::from_chars(c.data(), c.data()+c.size(), value);
        if (conversion.ec != std::errc() || value > kMaxAsciiColumnIndex)
          throw otswap::Error("the ASCII column index " + c +
                              " is too large; the maximum supported index is " +
                              std::to_string(kMaxAsciiColumnIndex));
        return value;
      }
      if (!namesFound)
        throw otswap::Error("the column '" + c + "' is requested by name, but " + file +
                            " has no line starting with ### naming its columns; give the column "
                            "as a 0-based index");
      const auto at = std::find(names.begin(), names.end(), c);
      if (at == names.end()) {
        std::string list;
        for (const std::string& n : names) list += (list.empty() ? "" : ", ") + n;
        throw otswap::Error("the column '" + c + "' was not found in " + file +
                            "; its ### line names: " + list);
      }
      return (std::size_t)(at - names.begin());
    };

    struct Slot { bool integer = false; std::size_t at = 0; bool requested = false; };
    std::vector<std::size_t> index;
    for (const std::string& c : columns) index.push_back(resolve(c));
    for (const std::string& c : integerColumns) index.push_back(resolve(c));

    std::size_t maxIndex = 0;
    for (const std::size_t c : index) maxIndex = std::max(maxIndex, c);

    std::vector<Slot> map(maxIndex+1);
    for (std::size_t i = 0; i < index.size(); ++i) {
      const std::string& name = i < columns.size() ? columns[i] : integerColumns[i - columns.size()];
      if (map[index[i]].requested)
        throw otswap::Error("the ASCII column '" + name + "' (index " + std::to_string(index[i]) +
                            ") is requested more than once");
      map[index[i]].requested = true;
      map[index[i]].integer = i >= columns.size();
      map[index[i]].at = i < columns.size() ? i : i - columns.size();
    }

    std::ifstream input(file);
    if (!input)
      throw otswap::Error("cannot open the input file: " + file);

    otswap::io::Table table;
    table.nColumns = columns.size();
    table.nIntegerColumns = integerColumns.size();

    std::vector<double> row(columns.size());
    std::vector<std::int64_t> integerRow(integerColumns.size());
    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(input, line)) {
      ++lineNumber;
      if (line.empty() || line[0] == comment) continue;

      const char* p = line.c_str();
      char* end = nullptr;
      std::size_t column = 0, fields = 0;
      bool complete = false;

      while (*p != '\0' && column <= maxIndex) {
        const double value = std::strtod(p, &end);
        if (end == p) { ++p; continue; }
        const Slot& slot = map[column];
        if (slot.requested && !slot.integer) row[slot.at] = value;
        if (slot.requested && slot.integer) {
          const char* start = p;
          while (start < end && std::isspace((unsigned char)*start)) ++start;
          std::int64_t integer = 0;
          const auto conversion = std::from_chars(start, (const char*)end, integer);
          if (conversion.ec != std::errc() || conversion.ptr != end)
            throw otswap::Error("the integer column '" + integerColumns[slot.at] + "' of " + file +
                                " holds " + std::string(start, (const char*)end) + " at line " +
                                std::to_string(lineNumber) +
                                (conversion.ec == std::errc::result_out_of_range
                                   ? ", outside the 64-bit integer range"
                                   : ", which is not an integer"));
          integerRow[slot.at] = integer;
        }
        ++fields;
        p = end;
        while (*p != '\0' && separator(*p, delimiter)) ++p;
        if (column == maxIndex) { complete = true; break; }
        ++column;
      }

      if (fields == 0) continue;

      if (!complete)
        throw otswap::Error("row too short in " + file + " at line " +
                            std::to_string(lineNumber) + ": found " +
                            std::to_string(fields) + " column(s), but column index " +
                            std::to_string(maxIndex) + " was requested");

      table.values.insert(table.values.end(), row.begin(), row.end());
      table.integers.insert(table.integers.end(), integerRow.begin(), integerRow.end());
      ++table.nRows;
    }

    return table;
  }

  otswap::io::Table read_fits (const std::string& file,
                               const std::vector<std::string>& columns,
                               const std::vector<std::string>& integerColumns)
  {
    fitsfile* fptr = nullptr;
    int status = 0;

    auto fail = [&] (const std::string& message) {
      int closeStatus = 0;
      if (fptr != nullptr) { fits_close_file(fptr, &closeStatus); fptr = nullptr; }
      throw otswap::Error(message);
    };

    if (fits_open_file(&fptr, file.c_str(), READONLY, &status))
      fail("cannot open the input file, or it is not a valid FITS file: " + file);

    int hdutype = 0;
    if (fits_movabs_hdu(fptr, 2, &hdutype, &status))
      fail("the first extension cannot be accessed in: " + file);

    if (hdutype != BINARY_TBL && hdutype != ASCII_TBL)
      fail("the first extension of " + file + " is not a table");

    long nrows = 0;
    if (fits_get_num_rows(fptr, &nrows, &status))
      fail("cannot read the number of rows of: " + file);

    std::vector<std::string> all = columns;
    all.insert(all.end(), integerColumns.begin(), integerColumns.end());
    const std::size_t ncolumns = all.size();
    std::vector<int> colnum(ncolumns, 0);
    std::vector<bool> floating(ncolumns, false);

    for (std::size_t i = 0; i < ncolumns; ++i) {
      std::vector<char> name(all[i].begin(), all[i].end());
      name.push_back('\0');
      if (fits_get_colnum(fptr, CASEINSEN, name.data(), &colnum[i], &status))
        fail("the column '" + all[i] + "' was not found in: " + file);
      for (std::size_t j = 0; j < i; ++j)
        if (colnum[j] == colnum[i])
          fail("the column '" + all[i] + "' of " + file + " is requested more than once");

      int typecode = 0, eqtype = 0;
      long repeat = 0, width = 0;
      if (fits_get_coltype(fptr, colnum[i], &typecode, &repeat, &width, &status) ||
          fits_get_eqcoltype(fptr, colnum[i], &eqtype, &repeat, &width, &status))
        fail("cannot read the type of the column '" + all[i] + "' in: " + file);

      if (repeat != 1)
        fail("the column '" + all[i] + "' of " + file + " has repeat count " +
             std::to_string(repeat) + ": vector columns are not supported");

      floating[i] = typecode == TFLOAT || typecode == TDOUBLE;
      const bool integer = eqtype == TBYTE || eqtype == TSBYTE || eqtype == TSHORT ||
                           eqtype == TUSHORT || eqtype == TINT || eqtype == TUINT ||
                           eqtype == TLONG || eqtype == TULONG || eqtype == TLONGLONG ||
                           eqtype == TULONGLONG;
      if (i >= columns.size() && !integer)
        fail("the column '" + all[i] + "' of " + file + " is requested as an integer column, "
             "but its values are not integers");
    }

    std::vector<std::vector<double>> byColumn(columns.size(), std::vector<double>(nrows));
    std::vector<std::vector<LONGLONG>> byIntegerColumn(integerColumns.size(),
                                                       std::vector<LONGLONG>(nrows));

    long chunk = 0;
    if (fits_get_rowsize(fptr, &chunk, &status)) { status = 0; chunk = 0; }
    if (chunk < 1) chunk = 1;

    double nulval = std::numeric_limits<double>::quiet_NaN();
    LONGLONG integerNulval = std::numeric_limits<LONGLONG>::min();
    std::vector<char> nullFlags(chunk);

    for (long first = 1; first <= nrows; first += chunk) {
      const long nread = std::min(chunk, nrows-first+1);
      for (std::size_t i = 0; i < ncolumns; ++i) {
        const bool integer = i >= columns.size();
        int anynul = 0;
        const int failed = integer
          ? fits_read_col(fptr, TLONGLONG, colnum[i], first, 1, nread, &integerNulval,
                          byIntegerColumn[i - columns.size()].data()+first-1, &anynul, &status)
          : fits_read_col(fptr, TDOUBLE, colnum[i], first, 1, nread, floating[i] ? nullptr : &nulval,
                          byColumn[i].data()+first-1, &anynul, &status);
        if (failed)
          fail("cannot read the column '" + all[i] + "' of " + file +
               (status == NUM_OVERFLOW ? ": a value lies outside the 64-bit integer range" : ""));

        if (anynul) {
          std::fill(nullFlags.begin(), nullFlags.begin()+nread, 0);
          int detailed = 0;
          std::vector<double> scratch(nread);
          if (fits_read_colnull(fptr, TDOUBLE, colnum[i], first, 1, nread,
                                scratch.data(), nullFlags.data(), &detailed, &status))
            fail("cannot locate the undefined value in the column '" + all[i] + "' of: " + file);

          const auto at = std::find_if(nullFlags.begin(), nullFlags.begin()+nread,
                                       [] (const char flag) { return flag != 0; });
          const long offset = at == nullFlags.begin()+nread ? 0 : std::distance(nullFlags.begin(), at);
          fail("undefined value in the column '" + all[i] + "', row " +
               std::to_string(first+offset) + " of: " + file);
        }
      }
    }

    int closeStatus = 0;
    const int closeResult = fits_close_file(fptr, &closeStatus);
    fptr = nullptr;
    if (closeResult)
      throw otswap::Error("cannot close the input FITS file: " + file);

    otswap::io::Table table;
    table.nRows = (std::size_t)nrows;
    table.nColumns = columns.size();
    table.nIntegerColumns = integerColumns.size();
    table.values.resize(table.nRows * table.nColumns);
    table.integers.resize(table.nRows * table.nIntegerColumns);
    for (std::size_t r = 0; r < table.nRows; ++r) {
      for (std::size_t c = 0; c < table.nColumns; ++c)
        table.values[r*table.nColumns + c] = byColumn[c][r];
      for (std::size_t c = 0; c < table.nIntegerColumns; ++c)
        table.integers[r*table.nIntegerColumns + c] = (std::int64_t)byIntegerColumn[c][r];
    }

    return table;
  }

  // --------------------------------------------------------------- writing

  void check_options (const std::vector<otswap::io::Column>& columns,
                      const otswap::io::WriteOptions& options, const std::string& file)
  {
    if (columns.empty())
      throw otswap::Error("no column was given for: " + file);

    for (const auto& c : columns)
      if (c.type != 'D' && c.type != 'J' && c.type != 'K')
        throw otswap::Error("the column '" + c.name + "' has type '" + std::string(1, c.type) +
                            "'; the types are 'D', 'J' and 'K'");

    if (options.precision < 1 || options.precision > kMaxPrecision)
      throw otswap::Error("the precision is " + std::to_string(options.precision) +
                          "; it must lie in [1, " + std::to_string(kMaxPrecision) + "]");

    for (const auto& k : options.keywords) {
      const bool valid = !k.name.empty() && k.name.size() <= 8 &&
        std::all_of(k.name.begin(), k.name.end(), [] (const char c) {
          return (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-' || c == '_';
        });
      if (!valid)
        throw otswap::Error("the keyword name '" + k.name + "' is not 1 to 8 characters of "
                            "A-Z, 0-9, '-' and '_'");
      if (k.value.size() > kMaxKeywordValue)
        throw otswap::Error("the value of the keyword " + k.name + " holds " +
                            std::to_string(k.value.size()) + " characters; at most " +
                            std::to_string(kMaxKeywordValue) + " fit a FITS card");
      for (const std::string* text : {&k.value, &k.comment})
        if (text->find_first_of("\n\r") != std::string::npos)
          throw otswap::Error("the keyword " + k.name + " holds a line break");
    }
  }

  void check_integer_value (const otswap::io::Column& column, const double value,
                            const std::size_t row)
  {
    if (!(value >= -2147483648. && value < 2147483648.) || value != std::floor(value)) {
      std::ostringstream text;
      text.imbue(std::locale::classic());
      text << std::setprecision(17) << value;
      throw otswap::Error("the 'J' column '" + column.name + "' holds " + text.str() + " at row " +
                          std::to_string(row) + "; a 'J' value must be an integer in "
                          "[-2147483648, 2147483647]");
    }
  }

  void write_ascii (const std::string& file,
                    const std::vector<otswap::io::Column>& columns,
                    const std::size_t nrows, const Fill& fillRow,
                    const otswap::io::WriteOptions& options)
  {
    std::ofstream out(file.c_str());
    if (!out)
      throw otswap::Error("cannot open the output file: " + file);
    out.imbue(std::locale::classic());

    for (const auto& k : options.keywords) {
      out << "## " << k.name << " = " << k.value;
      if (!k.comment.empty()) out << " / " << k.comment;
      out << '\n';
    }

    std::size_t width = 0;
    bool anyDescription = false;
    for (const auto& c : columns) {
      if (c.description.empty() && c.unit.empty()) continue;
      anyDescription = true;
      width = std::max(width, c.name.size());
    }

    for (const auto& c : columns) {
      if (c.description.empty() && c.unit.empty()) continue;
      out << "# " << c.name << std::string(width - c.name.size() + 3, ' ') << c.description;
      if (!c.unit.empty()) out << (c.description.empty() ? "[" : " [") << c.unit << "]";
      out << '\n';
    }
    if (anyDescription || !options.keywords.empty()) out << "#\n";

    out << "###";
    for (const auto& c : columns) out << "   " << c.name;
    out << '\n';

    out << std::defaultfloat << std::setprecision(options.precision);
    std::vector<double> values(columns.size());
    std::vector<std::int64_t> integers(columns.size());
    for (std::size_t row = 0; row < nrows; ++row) {
      fillRow(row, values, integers);
      for (std::size_t c = 0; c < columns.size(); ++c) {
        if (c > 0) out << ' ';
        if (columns[c].type == 'K') {
          out << integers[c];
        }
        else if (columns[c].type == 'J') {
          check_integer_value(columns[c], values[c], row);
          out << (long long)values[c];
        }
        else if (std::isnan(values[c])) out << "nan";
        else                             out << values[c];
      }
      out << '\n';
    }

    out.close();
    if (out.fail())
      throw otswap::Error("cannot complete writing the output file: " + file);
  }

  /// A keyword value written in FITS as a number: signs, digits, a point
  /// and an exponent only, read entirely by strtod as a finite value.
  bool numeric (const std::string& value)
  {
    if (value.empty() ||
        value.find_first_not_of("0123456789+-.eE") != std::string::npos) return false;
    char* end = nullptr;
    const double v = std::strtod(value.c_str(), &end);
    return end == value.c_str() + value.size() && std::isfinite(v);
  }

  void write_fits (const std::string& file,
                   const std::vector<otswap::io::Column>& columns,
                   const std::size_t nrows, const Fill& fillRow,
                   const otswap::io::WriteOptions& options)
  {
    fitsfile* fptr = nullptr;
    int status = 0;

    auto fail = [&] (const std::string& message) {
      int closeStatus = 0;
      if (fptr != nullptr) { fits_close_file(fptr, &closeStatus); fptr = nullptr; }
      throw otswap::Error(message);
    };

    const std::string overwrite = "!" + file;
    if (fits_create_file(&fptr, overwrite.c_str(), &status))
      fail("cannot create the output file: " + file);

    const int ncolumns = (int)columns.size();
    std::vector<std::string> nameStore(ncolumns), formStore(ncolumns);
    std::vector<char*> ttype(ncolumns), tform(ncolumns);
    for (int i = 0; i < ncolumns; ++i) {
      nameStore[i] = columns[i].name;
      formStore[i] = columns[i].type == 'J' ? "1J" : columns[i].type == 'K' ? "1K" : "1D";
      ttype[i] = &nameStore[i][0];
      tform[i] = &formStore[i][0];
    }

    if (nrows > (std::size_t)std::numeric_limits<long>::max())
      fail("too many rows for a FITS table in: " + file);

    const long fitsRows = (long)nrows;
    if (fits_create_tbl(fptr, BINARY_TBL, fitsRows, ncolumns, ttype.data(), tform.data(),
                        nullptr, nullptr, &status))
      fail("cannot create the binary table in: " + file);

    auto write_string = [&] (const std::string& key, const std::string& text,
                             const std::string& comment, const std::string& what) {
      std::vector<char> value(text.begin(), text.end());
      value.push_back('\0');
      if (fits_write_key(fptr, TSTRING, key.c_str(), value.data(),
                         comment.empty() ? nullptr : comment.c_str(), &status))
        fail("cannot write " + what + " in: " + file);
    };

    for (int i = 0; i < ncolumns; ++i) {
      const std::string what = "the description of the column '" + columns[i].name + "'";
      if (!columns[i].description.empty())
        write_string("TCOMM" + std::to_string(i+1), columns[i].description, "", what);
      if (!columns[i].unit.empty())
        write_string("TUNIT" + std::to_string(i+1), columns[i].unit, "", what);
    }

    for (const auto& k : options.keywords) {
      if (numeric(k.value)) {
        std::string literal = k.value;
        std::replace(literal.begin(), literal.end(), 'e', 'E');
        char card[FLEN_CARD];
        if (fits_make_key(k.name.c_str(), &literal[0], k.comment.c_str(), card, &status) ||
            fits_write_record(fptr, card, &status))
          fail("cannot write the keyword " + k.name + " in: " + file);
      }
      else write_string(k.name, k.value, k.comment, "the keyword " + k.name);
    }

    long chunk = 0;
    if (fits_get_rowsize(fptr, &chunk, &status)) { status = 0; chunk = 0; }
    if (chunk < 1) chunk = 1;

    std::vector<std::vector<double>> buffer(columns.size(), std::vector<double>(chunk));
    std::vector<std::vector<LONGLONG>> integerBuffer(columns.size());
    for (int i = 0; i < ncolumns; ++i)
      if (columns[i].type == 'K') integerBuffer[i].resize(chunk);
    std::vector<long> integers(chunk);
    std::vector<double> values(columns.size());
    std::vector<std::int64_t> rowIntegers(columns.size());

    try {
      for (long first = 1; first <= fitsRows; first += chunk) {
        const long nwrite = std::min(chunk, fitsRows-first+1);
        for (long row = 0; row < nwrite; ++row) {
          fillRow((std::size_t)(first-1+row), values, rowIntegers);
          for (std::size_t c = 0; c < columns.size(); ++c) {
            if (columns[c].type == 'K') integerBuffer[c][row] = (LONGLONG)rowIntegers[c];
            else                        buffer[c][row] = values[c];
          }
        }
        for (int i = 0; i < ncolumns; ++i) {
          int written = 0;
          if (columns[i].type == 'J') {
            for (long row = 0; row < nwrite; ++row) {
              check_integer_value(columns[i], buffer[i][row], (std::size_t)(first-1+row));
              integers[row] = (long)buffer[i][row];
            }
            written = fits_write_col(fptr, TLONG, i+1, first, 1, nwrite, integers.data(), &status);
          }
          else if (columns[i].type == 'K')
            written = fits_write_col(fptr, TLONGLONG, i+1, first, 1, nwrite,
                                     integerBuffer[i].data(), &status);
          else
            written = fits_write_col(fptr, TDOUBLE, i+1, first, 1, nwrite, buffer[i].data(), &status);
          if (written)
            fail("cannot write the column '" + columns[i].name + "' in: " + file);
        }
      }
    }
    catch (...) {
      int cleanup = 0;
      if (fptr != nullptr) { fits_close_file(fptr, &cleanup); fptr = nullptr; }
      throw;
    }

    int closeStatus = 0;
    const int closeResult = fits_close_file(fptr, &closeStatus);
    fptr = nullptr;
    if (closeResult)
      throw otswap::Error("cannot close the output FITS file: " + file);
  }

  void write_table (const std::string& file, const std::vector<otswap::io::Column>& columns,
                    const std::size_t nRows, const Fill& fillRow,
                    const otswap::io::WriteOptions& options)
  {
    check_options(columns, options, file);
    if (is_fits(file)) write_fits(file, columns, nRows, fillRow, options);
    else               write_ascii(file, columns, nRows, fillRow, options);
  }

}


// ============================================================================


otswap::io::Table otswap::io::read (const std::string& file,
                                    const std::vector<std::string>& columns,
                                    const char delimiter, const char comment)
{
  return read(file, columns, {}, delimiter, comment);
}


// ============================================================================


otswap::io::Table otswap::io::read (const std::string& file,
                                    const std::vector<std::string>& columns,
                                    const std::vector<std::string>& integerColumns,
                                    const char delimiter, const char comment)
{
  if (columns.empty() && integerColumns.empty())
    throw Error("no column was requested from: " + file);

  if (is_fits(file)) return read_fits(file, columns, integerColumns);
  return read_ascii(file, columns, integerColumns, delimiter, comment);
}


// ============================================================================


void otswap::io::write (const std::string& file, const std::vector<Column>& columns,
                        const WriteOptions& options)
{
  check_options(columns, options, file);

  auto length = [] (const Column& c) { return c.type == 'K' ? c.integers.size() : c.data.size(); };
  for (std::size_t i = 1; i < columns.size(); ++i)
    if (length(columns[i]) != length(columns[0]))
      throw Error("the column '" + columns[i].name + "' holds " +
                  std::to_string(length(columns[i])) + " entries, against " +
                  std::to_string(length(columns[0])) + " in '" + columns[0].name + "'");

  write_table(file, columns, length(columns[0]),
              [&columns] (const std::size_t row, std::vector<double>& values,
                          std::vector<std::int64_t>& integers) {
                for (std::size_t c = 0; c < columns.size(); ++c) {
                  if (columns[c].type == 'K') integers[c] = columns[c].integers[row];
                  else                        values[c] = columns[c].data[row];
                }
              }, options);
}


// ============================================================================


void otswap::io::write (const std::string& file, const std::vector<Column>& columns,
                        const std::size_t nRows,
                        const std::function<void(std::size_t, std::vector<double>&)>& fillRow,
                        const WriteOptions& options)
{
  for (const auto& c : columns)
    if (c.type == 'K')
      throw Error("the column '" + c.name + "' is of type 'K'; write it with the overload whose "
                  "fillRow also takes a buffer of integers");

  write_table(file, columns, nRows,
              [&fillRow] (const std::size_t row, std::vector<double>& values,
                          std::vector<std::int64_t>&) { fillRow(row, values); }, options);
}


// ============================================================================


void otswap::io::write (const std::string& file, const std::vector<Column>& columns,
                        const std::size_t nRows,
                        const std::function<void(std::size_t, std::vector<double>&,
                                                 std::vector<std::int64_t>&)>& fillRow,
                        const WriteOptions& options)
{
  write_table(file, columns, nRows, fillRow, options);
}


// ============================================================================


otswap::BiasTable otswap::io::readBiasTable (const std::string& file,
                                                 const char delimiter, const char comment)
{
  const std::vector<std::string> columns = is_fits(file)
    ? std::vector<std::string>{"REDSHIFT", "BIAS"}
    : std::vector<std::string>{"0", "1"};
  const Table table = read(file, columns, delimiter, comment);

  BiasTable bias;
  bias.redshift.resize(table.nRows);
  bias.bias.resize(table.nRows);
  for (std::size_t r = 0; r < table.nRows; ++r) {
    bias.redshift[r] = table.values[2*r];
    bias.bias[r] = table.values[2*r+1];
  }

  // Rows are counted from 1 over the data, comments and blank lines left out.
  internal::check_bias_table(bias.redshift, bias.bias, "the bias table " + file, "data row", 1);

  return bias;
}
