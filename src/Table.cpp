/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
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

#include <fitsio.h>

#include "internal.h"

namespace {

  // Largest 0-based column index an ASCII request may name. A larger one
  // is far more likely to be a mistyped column list than a real table.
  constexpr std::size_t kMaxAsciiColumnIndex = 4095;

  // Significant digits written to ASCII, the same as CosmoBolognaLib's
  // table writer, so that files from either are formatted alike. Nine
  // digits do not round-trip a double exactly.
  constexpr int kAsciiPrecision = 9;

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

  otswap::io::Table read_ascii (const std::string& file,
                                const std::vector<std::string>& columns,
                                const char delimiter, const char comment)
  {
    std::ifstream input(file);
    if (!input)
      throw otswap::Error("cannot open the input file: " + file);

    std::vector<std::size_t> index(columns.size());
    for (std::size_t i = 0; i < columns.size(); ++i) {
      const std::string& c = columns[i];
      std::size_t value = 0;
      const auto conversion = std::from_chars(c.data(), c.data()+c.size(), value);
      if (c.empty() || conversion.ec != std::errc() || conversion.ptr != c.data()+c.size())
        throw otswap::Error("invalid column index for an ASCII file: '" + c +
                            "'; a non-negative integer was expected");
      if (value > kMaxAsciiColumnIndex)
        throw otswap::Error("the ASCII column index " + std::to_string(value) +
                            " is too large; the maximum supported index is " +
                            std::to_string(kMaxAsciiColumnIndex));
      index[i] = value;
    }

    std::size_t maxIndex = 0;
    for (const std::size_t c : index) maxIndex = std::max(maxIndex, c);

    constexpr std::size_t notRequested = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> map(maxIndex+1, notRequested);
    for (std::size_t i = 0; i < index.size(); ++i) {
      if (map[index[i]] != notRequested)
        throw otswap::Error("the ASCII column index " + std::to_string(index[i]) +
                            " is requested more than once");
      map[index[i]] = i;
    }

    otswap::io::Table table;
    table.nColumns = columns.size();

    std::vector<double> row(columns.size());
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
        if (map[column] != notRequested) row[map[column]] = value;
        ++fields;
        p = end;
        while (*p != '\0' && (*p == delimiter || *p == ' ' || *p == '\t')) ++p;
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
      ++table.nRows;
    }

    return table;
  }

  otswap::io::Table read_fits (const std::string& file,
                               const std::vector<std::string>& columns)
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

    const std::size_t ncolumns = columns.size();
    std::vector<int> colnum(ncolumns, 0);

    for (std::size_t i = 0; i < ncolumns; ++i) {
      std::vector<char> name(columns[i].begin(), columns[i].end());
      name.push_back('\0');
      if (fits_get_colnum(fptr, CASEINSEN, name.data(), &colnum[i], &status))
        fail("the column '" + columns[i] + "' was not found in: " + file);

      int typecode = 0;
      long repeat = 0, width = 0;
      if (fits_get_coltype(fptr, colnum[i], &typecode, &repeat, &width, &status))
        fail("cannot read the type of the column '" + columns[i] + "' in: " + file);

      if (repeat != 1)
        fail("the column '" + columns[i] + "' of " + file + " has repeat count " +
             std::to_string(repeat) + ": vector columns are not supported");
    }

    std::vector<std::vector<double>> byColumn(ncolumns, std::vector<double>(nrows));

    long chunk = 0;
    if (fits_get_rowsize(fptr, &chunk, &status)) { status = 0; chunk = 0; }
    if (chunk < 1) chunk = 1;

    double nulval = std::numeric_limits<double>::quiet_NaN();
    std::vector<char> nullFlags(chunk);

    for (long first = 1; first <= nrows; first += chunk) {
      const long nread = std::min(chunk, nrows-first+1);
      for (std::size_t i = 0; i < ncolumns; ++i) {
        int anynul = 0;
        if (fits_read_col(fptr, TDOUBLE, colnum[i], first, 1, nread, &nulval,
                          byColumn[i].data()+first-1, &anynul, &status))
          fail("cannot read the column '" + columns[i] + "' of: " + file);

        if (anynul) {
          std::fill(nullFlags.begin(), nullFlags.begin()+nread, 0);
          int detailed = 0;
          if (fits_read_colnull(fptr, TDOUBLE, colnum[i], first, 1, nread,
                                byColumn[i].data()+first-1, nullFlags.data(), &detailed, &status))
            fail("cannot locate the undefined value in the column '" + columns[i] + "' of: " + file);

          const auto at = std::find_if(nullFlags.begin(), nullFlags.begin()+nread,
                                       [] (const char flag) { return flag != 0; });
          const long offset = at == nullFlags.begin()+nread ? 0 : std::distance(nullFlags.begin(), at);
          fail("undefined or non-finite value in the column '" + columns[i] + "', row " +
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
    table.nColumns = ncolumns;
    table.values.resize(table.nRows * ncolumns);
    for (std::size_t r = 0; r < table.nRows; ++r)
      for (std::size_t c = 0; c < ncolumns; ++c)
        table.values[r*ncolumns + c] = byColumn[c][r];

    return table;
  }

  void write_ascii (const std::string& file,
                    const std::vector<otswap::io::Column>& columns,
                    const std::size_t nrows,
                    const std::function<void(std::size_t, std::vector<double>&)>& fillRow)
  {
    std::ofstream out(file.c_str());
    if (!out)
      throw otswap::Error("cannot open the output file: " + file);

    std::size_t width = 0;
    bool anyDescription = false;
    for (const auto& c : columns) {
      if (c.description.empty()) continue;
      anyDescription = true;
      width = std::max(width, c.name.size());
    }

    if (anyDescription) {
      for (const auto& c : columns) {
        if (c.description.empty()) continue;
        out << "# " << c.name << std::string(width - c.name.size() + 3, ' ')
            << c.description << std::endl;
      }
      out << "#" << std::endl;
    }

    out << "###";
    for (const auto& c : columns) out << "   " << c.name;
    out << std::endl;

    out << std::defaultfloat << std::setprecision(kAsciiPrecision);
    std::vector<double> values(columns.size());
    for (std::size_t row = 0; row < nrows; ++row) {
      fillRow(row, values);
      for (std::size_t c = 0; c < columns.size(); ++c) {
        if (c > 0) out << " ";
        if (columns[c].type == 'J') out << (long long)values[c];
        else                        out << values[c];
      }
      out << std::endl;
    }

    out.close();
    if (out.fail())
      throw otswap::Error("cannot complete writing the output file: " + file);
  }

  void write_fits (const std::string& file,
                   const std::vector<otswap::io::Column>& columns,
                   const std::size_t nrows,
                   const std::function<void(std::size_t, std::vector<double>&)>& fillRow)
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
      formStore[i] = (columns[i].type == 'J') ? "1J" : "1D";
      ttype[i] = &nameStore[i][0];
      tform[i] = &formStore[i][0];
    }

    if (nrows > (std::size_t)std::numeric_limits<long>::max())
      fail("too many rows for a FITS table in: " + file);

    const long fitsRows = (long)nrows;
    if (fits_create_tbl(fptr, BINARY_TBL, fitsRows, ncolumns, ttype.data(), tform.data(),
                        nullptr, nullptr, &status))
      fail("cannot create the binary table in: " + file);

    for (int i = 0; i < ncolumns; ++i) {
      if (columns[i].description.empty()) continue;
      const std::string key = "TCOMM" + std::to_string(i+1);
      std::vector<char> value(columns[i].description.begin(), columns[i].description.end());
      value.push_back('\0');
      if (fits_write_key(fptr, TSTRING, key.c_str(), value.data(), nullptr, &status))
        fail("cannot write the description of the column '" + columns[i].name + "' in: " + file);
    }

    long chunk = 0;
    if (fits_get_rowsize(fptr, &chunk, &status)) { status = 0; chunk = 0; }
    if (chunk < 1) chunk = 1;

    std::vector<std::vector<double>> buffer(columns.size(), std::vector<double>(chunk));
    std::vector<long> integers(chunk);
    std::vector<double> values(columns.size());

    try {
      for (long first = 1; first <= fitsRows; first += chunk) {
        const long nwrite = std::min(chunk, fitsRows-first+1);
        for (long row = 0; row < nwrite; ++row) {
          fillRow((std::size_t)(first-1+row), values);
          for (std::size_t c = 0; c < columns.size(); ++c)
            buffer[c][row] = values[c];
        }
        for (int i = 0; i < ncolumns; ++i) {
          if (columns[i].type == 'J') {
            for (long row = 0; row < nwrite; ++row)
              integers[row] = (long)buffer[i][row];
            if (fits_write_col(fptr, TLONG, i+1, first, 1, nwrite, integers.data(), &status))
              fail("cannot write the column '" + columns[i].name + "' in: " + file);
          }
          else {
            if (fits_write_col(fptr, TDOUBLE, i+1, first, 1, nwrite, buffer[i].data(), &status))
              fail("cannot write the column '" + columns[i].name + "' in: " + file);
          }
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

}


// ============================================================================


otswap::io::Table otswap::io::read (const std::string& file,
                                    const std::vector<std::string>& columns,
                                    const char delimiter, const char comment)
{
  if (columns.empty())
    throw Error("no column was requested from: " + file);

  if (is_fits(file)) return read_fits(file, columns);
  return read_ascii(file, columns, delimiter, comment);
}


// ============================================================================


void otswap::io::write (const std::string& file, const std::vector<Column>& columns)
{
  if (columns.empty())
    throw Error("no column was given for: " + file);

  for (std::size_t i = 1; i < columns.size(); ++i)
    if (columns[i].data.size() != columns[0].data.size())
      throw Error("the column '" + columns[i].name + "' holds " +
                  std::to_string(columns[i].data.size()) + " entries, against " +
                  std::to_string(columns[0].data.size()) + " in '" + columns[0].name + "'");

  write(file, columns, columns[0].data.size(),
        [&columns] (const std::size_t row, std::vector<double>& values) {
          for (std::size_t c = 0; c < columns.size(); ++c)
            values[c] = columns[c].data[row];
        });
}


// ============================================================================


void otswap::io::write (const std::string& file, const std::vector<Column>& columns,
                        const std::size_t nRows,
                        const std::function<void(std::size_t, std::vector<double>&)>& fillRow)
{
  if (columns.empty())
    throw Error("no column was given for: " + file);

  if (is_fits(file)) write_fits(file, columns, nRows, fillRow);
  else               write_ascii(file, columns, nRows, fillRow);
}
