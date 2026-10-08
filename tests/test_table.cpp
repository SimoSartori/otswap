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
 *  @file tests/test_table.cpp
 *
 *  @brief The table layer: ASCII and FITS round trips, and the failures
 *  the header says are loud.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include <fitsio.h>

#include "otswap/OT.h"
#include "check.h"

using namespace otswap;

namespace {

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "otswap_test_" + name;
  }

  std::string text_of (const std::string& file)
  {
    std::ifstream in(file);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  }

  // The message of the Error f throws, or "" when it throws none.
  template <typename F>
  std::string message_of (F&& f)
  {
    try {
      f();
    }
    catch (const otswap::Error& e) {
      return e.what();
    }
    return "";
  }

  bool same_bits (const double a, const double b)
  {
    return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof(double)) == 0;
  }

  // A header keyword of HDU 2 of a FITS file, with its type code ('C'
  // string, given unquoted, 'I' integer, 'F' floating).
  std::string fits_keyword (const std::string& file, const std::string& name, char& type)
  {
    fitsfile* f = nullptr;
    int status = 0, hdutype = 0;
    char value[FLEN_VALUE] = "", comment[FLEN_COMMENT] = "";
    fits_open_file(&f, file.c_str(), READONLY, &status);
    fits_movabs_hdu(f, 2, &hdutype, &status);
    fits_read_keyword(f, name.c_str(), value, comment, &status);
    char dtype = ' ';
    if (status == 0) fits_get_keytype(value, &dtype, &status);
    type = dtype;
    if (status == 0 && dtype == 'C') fits_read_key(f, TSTRING, name.c_str(), value, nullptr, &status);
    int closing = 0;
    fits_close_file(f, &closing);
    return status == 0 ? std::string(value) : std::string("missing");
  }

  std::vector<io::Column> sample ()
  {
    return {
      {"alpha", 'D', "first column", "", {1.5, -2.25, 3.125, 0.}, {}},
      {"beta",  'D', "",             "", {10., 20., 30., 40.}, {}},
      {"gamma", 'J', "an integer",   "", {1., 2., 3., 4.}, {}}
    };
  }

}

int main ()
{
  const std::vector<io::Column> columns = sample();

  group("an ASCII table round-trips");
  {
    const std::string file = temporary("table.dat");
    io::write(file, columns);

    const io::Table back = io::read(file, {"0", "1", "2"});

    check(back.nRows == 4, "every row is read back");
    check(back.nColumns == 3, "every column is read back");
    check(back.values.size() == 12, "the flat array has one entry per cell");

    for (std::size_t r = 0; r < 4; ++r)
      for (std::size_t c = 0; c < 3; ++c)
        check_close(back.values[r*3 + c], columns[c].data[r], 1.e-8,
                    "the value survives the round trip");

    const io::Table subset = io::read(file, {"2", "0"});
    check(subset.nColumns == 2, "a subset of columns is read");
    check_close(subset.values[0], columns[2].data[0], 1.e-8,
                "the requested order is the returned order");
    check_close(subset.values[1], columns[0].data[0], 1.e-8,
                "and the second requested column follows");

    std::remove(file.c_str());
  }

  group("a FITS table round-trips, by column name");
  {
    const std::string file = temporary("table.fits");
    io::write(file, columns);

    const io::Table back = io::read(file, {"alpha", "beta", "gamma"});

    check(back.nRows == 4, "every row is read back");
    check(back.nColumns == 3, "every column is read back");

    for (std::size_t r = 0; r < 4; ++r)
      for (std::size_t c = 0; c < 3; ++c)
        check_close(back.values[r*3 + c], columns[c].data[r], 1.e-12,
                    "FITS keeps full precision");

    const io::Table insensitive = io::read(file, {"ALPHA"});
    check_close(insensitive.values[0], columns[0].data[0], 1.e-12,
                "column names are matched without regard to case");

    check_throws([&] { io::read(file, {"delta"}); }, "an absent column raises");

    std::remove(file.c_str());
  }

  group("the row-generating writer produces the same file");
  {
    const std::string direct = temporary("direct.dat");
    const std::string generated = temporary("generated.dat");

    io::write(direct, columns);
    io::write(generated, columns, 4,
              [&columns] (const std::size_t row, std::vector<double>& values) {
                for (std::size_t c = 0; c < columns.size(); ++c)
                  values[c] = columns[c].data[row];
              });

    std::ifstream a(direct), b(generated);
    const std::string textA((std::istreambuf_iterator<char>(a)), std::istreambuf_iterator<char>());
    const std::string textB((std::istreambuf_iterator<char>(b)), std::istreambuf_iterator<char>());
    check(textA == textB, "both writers agree byte for byte");

    std::remove(direct.c_str());
    std::remove(generated.c_str());
  }

  group("malformed requests and files are refused");
  {
    const std::string file = temporary("short.dat");
    {
      std::ofstream out(file);
      out << "1 2 3\n4 5\n7 8 9\n";
    }

    // A row with numbers but too few of them is an error, not a silently
    // shortened column.
    check_throws([&] { io::read(file, {"0", "1", "2"}); }, "a short row raises");
    std::remove(file.c_str());

    check_throws([&] { io::read("/nonexistent/otswap/missing.dat", {"0"}); },
                 "a missing file raises");
    check_throws([&] { io::read(temporary("x.dat"), {}); }, "an empty column list raises");
    check_throws([&] { io::read(temporary("x.dat"), {"notaninteger"}); },
                 "a name in a missing ASCII file raises");

    std::vector<io::Column> ragged = sample();
    ragged[1].data.pop_back();
    check_throws([&] { io::write(temporary("ragged.dat"), ragged); },
                 "columns of different length raise");

    check_throws([&] { io::write(temporary("empty.dat"), {}); },
                 "an empty column list raises on write");
  }

  group("in ASCII a NaN is written as nan whatever its sign");
  {
    const std::string file = temporary("nan.dat");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    io::write(file, {{"a", 'D', "", "", {nan, -nan, 1.5}, {}}});
    std::ifstream in(file);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    check(std::signbit(-nan) && text == "###   a\nnan\nnan\n1.5\n", "nan twice, then 1.5: " + text);
    std::remove(file.c_str());
  }

  group("comment and delimiter characters are honoured");
  {
    const std::string file = temporary("commented.dat");
    {
      std::ofstream out(file);
      out << "% a comment\n1,2,3\n% another\n4,5,6\n";
    }

    const io::Table back = io::read(file, {"0", "2"}, ',', '%');
    check(back.nRows == 2, "commented lines are skipped");
    check_close(back.values[0], 1., 1.e-12, "the first field is read");
    check_close(back.values[1], 3., 1.e-12, "the third field is read");

    std::remove(file.c_str());
  }

  group("integer ('K') columns write and read back exactly, near both ends of the range");
  {
    const std::int64_t lo = std::numeric_limits<std::int64_t>::min();
    const std::int64_t hi = std::numeric_limits<std::int64_t>::max();
    const std::vector<std::int64_t> ids {lo, lo + 1, -1, 0, 9007199254740993LL, hi - 1, hi};
    const std::vector<double> x {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
    std::vector<io::Column> table {
      {"x", 'D', "", "", x, {}},
      {"id", 'K', "an identifier", "", {}, ids}
    };

    for (const std::string name : {"ints.dat", "ints.fits"}) {
      const std::string file = temporary(name);
      io::write(file, table);
      const io::Table back = io::read(file, {"x"}, {"id"});
      check(back.nRows == ids.size() && back.nColumns == 1 && back.nIntegerColumns == 1,
            name + ": the sizes of both parts");
      bool exact = back.integers.size() == ids.size();
      for (std::size_t r = 0; exact && r < ids.size(); ++r)
        exact = back.integers[r] == ids[r] && back.values[r] == x[r];
      check(exact, name + ": every integer and value read back exactly");

      const io::Table onlyIntegers = io::read(file, {}, {"id"});
      check(onlyIntegers.nColumns == 0 && onlyIntegers.values.empty() &&
            onlyIntegers.integers == ids, name + ": an integer column alone");
      std::remove(file.c_str());
    }

    const std::string file = temporary("ints.dat");
    io::write(file, table);
    check(text_of(file) ==
          "# id   an identifier\n#\n###   x   id\n"
          "0.5 -9223372036854775808\n1.5 -9223372036854775807\n2.5 -1\n3.5 0\n"
          "4.5 9007199254740993\n5.5 9223372036854775806\n6.5 9223372036854775807\n",
          "the ASCII text of a 'K' column: " + text_of(file));

    const io::Table asDouble = io::read(file, {"1"});
    check(asDouble.values[4] == 9007199254740992., "read as double, a large integer rounds");
    std::remove(file.c_str());

    table[1].integers.pop_back();
    check(message_of([&] { io::write(file, table); }).find("'id' holds 6") != std::string::npos,
          "a 'K' column of another length is refused, by its integers");

    std::vector<io::Column> withData {{"id", 'K', "", "", {1., 2.}, {7, 8, 9}}, {"x", 'D', "", "", {1., 2., 3.}, {}}};
    io::write(file, withData);
    check(io::read(file, {}, {"id"}).integers == std::vector<std::int64_t>{7, 8, 9},
          "the data of a 'K' column is ignored");
    std::remove(file.c_str());
  }

  group("the row writer with an integer buffer gives the same file as the direct one");
  {
    const std::vector<io::Column> table {
      {"x", 'D', "a value", "Mpc/h", {1.25, -2.5, 3.}, {}},
      {"n", 'J', "", "", {3., 2., 1.}, {}},
      {"id", 'K', "", "", {}, {10, 20, 30}}
    };
    for (const std::string name : {"rows.dat", "rows.fits"}) {
      const std::string direct = temporary("direct_" + name), generated = temporary("generated_" + name);
      io::write(direct, table);
      io::write(generated, table, 3,
                [&table] (const std::size_t row, std::vector<double>& values,
                          std::vector<std::int64_t>& integers) {
                  values[0] = table[0].data[row];
                  values[1] = table[1].data[row];
                  values[2] = -1.;
                  integers[0] = -1;
                  integers[1] = -1;
                  integers[2] = table[2].integers[row];
                });
      if (name == "rows.dat")
        check(text_of(direct) == text_of(generated), "both ASCII writers agree byte for byte");
      const io::Table a = io::read(direct, {"x", "n"}, {"id"});
      const io::Table b = io::read(generated, {"x", "n"}, {"id"});
      check(a.values == b.values && a.integers == b.integers, name + ": the same table back");
      std::remove(direct.c_str());
      std::remove(generated.c_str());
    }

    check(message_of([&] {
            io::write(temporary("k.dat"), table, 3,
                      [] (const std::size_t, std::vector<double>&) {});
          }).find("'id' is of type 'K'") != std::string::npos,
          "a 'K' column is refused by the row writer without an integer buffer");
  }

  group("a 'J' value that is not a 32-bit integer is refused, naming the column and the row");
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const std::string name : {"j.dat", "j.fits"}) {
      const std::string file = temporary(name);
      for (const double bad : {nan, 1.5, 2147483648., -2147483649., std::numeric_limits<double>::infinity()}) {
        const std::string m = message_of([&] {
          io::write(file, {{"count", 'J', "", "", {1., bad, 3.}, {}}});
        });
        check(m.find("'J' column 'count'") != std::string::npos && m.find("at row 1") != std::string::npos,
              name + ": refused: " + m);
      }
      io::write(file, {{"count", 'J', "", "", {-2147483648., 2147483647., 0.}, {}}});
      const io::Table back = io::read(file, {}, {name == "j.dat" ? "0" : "count"});
      check(back.integers == std::vector<std::int64_t>{-2147483648LL, 2147483647LL, 0},
            name + ": both ends of the 32-bit range are written");
      std::remove(file.c_str());
    }
  }

  group("ASCII integer fields that are not integers are refused, naming the line");
  {
    const std::string file = temporary("badints.dat");
    for (const std::string bad : {"1.5", "+3", "nan", "1e3", "9223372036854775808"}) {
      {
        std::ofstream out(file);
        out << "# a comment\n1 2\n3 " << bad << "\n";
      }
      const std::string m = message_of([&] { io::read(file, {"0"}, {"1"}); });
      check(m.find("holds " + bad + " at line 3") != std::string::npos, "refused: " + m);
    }
    {
      std::ofstream out(file);
      out << "1 -9223372036854775808\n  2   -0\n";
    }
    check(io::read(file, {"0"}, {"1"}).integers ==
          std::vector<std::int64_t>{std::numeric_limits<std::int64_t>::min(), 0},
          "the most negative integer and -0 are read");
    std::remove(file.c_str());
  }

  group("FITS: integer requests need an integer column, and a null integer is an error");
  {
    const std::string file = temporary("nulls.fits");
    io::write(file, {{"x", 'D', "", "", {1., 2.}, {}}});
    check(message_of([&] { io::read(file, {}, {"x"}); }).find("not integers") != std::string::npos,
          "a floating column is refused as an integer column");

    fitsfile* f = nullptr;
    int status = 0;
    char ttype[] = "n", tform[] = "1J";
    char* types[] = {ttype};
    char* forms[] = {tform};
    const std::string overwrite = "!" + file;
    fits_create_file(&f, overwrite.c_str(), &status);
    fits_create_tbl(f, BINARY_TBL, 3, 1, types, forms, nullptr, nullptr, &status);
    long null = -99;
    fits_write_key(f, TLONG, "TNULL1", &null, nullptr, &status);
    long values[] = {1, -99, 3};
    fits_write_col(f, TLONG, 1, 1, 1, 3, values, &status);
    fits_close_file(f, &status);
    check(status == 0, "the file with a null integer is written");

    const std::string asInteger = message_of([&] { io::read(file, {}, {"n"}); });
    const std::string asDouble = message_of([&] { io::read(file, {"n"}); });
    check(asInteger.find("undefined value in the column 'n', row 2") != std::string::npos,
          "read as an integer, the null is an error: " + asInteger);
    check(asDouble.find("undefined value in the column 'n', row 2") != std::string::npos,
          "read as double, too: " + asDouble);
    std::remove(file.c_str());
  }

  group("NaN and infinities read back as stored, in both formats");
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const std::vector<double> x {nan, -nan, inf, -inf, 0.1};
    for (const std::string name : {"special.dat", "special.fits"}) {
      const std::string file = temporary(name);
      io::write(file, {{"x", 'D', "", "", x, {}}}, {17, {}});
      const io::Table back = io::read(file, {name == "special.dat" ? "0" : "x"});
      bool same = back.nRows == x.size();
      for (std::size_t r = 0; same && r < x.size(); ++r)
        same = r < 2 ? std::isnan(back.values[r]) : back.values[r] == x[r];
      check(same, name + ": NaN, +inf, -inf and 0.1 read back");
      std::remove(file.c_str());
    }
  }

  group("17 significant digits read back as the same double in ASCII");
  {
    std::vector<double> x;
    double v = 0.1;
    for (int i = 0; i < 2000; ++i) {
      v = v * 3.7 + 0.123456789;
      if (v > 1.e6) v = v / 9.87654321e5 - 0.5;
      x.push_back(i % 3 == 0 ? -v : v);
    }
    x.push_back(std::numeric_limits<double>::denorm_min());
    x.push_back(std::numeric_limits<double>::max());
    x.push_back(-0.);
    const std::string file = temporary("exact.dat");
    io::write(file, {{"x", 'D', "", "", x, {}}}, {17, {}});
    const io::Table back = io::read(file, {"x"});
    bool exact = back.nRows == x.size();
    for (std::size_t r = 0; exact && r < x.size(); ++r) exact = same_bits(back.values[r], x[r]);
    check(exact, "every double, -0 included, is the same bits");

    io::write(file, {{"x", 'D', "", "", {1. / 3.}, {}}});
    check(text_of(file) == "###   x\n0.333333333\n", "the default is 9 digits");
    io::write(file, {{"x", 'D', "", "", {1. / 3.}, {}}}, {4, {}});
    check(text_of(file) == "###   x\n0.3333\n", "the precision is honoured");

    for (const int p : {0, 18, -1})
      check_throws([&] { io::write(file, {{"x", 'D', "", "", {1.}, {}}}, {p, {}}); },
                   "a precision outside [1, 17] is refused");
    std::remove(file.c_str());
  }

  group("units and keywords in the ASCII header");
  {
    const std::string file = temporary("header.dat");
    io::WriteOptions options;
    options.keywords = {{"PRODUCT", "displacements", "what the file holds"},
                        {"NOBJECTS", "3", ""},
                        {"MPS", "10.817912345678901", "Mpc/h"}};
    io::write(file, {{"tracX", 'D', "position", "Mpc/h", {1., 2., 3.}, {}},
                     {"n", 'J', "", "", {1., 2., 3.}, {}},
                     {"z", 'D', "", "", {0.5, 0.6, 0.7}, {}},
                     {"lagrRA", 'D', "", "deg", {10., 20., 30.}, {}}}, options);
    check(text_of(file) ==
          "## PRODUCT = displacements / what the file holds\n"
          "## NOBJECTS = 3\n"
          "## MPS = 10.817912345678901 / Mpc/h\n"
          "# tracX    position [Mpc/h]\n"
          "# lagrRA   [deg]\n"
          "#\n"
          "###   tracX   n   z   lagrRA\n"
          "1 1 0.5 10\n2 2 0.6 20\n3 3 0.7 30\n",
          "the header: " + text_of(file));
    const io::Table back = io::read(file, {"lagrRA", "z"});
    check(back.nRows == 3 && back.values[0] == 10. && back.values[1] == 0.5,
          "a file with keywords and units reads back by name");

    io::write(file, {{"a", 'D', "", "", {1.}, {}}}, {9, {{"SEED", "12345", ""}}});
    check(text_of(file) == "## SEED = 12345\n#\n###   a\n1\n",
          "keywords alone are followed by the # line: " + text_of(file));

    for (const std::string bad : {"", "TOOLONGNAME", "lower", "A B", "A.B"})
      check_throws([&] { io::write(file, {{"a", 'D', "", "", {1.}, {}}}, {9, {{bad, "1", ""}}}); },
                   "the keyword name '" + bad + "' is refused");
    check_throws([&] { io::write(file, {{"a", 'D', "", "", {1.}, {}}}, {9, {{"A", "x\ny", ""}}}); },
                 "a line break in a keyword is refused");
    check_throws([&] { io::write(file, {{"a", 'D', "", "", {1.}, {}}},
                                 {9, {{"A", std::string(69, 'v'), ""}}}); },
                 "a keyword value too long for a FITS card is refused");
    check_throws([&] { io::write(file, {{"a", 'X', "", "", {1.}, {}}}); },
                 "an unknown column type is refused");
    std::remove(file.c_str());
  }

  group("units and keywords in FITS: TUNITn, and numbers written as numbers");
  {
    const std::string file = temporary("header.fits");
    io::WriteOptions options;
    options.keywords = {{"PRODUCT", "displacements", "what the file holds"},
                        {"NOBJECTS", "29579", ""},
                        {"MPS", "1.0817912345678901e+01", "Mpc/h"},
                        {"ZCUTMIN", "-0.5", ""},
                        {"OTSWAPV", "0.1.0", ""}};
    io::write(file, {{"tracX", 'D', "position", "Mpc/h", {1., 2.}, {}},
                     {"id", 'K', "", "", {}, {5, 6}}}, options);
    char type = ' ';
    check(fits_keyword(file, "TUNIT1", type) == "Mpc/h" && type == 'C', "the unit is TUNIT1");
    check(fits_keyword(file, "TCOMM1", type) == "position", "the description is TCOMM1");
    check(fits_keyword(file, "TUNIT2", type) == "missing", "no unit, no TUNIT2");
    check(fits_keyword(file, "PRODUCT", type) == "displacements" && type == 'C', "a string");
    check(fits_keyword(file, "NOBJECTS", type) == "29579" && type == 'I', "an integer");
    const std::string mps = fits_keyword(file, "MPS", type);
    check(type == 'F' && std::strtod(mps.c_str(), nullptr) == 10.817912345678901,
          "a number, which reads back as the same double: " + mps);
    check(fits_keyword(file, "ZCUTMIN", type) == "-0.5" && type == 'F', "a negative number");
    check(fits_keyword(file, "OTSWAPV", type) == "0.1.0" && type == 'C',
          "a version is not a number");
    const io::Table back = io::read(file, {"tracX"}, {"id"});
    check(back.values == std::vector<double>{1., 2.} && back.integers == std::vector<std::int64_t>{5, 6},
          "the table reads back");
    std::remove(file.c_str());
  }

  group("ASCII columns by name, from the ### line");
  {
    const std::string file = temporary("named.dat");
    {
      std::ofstream out(file);
      out << "### old names\n# a comment\n###   a   b   c\n1 2 3\n4 5 6\n### later\n7 8 9\n";
    }
    const io::Table back = io::read(file, {"c", "0"});
    check(back.nRows == 3 && back.values[0] == 3. && back.values[1] == 1. && back.values[4] == 9.,
          "names and indices mix, from the last ### line before the data");
    const std::string m = message_of([&] { io::read(file, {"d"}); });
    check(m.find("'d' was not found") != std::string::npos && m.find("names: a, b, c") != std::string::npos,
          "an absent name lists the names found: " + m);
    check(message_of([&] { io::read(file, {"a", "0"}); }).find("more than once") != std::string::npos,
          "a column requested twice, by name and index, is refused");
    check(message_of([&] { io::read(file, {"a"}, {"a"}); }).find("more than once") != std::string::npos,
          "as double and as integer, too");

    {
      std::ofstream out(file);
      out << "# x y\n1 2\n";
    }
    check(message_of([&] { io::read(file, {"x"}); }).find("no line starting with ###") != std::string::npos,
          "a name without a ### line is refused");
    check_throws([&] { io::read(file, {"4096"}); }, "an index above 4095 is refused");
    std::remove(file.c_str());

    const std::string fits = temporary("named.fits");
    io::write(fits, sample());
    check(message_of([&] { io::read(fits, {"alpha", "ALPHA"}); }).find("more than once") != std::string::npos,
          "FITS: a column requested twice is refused");
    std::remove(fits.c_str());
  }

  return report("test_table");
}
