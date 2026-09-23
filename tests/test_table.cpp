/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
 ********************************************************************/

/**
 *  @file tests/test_table.cpp
 *
 *  @brief The table layer: ASCII and FITS round trips, and the failures
 *  the header says are loud.
 */

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

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

  std::vector<io::Column> sample ()
  {
    return {
      {"alpha", 'D', "first column", {1.5, -2.25, 3.125, 0.}},
      {"beta",  'D', "",             {10., 20., 30., 40.}},
      {"gamma", 'J', "an integer",   {1., 2., 3., 4.}}
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
                 "a non-integer ASCII column raises");

    std::vector<io::Column> ragged = sample();
    ragged[1].data.pop_back();
    check_throws([&] { io::write(temporary("ragged.dat"), ragged); },
                 "columns of different length raise");

    check_throws([&] { io::write(temporary("empty.dat"), {}); },
                 "an empty column list raises on write");
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

  return report("test_table");
}
