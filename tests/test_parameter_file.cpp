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
 *  @file tests/test_parameter_file.cpp
 *
 *  @brief The parameter-file reader: the syntax it accepts and refuses,
 *  the typed reads and their messages, and the marks behind unread().
 */

#include <clocale>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "otswap/ParameterFile.h"
#include "check.h"

using otswap::io::ParameterFile;

namespace {

  using Strings = std::vector<std::string>;

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "otswap_test_parameter_file_" + name;
  }

  /// Writes the bytes of @p text, unchanged, to a temporary file of that
  /// name, and returns its path.
  std::string file_with (const std::string& name, const std::string& text)
  {
    const std::string path = temporary(name);
    std::ofstream out(path, std::ios::binary);
    out << text;
    return path;
  }

  /// @p f must raise an otswap::Error whose message is exactly @p text.
  template <typename F>
  void check_message (F&& f, const std::string& text, const std::string& what)
  {
    try {
      f();
    }
    catch (const otswap::Error& e) {
      if (std::string(e.what()) != text) {
        std::cerr << "FAILED: " << what << " (the message is \"" << e.what() << "\", not \""
                  << text << "\")" << std::endl;
        ++failures;
      }
      return;
    }
    catch (const std::exception& e) {
      std::cerr << "FAILED: " << what << " (threw the wrong type: " << e.what() << ")"
                << std::endl;
      ++failures;
      return;
    }
    std::cerr << "FAILED: " << what << " (did not throw)" << std::endl;
    ++failures;
  }

  /// The constructor on a file of @p text must raise an Error whose message
  /// is the file's path followed by @p rest.
  void check_refused (const std::string& name, const std::string& text, const std::string& rest,
                      const std::string& what)
  {
    const std::string path = file_with(name, text);
    check_message([&] { ParameterFile p(path); }, path + rest, what);
  }

  /// Reading @p key of a file holding the one line "k = <value>" as a T
  /// must raise an Error whose message is "<path>:1: key 'k': " + @p rest.
  template <typename T>
  void check_read_refused (const std::string& value, const std::string& rest, const std::string& what)
  {
    const std::string path = file_with("read_refused", "k = " + value + "\n");
    const ParameterFile p(path);
    check_message([&] { p.get<T>("k"); }, path + ":1: key 'k': " + rest, what);
  }

  bool same_bits (const double a, const double b)
  {
    return std::memcmp(&a, &b, sizeof a) == 0;
  }

  bool same_entries (const std::vector<ParameterFile::Entry>& entries,
                     const std::vector<std::pair<std::string, std::size_t>>& expected)
  {
    if (entries.size() != expected.size()) return false;
    for (std::size_t i = 0; i < entries.size(); ++i)
      if (entries[i].key != expected[i].first || entries[i].line != expected[i].second) return false;
    return true;
  }

  /// The elements of @p key, read as a std::vector<std::string>.
  Strings elements (const ParameterFile& p, const std::string& key)
  {
    return p.get<Strings>(key);
  }

}

int main ()
{
  group("the accepted syntax: blanks, comments, tabs, CRLF, a BOM, no final newline");
  {
    const std::string text =
      "\xEF\xBB\xBF# a comment on line 1, after the BOM\r\n"
      "\r\n"
      "   # an indented comment\r\n"
      "\t\r\n"
      "a = 1 # a comment after the value\r\n"
      "b=2\r\n"
      "\tc\t=\t3\t4\t\r\n"
      "d = x";
    const ParameterFile p(file_with("accepted", text));
    check(elements(p, "a") == Strings{"1"}, "a comment after a value is not part of it");
    check(elements(p, "b") == Strings{"2"}, "k=v without spaces is read");
    check(elements(p, "c") == Strings{"3", "4"}, "tabs separate as spaces do, also around '='");
    check(elements(p, "d") == Strings{"x"}, "a last line without a newline is read");
    check(p.get<int>("a") == 1, "the trailing \\r is not part of the value");
    const ParameterFile q(file_with("accepted_lines", "\xEF\xBB\xBF" "first = 1\n\n# c\nsecond = 2\n"));
    check(same_entries(q.unread(), {{"first", 1}, {"second", 4}}),
          "a BOM before the first key is removed, and lines count from 1");
  }

  group("an empty file and a file of comments only are valid");
  {
    for (const std::string& text : {std::string(""), std::string("# one\n\n  # two\n")}) {
      const ParameterFile p(file_with("empty", text));
      check(!p.has("k"), "no key is set");
      check(p.unread().empty(), "nothing is left unread");
    }
  }

  group("the elements of a value, the table of the signatures note");
  {
    const std::string text =
      "Nrec = 8\n"
      "column_tracer = 0 1 2\n"
      "braces = {0, 1, 2}\n"
      "name = \"a b\" c\n"
      "hash = \"#\"\n"
      "comment = #\n"
      "space = \" \"\n"
      "single = ' '\n"
      "tab = \\t\n"
      "none =\n"
      "empty = \"\"\n"
      "empties = \"\" \"\"\n"
      "output_dir = ../output/ \n"
      "glued_comment = \"a\"#c\n"
      "unquoted_hash = a#b\n";
    const ParameterFile p(file_with("elements", text));
    check(elements(p, "Nrec") == Strings{"8"}, "one element");
    check(elements(p, "column_tracer") == Strings{"0", "1", "2"}, "three elements");
    check(elements(p, "braces") == Strings{"{0,", "1,", "2}"}, "braces and commas are ordinary");
    check(elements(p, "name") == Strings{"a b", "c"}, "quotes group spaces");
    check(elements(p, "hash") == Strings{"#"}, "a quoted # is a character");
    check(elements(p, "comment").empty(), "an unquoted # opens a comment");
    check(elements(p, "space") == Strings{" "}, "a quoted space");
    check(elements(p, "single") == Strings{"'", "'"}, "single quotes are ordinary");
    check(elements(p, "tab") == Strings{"\\t"}, "\\t is two characters as a string");
    check(elements(p, "none").empty(), "k = is an empty list");
    check(elements(p, "empty") == Strings{""}, "k = \"\" is one empty element");
    check(elements(p, "empties") == Strings{"", ""}, "two empty elements");
    check(elements(p, "output_dir") == Strings{"../output/"}, "a trailing space is not kept");
    check(elements(p, "glued_comment") == Strings{"a"}, "a comment may follow a quote");
    check(elements(p, "unquoted_hash") == Strings{"a"}, "# ends an unquoted element");
    check(p.get<char>("tab") == '\t', "\\t is a tab as a char");
    check(p.get<std::string>("tab") == "\\t", "and two characters as a std::string");
  }

  group("'=' after the first belongs to the value");
  {
    const ParameterFile p(file_with("equals", "k = a=b\nj = = x\n"));
    check(elements(p, "k") == Strings{"a=b"}, "k = a=b");
    check(elements(p, "j") == Strings{"=", "x"}, "k = = x");
  }

  group("quotes group and never change the meaning");
  {
    const ParameterFile p(file_with("quotes",
      "a = \"a b\"\nb = \"'\"\nc = \"\t\"\nd = \"5\"\ne = \"true\"\nf = \"1.5\"\n"));
    check(p.get<std::string>("a") == "a b", "a quoted space");
    check(p.get<char>("b") == '\'', "a quoted single quote");
    check(p.get<char>("c") == '\t', "a quoted tab");
    check(p.get<int>("d") == 5, "\"5\" is 5");
    check(p.get<bool>("e"), "\"true\" is true");
    check(p.get<double>("f") == 1.5, "\"1.5\" is 1.5");
  }

  group("keys: letters, digits and _, case-sensitive");
  {
    const ParameterFile p(file_with("keys", "_x = 1\nx_1 = 2\nX = 3\nx = 4\ny = 5\n"));
    check(p.get<int>("_x") == 1 && p.get<int>("x_1") == 2 && p.get<int>("y") == 5,
          "keys of _, digits and one letter");
    check(p.get<int>("X") == 3 && p.get<int>("x") == 4, "X and x are two keys");
    check(!p.has("my key") && !p.has("X "), "a name that is not a key is not found");
    check_message([&] { p.get<int>("my key"); },
                  temporary("keys") + ": key 'my key' is missing", "nor read");
  }

  group("files that cannot be opened or read");
  {
    const std::string missing = temporary("does_not_exist");
    std::remove(missing.c_str());
    check_message([&] { ParameterFile p(missing); }, missing + ": cannot open the parameter file",
                  "a missing file");
    std::string directory = temporary("");
    directory = directory.substr(0, directory.rfind('/'));
    if (std::FILE* f = std::fopen(directory.c_str(), "rb")) {
      std::fclose(f);
      check_message([&] { ParameterFile p(directory); }, directory + ": cannot read the parameter file",
                    "a directory, which this platform opens");
    }
    else
      std::cout << "    (this platform does not open a directory as a file; not checked)" << std::endl;
  }

  group("the refused syntax, with the whole message");
  {
    const std::string valid = ": not a valid key; a key holds letters, digits and _, and does not "
                              "start with a digit";
    const std::string glued = ": text glued to a quote; separate the elements with a space";
    const std::string expected = ":2: expected 'key = value', a comment or a blank line";
    check_refused("e2a", "a = 1\nseed 5\n", expected, "a line without '='");
    check_refused("e2b", "a = 1\nseed # = 5\n", expected, "'=' inside a comment");
    check_refused("e2c", "a = 1\n\"\n", expected, "a line of one quote");
    check_refused("e3a", "= 5\n", ":1: no key before '='", "no key");
    check_refused("e3b", "  = 5\n", ":1: no key before '='", "only blanks before '='");
    check_refused("e4a", "my key = 1\n", ":1: key 'my key'" + valid, "a key with a space");
    check_refused("e4b", "2pt = 1\n", ":1: key '2pt'" + valid, "a key starting with a digit");
    check_refused("e4c", "k\" = 1\n", ":1: key 'k\"'" + valid, "a key with a quote");
    check_refused("e4d", "a = 1\n\xEF\xBB\xBF" "b = 2\n", ":2: key '\xEF\xBB\xBF" "b'" + valid,
                  "a BOM on line 2 is part of the key");
    check_refused("e4e", "k\r = 1\n", ":1: key 'k\r'" + valid, "a \\r inside the key");
    check_refused("e5a", "k = \"abc\n", ":1: key 'k': unterminated quote", "a quote not closed");
    check_refused("e5b", "k = a \"b\n", ":1: key 'k': unterminated quote", "a second element not closed");
    check_refused("e6a", "k = \"a\"b\n", ":1: key 'k'" + glued, "text after a quote");
    check_refused("e6b", "k = \"a\"\"b\"\n", ":1: key 'k'" + glued, "two quoted parts glued");
    check_refused("e6c", "k = a\"b\"\n", ":1: key 'k'" + glued, "a quote inside an element");
    check_refused("e6d", "k = a\"\n", ":1: key 'k'" + glued, "a quote ending an element");
    check_refused("e7", "seed = 1\nx = 2\nseed = 3\nseed = 4\n",
                  ":3: key 'seed': repeated; first set at line 1", "a key set twice");
    check_refused("e7syntax", "seed = 1\nseed = \"x\n", ":2: key 'seed': unterminated quote",
                  "a repeated key with a broken value reports the value");
    check_refused("first", "a = 1\nb = 2\nc = \"x\"y\nd = 4\ne 5\n", ":3: key 'c'" + glued,
                  "the first error in file order wins");
  }

  group("each scalar type on valid values");
  {
    const std::string text =
      "s = hello\n"
      "t = true\nf = false\n"
      "i0 = 0\nim0 = -0\ni7 = 007\nimin = -2147483648\nimax = 2147483647\n"
      "u0 = 0\numax = 4294967295\n"
      "lmin = -9223372036854775808\nlmax = 9223372036854775807\nl53 = 9007199254740993\n"
      "d1 = 1.\nd2 = .5\nd3 = 1.e-2\nd4 = 2.194e-9\nd5 = -0\nd6 = 1e308\nd7 = 1E+2\n"
      "d8 = 4.9406564584124654e-324\nd9 = 1.7976931348623157e308\nd10 = 0e-999\n"
      "ca = a\ncs = \" \"\nch = \"#\"\ncq = '\nct = \\t\n";
    const ParameterFile p(file_with("scalars", text));
    check(p.get<std::string>("s") == "hello", "std::string");
    check(p.get<bool>("t") && !p.get<bool>("f"), "bool");
    check(p.get<int>("i0") == 0 && p.get<int>("im0") == 0 && p.get<int>("i7") == 7, "int: 0, -0, 007");
    check(p.get<int>("imin") == std::numeric_limits<int>::min() &&
          p.get<int>("imax") == std::numeric_limits<int>::max(), "int: the two limits");
    check(p.get<unsigned>("u0") == 0 && p.get<unsigned>("umax") == std::numeric_limits<unsigned>::max(),
          "unsigned: 0 and the max");
    check(p.get<std::int64_t>("lmin") == std::numeric_limits<std::int64_t>::min() &&
          p.get<std::int64_t>("lmax") == std::numeric_limits<std::int64_t>::max(),
          "int64_t: the two limits");
    check(p.get<std::int64_t>("l53") == 9007199254740993, "int64_t: 2^53 + 1 exactly");
    check(same_bits(p.get<double>("d1"), 1.) && same_bits(p.get<double>("d2"), .5) &&
          same_bits(p.get<double>("d3"), 1.e-2) && same_bits(p.get<double>("d4"), 2.194e-9) &&
          same_bits(p.get<double>("d5"), -0.) && same_bits(p.get<double>("d6"), 1e308) &&
          same_bits(p.get<double>("d7"), 100.), "double, bit for bit");
    check(same_bits(p.get<double>("d8"), 4.9406564584124654e-324) &&
          same_bits(p.get<double>("d9"), 1.7976931348623157e308) &&
          same_bits(p.get<double>("d10"), 0.),
          "double: the smallest subnormal, the largest double, a zero with a large exponent");
    check(p.get<char>("ca") == 'a' && p.get<char>("cs") == ' ' && p.get<char>("ch") == '#' &&
          p.get<char>("cq") == '\'' && p.get<char>("ct") == '\t', "char");
  }

  group("each vector type: several elements, one, none");
  {
    const std::string text =
      "s = a \"b c\"\nb = true false\ni = -1 2\nu = 1 2 3\nl = 9223372036854775807 -5\n"
      "d = 1.5 -2e3\nc = a \\t \" \"\n"
      "one = 7\nnone =\n";
    const ParameterFile p(file_with("vectors", text));
    check(p.get<Strings>("s") == Strings{"a", "b c"}, "std::vector<std::string>");
    check(p.get<std::vector<bool>>("b") == std::vector<bool>{true, false}, "std::vector<bool>");
    check(p.get<std::vector<int>>("i") == std::vector<int>{-1, 2}, "std::vector<int>");
    check(p.get<std::vector<unsigned>>("u") == std::vector<unsigned>{1, 2, 3}, "std::vector<unsigned>");
    check(p.get<std::vector<std::int64_t>>("l") ==
          std::vector<std::int64_t>{9223372036854775807, -5}, "std::vector<std::int64_t>");
    check(p.get<std::vector<double>>("d") == std::vector<double>{1.5, -2e3}, "std::vector<double>");
    check(p.get<std::vector<char>>("c") == std::vector<char>{'a', '\t', ' '}, "std::vector<char>");
    check(p.get<std::vector<int>>("one") == std::vector<int>{7} &&
          p.get<std::vector<double>>("one") == std::vector<double>{7.}, "one element");
    check(p.get<std::vector<int>>("none").empty() && p.get<Strings>("none").empty() &&
          p.get<std::vector<bool>>("none").empty(), "an empty list is an empty vector");
  }

  group("the refused values, with the whole message");
  {
    for (const char* v : {"True", "1", "yes", "\"\""})
      check_read_refused<bool>(v, std::string("'") + (v[0] == '"' ? "" : v) + "' is not true or false",
                               std::string("bool ") + v);

    for (const char* v : {"+5", "1e3", "5.0", "5.", "0x10", "1_000", "-"})
      check_read_refused<int>(v, std::string("'") + v + "' is not an integer", std::string("int ") + v);
    check_read_refused<int>("\"\"", "'' is not an integer", "int, an empty element");
    check_read_refused<unsigned>("+5", "'+5' is not an integer", "unsigned +5");

    const std::string intRange = " is out of the range of int (-2147483648 to 2147483647)";
    const std::string unsignedRange = " is out of the range of unsigned (0 to 4294967295)";
    const std::string int64Range =
      " is out of the range of int64 (-9223372036854775808 to 9223372036854775807)";
    check_read_refused<int>("2147483648", "'2147483648'" + intRange, "int, above the max");
    check_read_refused<int>("-2147483649", "'-2147483649'" + intRange, "int, below the min");
    check_read_refused<unsigned>("4294967296", "'4294967296'" + unsignedRange, "unsigned, above the max");
    check_read_refused<unsigned>("-1", "'-1'" + unsignedRange, "unsigned, negative");
    check_read_refused<unsigned>("-0", "'-0'" + unsignedRange, "unsigned, -0");
    check_read_refused<std::int64_t>("9223372036854775808", "'9223372036854775808'" + int64Range,
                                     "int64_t, above the max");
    check_read_refused<std::int64_t>("-9223372036854775809", "'-9223372036854775809'" + int64Range,
                                     "int64_t, below the min");

    for (const char* v : {"+1", "1e", "e3", "0x1p3", "1,5", "-", ".", "1.5.2", "1e+", "--1", "+inf"})
      check_read_refused<double>(v, std::string("'") + v + "' is not a number", std::string("double ") + v);
    check_read_refused<double>("\"\"", "'' is not a number", "double, an empty element");
    for (const char* v : {"nan", "inf", "infinity", "NaN", "INF", "Infinity", "-inf", "-nan"})
      check_read_refused<double>(v, std::string("'") + v + "' is not a finite number",
                                 std::string("double ") + v);
    for (const char* v : {"1e999", "-1e999", "1.7976931348623159e308", "1e-400", "-1e-400",
                          "2.4703282292062327e-324"})
      check_read_refused<double>(v, std::string("'") + v + "' is out of the range of double",
                                 std::string("double ") + v);

    const std::string notChar = " is not one character (write \\t for a tab)";
    check_read_refused<char>("ab", "'ab'" + notChar, "char ab");
    check_read_refused<char>("\"\"", "''" + notChar, "char, an empty element");
    check_read_refused<char>("\xC3\xA9", "'\xC3\xA9'" + notChar, "char, a non-ASCII letter");
    check_read_refused<char>("\\\\", "'\\\\'" + notChar, "char, two backslashes");

    check_read_refused<int>("", "no value", "a scalar read of an empty list");
    check_read_refused<std::string>("", "no value", "a std::string read of an empty list");
    check_read_refused<int>("1 2", "expected one value, found 2", "a scalar read of two elements");
    check_read_refused<std::string>("a b c", "expected one value, found 3; quote a value that holds spaces",
                                    "a std::string read of three elements");

    check_read_refused<std::vector<int>>("{0, 1, 2}", "element 1, '{0,', is not an integer",
                                         "an int list in braces");
    check_read_refused<std::vector<unsigned>>("1 -1", "element 2, '-1'," + unsignedRange,
                                              "the second element of an unsigned list");
    check_read_refused<std::vector<bool>>("true 1", "element 2, '1', is not true or false",
                                          "the second element of a bool list");
    check_read_refused<std::vector<double>>("1 nan", "element 2, 'nan', is not a finite number",
                                            "the second element of a double list");
    check_read_refused<std::vector<char>>("a bc", "element 2, 'bc'," + notChar,
                                          "the second element of a char list");

    const std::string path = file_with("missing_key", "a = 1\n");
    const ParameterFile p(path);
    check_message([&] { p.get<int>("threshold"); }, path + ": key 'threshold' is missing",
                  "a missing key");
    check_message([&] { p.get<Strings>("threshold"); }, path + ": key 'threshold' is missing",
                  "a missing key, as a list");
  }

  group("numbers do not depend on the global locale");
  {
    bool tried = false;
    for (const char* name : {"de_DE.UTF-8", "it_IT.UTF-8"}) {
      if (std::setlocale(LC_ALL, name) == nullptr) continue;
      bool cpp = true;
      try {
        std::locale::global(std::locale(name));
      }
      catch (const std::runtime_error&) {
        cpp = false;
        std::setlocale(LC_ALL, name);
      }
      if (std::strcmp(std::localeconv()->decimal_point, ",") != 0) {
        std::locale::global(std::locale::classic());
        std::setlocale(LC_ALL, "C");
        continue;
      }
      tried = true;
      std::cout << "    (the C locale " << name
                << (cpp ? " and the C++ one" : "; the C++ library has no named locale") << ")" << std::endl;
      const std::string path = file_with("locale", "a = 1.5\nb = 1,5\n");
      const ParameterFile p(path);
      check(p.get<double>("a") == 1.5, "1.5 is read with a decimal point");
      check_message([&] { p.get<double>("b"); }, path + ":2: key 'b': '1,5' is not a number",
                    "and 1,5 is refused");
      std::locale::global(std::locale::classic());
      std::setlocale(LC_ALL, "C");
      break;
    }
    if (!tried)
      std::cout << "    (neither de_DE.UTF-8 nor it_IT.UTF-8 exists with a decimal comma; not checked)"
                << std::endl;
  }

  group("the read marks");
  {
    const std::string path = file_with("marks", "a = 1\nb = 2\nc = x\nd = 4\ne = 5\n");
    const ParameterFile p(path);
    check(p.has("a") && p.has("b") && p.has("c"), "has finds the keys");
    check(p.unread().size() == 5, "has marks nothing");
    p.get<int>("a");
    check(p.get<int>("b", 7) == 2, "a default read of a present key gives its value");
    check(p.get<int>("z", 7) == 7, "and of a missing key the default");
    check(p.get<std::string>("missing", "none") == "none", "a std::string default");
    check(p.get<std::vector<int>>("missing", {0, 1, 2}) == std::vector<int>{0, 1, 2}, "a list default");
    check_throws([&] { p.get<int>("c"); }, "a read that throws");
    check(p.get<int>("d") == 4 && p.get<std::string>("d") == "4", "one key read as int and as std::string");
    check(same_entries(p.unread(), {{"e", 5}}), "get marks, the default overload marks, a throwing get marks");
  }

  group("unread() lists the keys in file order with their lines");
  {
    const ParameterFile p(file_with("order", "# header\nzeta = 1\n\nalpha = 2\nmid = 3\nbeta = 4\nomega = 5\n"));
    p.get<int>("omega");
    p.get<int>("alpha");
    check(same_entries(p.unread(), {{"zeta", 2}, {"mid", 5}, {"beta", 6}}), "file order, not name order");
    check(same_entries(p.unread(), {{"zeta", 2}, {"mid", 5}, {"beta", 6}}), "unread changes nothing");
  }

  group("a default never hides a broken value");
  {
    const std::string path = file_with("default_broken", "k = x\nj =\n");
    const ParameterFile p(path);
    check_message([&] { p.get<int>("k", 3); }, path + ":1: key 'k': 'x' is not an integer",
                  "a value that is not a T is an error, not the default");
    check_message([&] { p.get<int>("j", 3); }, path + ":2: key 'j': no value",
                  "an empty list read as a scalar with a default is an error");
  }

  group("nothing is printed");
  {
    std::ostringstream out, log, err;
    std::streambuf* const outBuffer = std::cout.rdbuf(out.rdbuf());
    std::streambuf* const logBuffer = std::clog.rdbuf(log.rdbuf());
    std::streambuf* const errBuffer = std::cerr.rdbuf(err.rdbuf());
    const ParameterFile p(file_with("quiet", "a = 1\n"));
    const double value = p.get<double>("absent", 2.5);
    p.unread();
    std::cout.rdbuf(outBuffer);
    std::clog.rdbuf(logBuffer);
    std::cerr.rdbuf(errBuffer);
    check(value == 2.5, "the default is returned");
    check(out.str().empty() && log.str().empty() && err.str().empty(),
          "a default read of a missing key prints nothing");
  }

  group("a const ParameterFile& can be read, and its reads count");
  {
    const ParameterFile file(file_with("const", "a = 1\nb = 2\n"));
    auto readA = [] (const ParameterFile& p) { return p.get<int>("a"); };
    check(readA(file) == 1, "read through a const reference");
    check(same_entries(file.unread(), {{"b", 2}}), "the read shows in unread()");
  }

  return report("test_parameter_file");
}
