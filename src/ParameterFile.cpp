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
 *  @file src/ParameterFile.cpp
 *
 *  @brief The parameter-file reader: the syntax check of each line and the
 *  typed reads.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#include "otswap/ParameterFile.h"

namespace {

  bool blank (const char c)
  {
    return c == ' ' || c == '\t';
  }

  bool digit (const char c)
  {
    return c >= '0' && c <= '9';
  }

  bool valid_key (const std::string& key)
  {
    auto letter = [] (const char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_'; };
    if (key.empty() || !letter(key[0])) return false;
    for (const char c : key)
      if (!letter(c) && !digit(c)) return false;
    return true;
  }

  /// The whole file, read in binary mode. A read that fails after the
  /// file opened, as on a directory, is an error of its own.
  std::string contents (const std::string& file)
  {
    std::unique_ptr<std::FILE, int (*)(std::FILE*)> input(std::fopen(file.c_str(), "rb"), &std::fclose);
    if (!input)
      throw otswap::Error(file + ": cannot open the parameter file");
    std::string text;
    char buffer[4096];
    std::size_t n;
    while ((n = std::fread(buffer, 1, sizeof buffer, input.get())) > 0)
      text.append(buffer, n);
    if (std::ferror(input.get()))
      throw otswap::Error(file + ": cannot read the parameter file");
    return text;
  }

  // ==========================================================================
  // Typing one element. Each returns the end of the message, after the
  // element, when the element is not of the type, and an empty string when
  // it is.
  // ==========================================================================

  std::string type_element (const std::string& element, std::string& value)
  {
    value = element;
    return "";
  }

  std::string type_element (const std::string& element, bool& value)
  {
    if (element == "true") value = true;
    else if (element == "false") value = false;
    else return "is not true or false";
    return "";
  }

  std::string type_element (const std::string& element, char& value)
  {
    if (element.size() == 1) value = element[0];
    else if (element == "\\t") value = '\t';
    else return "is not one character (write \\t for a tab)";
    return "";
  }

  template <typename I> const char* integer_name ();
  template <> const char* integer_name<int> () { return "int"; }
  template <> const char* integer_name<unsigned> () { return "unsigned"; }
  template <> const char* integer_name<std::int64_t> () { return "int64"; }

  /// A decimal integer: an optional '-' and one or more digits, nothing
  /// else, converted by std::from_chars, which needs no locale.
  template <typename I>
  std::string type_integer (const std::string& element, I& value)
  {
    const bool negative = !element.empty() && element[0] == '-';
    const std::size_t first = negative ? 1 : 0;
    bool integer = element.size() > first;
    for (std::size_t i = first; integer && i < element.size(); ++i) integer = digit(element[i]);
    if (!integer) return "is not an integer";

    const std::string range = std::string("is out of the range of ") + integer_name<I>() + " (" +
      std::to_string(std::numeric_limits<I>::min()) + " to " +
      std::to_string(std::numeric_limits<I>::max()) + ")";
    if (negative && !std::numeric_limits<I>::is_signed) return range;
    const auto conversion = std::from_chars(element.data(), element.data()+element.size(), value);
    if (conversion.ec == std::errc::result_out_of_range) return range;
    if (conversion.ec != std::errc() || conversion.ptr != element.data()+element.size())
      return "is not an integer";
    return "";
  }

  std::string type_element (const std::string& element, int& value) { return type_integer(element, value); }
  std::string type_element (const std::string& element, unsigned& value) { return type_integer(element, value); }
  std::string type_element (const std::string& element, std::int64_t& value) { return type_integer(element, value); }

  bool equal_ignoring_case (const std::string& a, const char* b)
  {
    std::size_t i = 0;
    for (; i < a.size() && b[i] != '\0'; ++i) {
      const char c = (a[i] >= 'A' && a[i] <= 'Z') ? char(a[i]-'A'+'a') : a[i];
      if (c != b[i]) return false;
    }
    return i == a.size() && b[i] == '\0';
  }

  /// A double: the grammar is checked here, so that what is accepted does
  /// not depend on the standard library, then the element is converted by
  /// a stream in the classic locale, which ignores the global one. The
  /// libraries flag a range error differently (some flag subnormals too, and
  /// store either infinity or the largest double on overflow), so the
  /// value is judged here: an overflow, or a zero from digits that are not
  /// all zero, is out of range; a subnormal is kept.
  std::string type_element (const std::string& element, double& value)
  {
    const std::size_t first = (!element.empty() && element[0] == '-') ? 1 : 0;
    const std::string unsigned_part = element.substr(first);
    if (equal_ignoring_case(unsigned_part, "nan") || equal_ignoring_case(unsigned_part, "inf") ||
        equal_ignoring_case(unsigned_part, "infinity"))
      return "is not a finite number";

    std::size_t i = first;
    std::size_t digits = 0;
    bool nonzero = false;
    for (; i < element.size() && digit(element[i]); ++i, ++digits) nonzero |= element[i] != '0';
    if (i < element.size() && element[i] == '.') {
      const std::size_t before = digits;
      for (++i; i < element.size() && digit(element[i]); ++i, ++digits) nonzero |= element[i] != '0';
      if (before == 0 && digits == 0) return "is not a number";
    }
    if (digits == 0) return "is not a number";
    if (i < element.size() && (element[i] == 'e' || element[i] == 'E')) {
      ++i;
      if (i < element.size() && (element[i] == '+' || element[i] == '-')) ++i;
      std::size_t exponent = 0;
      for (; i < element.size() && digit(element[i]); ++i) ++exponent;
      if (exponent == 0) return "is not a number";
    }
    if (i != element.size()) return "is not a number";

    std::istringstream input(element);
    input.imbue(std::locale::classic());
    double converted = 0;
    input >> converted;
    if (input.fail() && !(std::fabs(converted) < std::numeric_limits<double>::max()))
      return "is out of the range of double";
    if (converted == 0 && nonzero) return "is out of the range of double";
    value = converted;
    return "";
  }

  template <typename T> struct Scalar { using type = T; static constexpr bool vector = false; };
  template <typename T> struct Scalar<std::vector<T>> { using type = T; static constexpr bool vector = true; };

}


// ============================================================================


otswap::io::ParameterFile::ParameterFile (const std::string& file)
  : m_file(file)
{
  const std::string text = contents(file);

  std::size_t number = 0;
  for (std::size_t start = 0; start < text.size(); ) {
    std::size_t end = text.find('\n', start);
    if (end == std::string::npos) end = text.size();
    std::string line = text.substr(start, end-start);
    start = end+1;
    ++number;
    const std::string at = file + ":" + std::to_string(number) + ": ";

    if (number == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
    if (!line.empty() && line.back() == '\r') line.pop_back();

    // The key part, up to the first '=' or '#'.
    const std::size_t cut = line.find_first_of("=#");
    if (cut == std::string::npos || line[cut] == '#') {
      for (std::size_t i = 0; i < (cut == std::string::npos ? line.size() : cut); ++i)
        if (!blank(line[i]))
          throw otswap::Error(at + "expected 'key = value', a comment or a blank line");
      continue;
    }
    std::size_t keyStart = 0, keyEnd = cut;
    while (keyStart < keyEnd && blank(line[keyStart])) ++keyStart;
    while (keyEnd > keyStart && blank(line[keyEnd-1])) --keyEnd;
    const std::string key = line.substr(keyStart, keyEnd-keyStart);
    if (key.empty())
      throw otswap::Error(at + "no key before '='");
    const std::string atKey = at + "key '" + key + "': ";
    if (!valid_key(key))
      throw otswap::Error(atKey + "not a valid key; a key holds letters, digits and _, and does not "
                          "start with a digit");

    // The value part, everything after the first '=', cut into elements.
    std::vector<std::string> elements;
    std::size_t i = cut+1;
    while (true) {
      while (i < line.size() && blank(line[i])) ++i;
      if (i == line.size() || line[i] == '#') break;
      if (line[i] == '"') {
        const std::size_t close = line.find('"', i+1);
        if (close == std::string::npos)
          throw otswap::Error(atKey + "unterminated quote");
        elements.push_back(line.substr(i+1, close-i-1));
        i = close+1;
        if (i < line.size() && !blank(line[i]) && line[i] != '#')
          throw otswap::Error(atKey + "text glued to a quote; separate the elements with a space");
      }
      else {
        const std::size_t begin = i;
        for (; i < line.size() && !blank(line[i]) && line[i] != '#'; ++i)
          if (line[i] == '"')
            throw otswap::Error(atKey + "text glued to a quote; separate the elements with a space");
        elements.push_back(line.substr(begin, i-begin));
      }
    }

    const auto found = m_index.find(key);
    if (found != m_index.end())
      throw otswap::Error(atKey + "repeated; first set at line " +
                          std::to_string(m_values[found->second].line));
    m_index.emplace(key, m_keys.size());
    m_keys.push_back(key);
    m_values.push_back(Value{std::move(elements), number});
  }

  m_read.assign(m_keys.size(), false);
}


// ============================================================================


bool otswap::io::ParameterFile::has (const std::string& key) const
{
  return m_index.count(key) != 0;
}


// ============================================================================


std::vector<otswap::io::ParameterFile::Entry> otswap::io::ParameterFile::unread () const
{
  std::vector<Entry> entries;
  for (std::size_t k = 0; k < m_keys.size(); ++k)
    if (!m_read[k]) entries.push_back(Entry{m_keys[k], m_values[k].line});
  return entries;
}


// ============================================================================


template <typename T>
T otswap::io::ParameterFile::read (const std::string& key, const T* defaultValue) const
{
  const auto found = m_index.find(key);
  if (found == m_index.end()) {
    if (defaultValue != nullptr) return *defaultValue;
    throw otswap::Error(m_file + ": key '" + key + "' is missing");
  }
  m_read[found->second] = true;

  const Value& value = m_values[found->second];
  const std::string at = m_file + ":" + std::to_string(value.line) + ": key '" + key + "': ";
  using Element = typename Scalar<T>::type;

  if constexpr (Scalar<T>::vector) {
    T result;
    result.reserve(value.elements.size());
    for (std::size_t e = 0; e < value.elements.size(); ++e) {
      Element element{};
      const std::string failure = type_element(value.elements[e], element);
      if (!failure.empty())
        throw otswap::Error(at + "element " + std::to_string(e+1) + ", '" + value.elements[e] +
                            "', " + failure);
      result.push_back(element);
    }
    return result;
  }
  else {
    if (value.elements.empty())
      throw otswap::Error(at + "no value");
    if (value.elements.size() > 1)
      throw otswap::Error(at + "expected one value, found " + std::to_string(value.elements.size()) +
                          (std::is_same<T, std::string>::value ? "; quote a value that holds spaces" : ""));
    T result{};
    const std::string failure = type_element(value.elements[0], result);
    if (!failure.empty())
      throw otswap::Error(at + "'" + value.elements[0] + "' " + failure);
    return result;
  }
}

namespace otswap {
  namespace io {

    template std::string ParameterFile::read (const std::string&, const std::string*) const;
    template bool ParameterFile::read (const std::string&, const bool*) const;
    template int ParameterFile::read (const std::string&, const int*) const;
    template unsigned ParameterFile::read (const std::string&, const unsigned*) const;
    template std::int64_t ParameterFile::read (const std::string&, const std::int64_t*) const;
    template double ParameterFile::read (const std::string&, const double*) const;
    template char ParameterFile::read (const std::string&, const char*) const;
    template std::vector<std::string> ParameterFile::read (const std::string&, const std::vector<std::string>*) const;
    template std::vector<bool> ParameterFile::read (const std::string&, const std::vector<bool>*) const;
    template std::vector<int> ParameterFile::read (const std::string&, const std::vector<int>*) const;
    template std::vector<unsigned> ParameterFile::read (const std::string&, const std::vector<unsigned>*) const;
    template std::vector<std::int64_t> ParameterFile::read (const std::string&, const std::vector<std::int64_t>*) const;
    template std::vector<double> ParameterFile::read (const std::string&, const std::vector<double>*) const;
    template std::vector<char> ParameterFile::read (const std::string&, const std::vector<char>*) const;

  }
}
