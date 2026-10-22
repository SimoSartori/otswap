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
 *  @file include/otswap/ParameterFile.h
 *
 *  @brief Reading a parameter file of `key = value` lines, with typed
 *  reads and the list of the keys never read.
 *
 *  Not included by otswap/OT.h or otswap/io.h.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_PARAMETERFILE_H
#define OTSWAP_PARAMETERFILE_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "otswap/OT.h"

namespace otswap {

  namespace io {

    /// @cond
    namespace detail {

      template <typename T> struct NonDeduced { using type = T; };

      template <typename T> struct IsParameterScalar : std::integral_constant<bool,
        std::is_same<T, std::string>::value || std::is_same<T, bool>::value ||
        std::is_same<T, int>::value || std::is_same<T, unsigned>::value ||
        std::is_same<T, std::int64_t>::value || std::is_same<T, double>::value ||
        std::is_same<T, char>::value> {};

      template <typename T> struct IsParameterType : IsParameterScalar<T> {};

      template <typename T> struct IsParameterType<std::vector<T>> : IsParameterScalar<T> {};

    }
    /// @endcond

    /**
     *  @brief A parameter file: one `key = value` per line, read once by the
     *  constructor and typed at each read.
     *
     *  Syntax:
     *  - one `key = value` per line; a key holds letters, digits and `_`, does
     *    not start with a digit, and is case-sensitive; each key at most once;
     *  - outside double quotes, `#` opens a comment that runs to the end of the
     *    line, also after a value; blank and comment lines are skipped;
     *  - a value is a list of elements separated by spaces or tabs; double
     *    quotes group an element that holds spaces, tabs or `#`, and do not
     *    change its meaning (`"5"` is 5, `"true"` is true); there is no escape,
     *    so an element cannot hold a double quote; single quotes, braces and
     *    commas are ordinary characters;
     *  - `key =` gives an empty list, `key = ""` one empty element.
     *
     *  A scalar read requires exactly one element; a std::vector read takes
     *  every element. The types are std::string, bool (`true` or `false`), int,
     *  unsigned, std::int64_t (decimal integers), double (finite, with an
     *  optional exponent), char (one character, or `\t` for a tab), and
     *  std::vector of each.
     *
     *  Every error is an otswap::Error whose message starts with the file and,
     *  when there is one, the line and the key: `run.ini:7: key 'seed': ...`.
     *
     *  The marks behind unread() are the only state a read changes; reading
     *  one ParameterFile from several threads at once is not safe.
     */
    class ParameterFile {
    public:

      /// A key of the file and the line, from 1, that sets it.
      struct Entry {
        std::string key;
        std::size_t line;
      };

      /// Reads and checks the syntax of the whole file.
      /// @param file the path; it appears as given in every error message
      /// @throw otswap::Error if the file cannot be opened or read, or at the
      /// first line that breaks the syntax
      explicit ParameterFile (const std::string& file);

      /// True if the file sets the key. Does not mark it read.
      bool has (const std::string& key) const;

      /// The value of a key the file must set, as a T. Marks the key read.
      /// @throw otswap::Error if the key is missing or its value is not a T
      template <typename T>
      T get (const std::string& key) const
      {
        static_assert(detail::IsParameterType<T>::value,
                      "ParameterFile::get reads std::string, bool, int, unsigned, std::int64_t, "
                      "double, char, or a std::vector of one of them");
        return read<T>(key, nullptr);
      }

      /// The value of a key, as a T, or defaultValue if the file does not set
      /// it. Marks the key read if it is set. T is always given explicitly:
      /// get<double>("sigma", 1).
      /// @throw otswap::Error if the key is set and its value is not a T
      template <typename T>
      T get (const std::string& key, const typename detail::NonDeduced<T>::type& defaultValue) const
      {
        static_assert(detail::IsParameterType<T>::value,
                      "ParameterFile::get reads std::string, bool, int, unsigned, std::int64_t, "
                      "double, char, or a std::vector of one of them");
        return read<T>(key, &defaultValue);
      }

      /// The keys no get has read, in file order, each with its line.
      std::vector<Entry> unread () const;

    private:

      struct Value {
        std::vector<std::string> elements;
        std::size_t line;
      };

      // Defined, and instantiated for the 14 types, in ParameterFile.cpp.
      template <typename T>
      T read (const std::string& key, const T* defaultValue) const;

      std::string m_file;
      std::vector<std::string> m_keys;                       // file order
      std::unordered_map<std::string, std::size_t> m_index;  // key -> position in m_keys
      std::vector<Value> m_values;                           // same order as m_keys
      mutable std::vector<bool> m_read;
    };

  }

}

#endif
