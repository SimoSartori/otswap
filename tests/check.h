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
 *  @file tests/check.h
 *
 *  @brief The assertion the test programs share. Each test prints the
 *  guarantee it checks and returns non-zero on the first failure.
 */

#ifndef OTSWAP_TESTS_CHECK_H
#define OTSWAP_TESTS_CHECK_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

  int failures = 0;

  void check (const bool condition, const std::string& what)
  {
    if (condition) return;
    std::cerr << "FAILED: " << what << std::endl;
    ++failures;
  }

  [[maybe_unused]] void check_close (const double a, const double b, const double tolerance,
                                     const std::string& what)
  {
    check(std::fabs(a-b) <= tolerance,
          what + " (got " + std::to_string(a) + ", expected " + std::to_string(b) +
          ", tolerance " + std::to_string(tolerance) + ")");
  }

  template <typename F>
  void check_throws (F&& f, const std::string& what)
  {
    try {
      f();
    }
    catch (const otswap::Error&) {
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

  void group (const std::string& what)
  {
    std::cout << "  " << what << std::endl;
  }

  int report (const std::string& name)
  {
    if (failures == 0) {
      std::cout << name << ": all groups passed" << std::endl;
      return 0;
    }
    std::cerr << name << ": " << failures << " failure(s)" << std::endl;
    return 1;
  }

}

#endif
