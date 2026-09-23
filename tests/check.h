/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
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

  void check_close (const double a, const double b, const double tolerance,
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
