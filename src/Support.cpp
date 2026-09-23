/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
 ********************************************************************/

/**
 *  @file src/Support.cpp
 *
 *  @brief Random numbers, input validation and coordinate conversion.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#include "internal.h"

namespace {

  std::uint64_t splitmix64 (std::uint64_t x)
  {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
  }

}


// ============================================================================


unsigned otswap::internal::uniform_int (std::mt19937& rng, const unsigned bound)
{
  if (bound == 0) return 0;

  const std::uint64_t range = bound;
  std::uint64_t product = (std::uint64_t)rng() * range;
  std::uint32_t low = (std::uint32_t)product;

  if (low < range) {
    const std::uint32_t threshold = (std::uint32_t)(0u - (std::uint32_t)range) % (std::uint32_t)range;
    while (low < threshold) {
      product = (std::uint64_t)rng() * range;
      low = (std::uint32_t)product;
    }
  }

  return (unsigned)(product >> 32);
}


// ============================================================================


double otswap::internal::uniform_real (std::mt19937& rng, const double lo, const double hi)
{
  const std::uint64_t high = (std::uint64_t)(rng() >> 5);
  const std::uint64_t low = (std::uint64_t)(rng() >> 6);
  const double unit = (double)((high << 26) + low) / 9007199254740992.;
  return lo + unit * (hi - lo);
}


// ============================================================================


std::mt19937 otswap::internal::stream (const unsigned seed, const unsigned realization,
                                      const Stream id)
{
  const std::uint64_t inner =
    splitmix64(((std::uint64_t)(std::uint32_t)seed << 32) | (std::uint32_t)realization);
  const std::uint64_t key = splitmix64(inner ^ (std::uint64_t)id);

  std::seed_seq words {(std::uint32_t)key, (std::uint32_t)(key >> 32)};
  return std::mt19937(words);
}


// ============================================================================


std::size_t otswap::internal::check_coordinates (const std::vector<double>& a,
                                                 const std::string& name,
                                                 const bool allowEmpty)
{
  if (a.empty()) {
    if (allowEmpty) return 0;
    throw Error(name + " is empty");
  }

  if (a.size() % 3 != 0)
    throw Error(name + " holds " + std::to_string(a.size()) +
                " entries, which is not a multiple of three");

  for (std::size_t i = 0; i < a.size(); ++i)
    if (!std::isfinite(a[i]))
      throw Error(name + " holds a non-finite value at entry " +
                  std::to_string(i) + " (object " + std::to_string(i/3) +
                  ", component " + std::to_string(i%3) + ")");

  return a.size() / 3;
}


// ============================================================================


unsigned otswap::internal::check_config (const Config& config)
{
  if (config.nRealizations == 0)
    throw Error("config.nRealizations is zero; at least one realization is required");

  if (!std::isfinite(config.convergence) || config.convergence <= 0. || config.convergence >= 1.)
    throw Error("config.convergence is " + std::to_string(config.convergence) +
                "; it must lie strictly between 0 and 1");

  if (!std::isfinite(config.cellSize) || config.cellSize <= 0.)
    throw Error("config.cellSize is " + std::to_string(config.cellSize) +
                "; it must be positive");

  if (config.seed != 0) return config.seed;

  std::random_device rd;
  unsigned drawn = rd();
  while (drawn == 0) drawn = rd();
  return drawn;
}


// ============================================================================


void otswap::internal::check_random_supply (const std::size_t nRandoms,
                                            const std::size_t nObjects,
                                            const unsigned nRealizations)
{
  const std::size_t needed = (std::size_t)nRealizations * nObjects;
  if (nRandoms >= needed) return;

  const double ratio = nObjects > 0 ? (double)nRandoms / (double)nObjects : 0.;
  throw Error("the randoms are not dense enough: " + std::to_string(nRandoms) +
              " randoms for " + std::to_string(nObjects) + " tracers is a ratio of " +
              std::to_string(ratio) + ", below the " + std::to_string(nRealizations) +
              " realizations requested; " + std::to_string(needed) + " randoms are needed");
}


// ============================================================================


void otswap::internal::to_cartesian (const double ra, const double dec, const double distance,
                                     double& x, double& y, double& z)
{
  x = distance * std::cos(dec) * std::cos(ra);
  y = distance * std::cos(dec) * std::sin(ra);
  z = distance * std::sin(dec);
}


// ============================================================================


void otswap::internal::to_sky (const double x, const double y, const double z,
                               double& ra, double& dec, double& distance)
{
  distance = std::sqrt(x*x + y*y + z*z);

  if (distance == 0.) {
    ra = 0.;
    dec = 0.;
    return;
  }

  ra = normalize_ra(std::atan2(y, x));
  dec = std::asin(std::max(-1., std::min(1., z/distance)));
}


// ============================================================================


double otswap::internal::normalize_ra (const double ra)
{
  constexpr double twopi = 2. * 3.14159265358979323846;
  double normalized = std::fmod(ra, twopi);
  if (normalized < 0.) normalized += twopi;
  return normalized;
}


// ============================================================================


double otswap::internal::profile_at (const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     const double at)
{
  if (x.size() == 1) return y[0];

  std::size_t i = (std::size_t)(std::lower_bound(x.begin(), x.end(), at) - x.begin());
  if (i == 0) i = 1;
  if (i == x.size()) i = x.size() - 1;
  if (at == x[i]) return y[i];

  const double f = (at - x[i-1]) / (x[i] - x[i-1]);
  return y[i-1] + f * (y[i] - y[i-1]);
}
