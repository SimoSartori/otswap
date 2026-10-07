// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
//
// A program built against otswap the way a consumer builds it, from its
// public header alone: a mask from a map in memory, a small box
// reconstruction, and otswap::Error on invalid input. Returns non-zero on
// any mismatch.

#include <cmath>
#include <cstdio>
#include <vector>

#include "otswap/OT.h"

int main ()
{
  int failures = 0;
  auto check = [&failures] (const bool ok, const char* what) {
    if (!ok) {
      std::fprintf(stderr, "consumer: %s\n", what);
      ++failures;
    }
  };

  // NSIDE 1, RING: pixels 1 and 7 unobserved.
  const std::vector<double> values = {1., 0., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1.};
  const otswap::Mask mask(values, otswap::PixelOrdering::Ring);
  const double pi = 3.14159265358979323846;
  const double sphere = 4. * pi * (180. / pi) * (180. / pi);
  check(mask.nside() == 1, "the mask has NSIDE 1");
  check(std::fabs(mask.skyAreaDeg2() - sphere * 10. / 12.) < 1.e-9 * sphere,
        "the mask covers 10 of its 12 pixels");
  check(mask.allows(0.1, 1.5) && mask.allows(0., 0.) && !mask.allows(2.3, 1.) &&
        !mask.allows(4.712, 0.), "the mask allows the observed pixels only");

  // A jittered lattice of 10^3 tracers 8 apart.
  std::vector<double> tracers;
  unsigned state = 1;
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        for (const int c : {i, j, k}) {
          state = state * 1103515245u + 12345u;
          tracers.push_back(8. * c + 2. * ((double)(state >> 16) / 65536. - 0.5));
        }

  otswap::Config config;
  config.seed = 7;
  config.nRealizations = 2;
  const otswap::Result result = otswap::reconstructBox(tracers, {}, 8., config);
  check(result.nObjects == 1000 && result.nRealizations == 2, "one row per tracer and realization");
  bool finite = true;
  for (const double d : result.displacement) finite = finite && std::isfinite(d);
  check(finite && result.displacement.size() == 6000, "every displacement is finite");
  check(result.selection.message().empty() && result.outsideMask.size() == 1000,
        "a box result has no selection");

  bool threw = false;
  try {
    otswap::reconstructBox({1., 2.}, {}, 8., config);
  }
  catch (const otswap::Error&) {
    threw = true;
  }
  check(threw, "invalid input raises otswap::Error");

  if (failures == 0) std::printf("consumer: all checks passed\n");
  return failures == 0 ? 0 : 1;
}
