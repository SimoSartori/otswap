/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
 ********************************************************************/

/**
 *  @file examples/example_box.cpp
 *
 *  @brief Reconstructing the displacement field of a cubic box.
 */

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "otswap/OT.h"

int main ()
{
  // A cubic box of tracers. Coordinates are flat and row-major: one
  // object after the other, x then y then z. Nothing is nested, so a
  // numpy array of shape (N, 3) maps onto this without a copy.
  const unsigned nSide = 10;
  const double spacing = 8.;

  std::vector<double> tracers;
  tracers.reserve(3ull * nSide * nSide * nSide);

  std::mt19937 rng(1);
  std::uniform_real_distribution<double> jitter(-2., 2.);

  for (unsigned i = 0; i < nSide; ++i)
    for (unsigned j = 0; j < nSide; ++j)
      for (unsigned k = 0; k < nSide; ++k) {
        tracers.push_back(i * spacing + jitter(rng));
        tracers.push_back(j * spacing + jitter(rng));
        tracers.push_back(k * spacing + jitter(rng));
      }

  const std::size_t nObjects = tracers.size() / 3;

  // The mean particle separation of the box. In box geometry it is a
  // single number the caller supplies; the lightcone entry points measure
  // it from the tracers instead.
  const double mps = spacing;

  otswap::Config config;
  config.nRealizations = 4;     // four independent matchings, then averaged
  config.convergence = 1.e-3;   // a sweep changing fewer pairs than this ends it
  config.seed = 12345;          // a fixed seed makes the run reproducible
  config.cellSize = 4.;         // grid cells, in units of mps: speed only

  // Passing an empty random array asks for the randoms to be drawn
  // uniformly in the bounding box of the tracers, one disjoint set per
  // realization. Supply your own array instead when the randoms have to
  // carry a selection function.
  const otswap::Result result = otswap::reconstructBox(tracers, {}, mps, config);

  std::cout << "objects:      " << result.nObjects << "\n"
            << "realizations: " << result.nRealizations << "\n";

  // The Lagrangian position of object i is its coordinate plus its mean
  // displacement.
  double sum = 0., largest = 0.;
  for (std::size_t i = 0; i < nObjects; ++i) {
    double squared = 0.;
    for (int c = 0; c < 3; ++c)
      squared += result.meanDisplacement[3*i+c] * result.meanDisplacement[3*i+c];
    const double length = std::sqrt(squared);
    sum += length;
    if (length > largest) largest = length;
  }

  std::cout << "mean |displacement|:    " << sum / (double)nObjects << "\n"
            << "largest |displacement|: " << largest << "\n";

  // Writing the field out. The table layer picks FITS or ASCII from the
  // extension, and generates the rows on demand so that nothing the size
  // of the catalog is held twice.
  std::vector<otswap::io::Column> columns;
  for (const char* name : {"tracX", "tracY", "tracZ",
                           "lagrX", "lagrY", "lagrZ",
                           "displX", "displY", "displZ"})
    columns.push_back({name, 'D', "", {}});

  otswap::io::write("displacement_box.dat", columns, nObjects,
                    [&] (const std::size_t i, std::vector<double>& row) {
                      for (int c = 0; c < 3; ++c) {
                        row[c]   = tracers[3*i+c];
                        row[c+3] = tracers[3*i+c] + result.meanDisplacement[3*i+c];
                        row[c+6] = result.meanDisplacement[3*i+c];
                      }
                    });

  std::cout << "written to displacement_box.dat" << std::endl;

  return 0;
}
