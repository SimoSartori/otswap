// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// =================================================================
// Example code: how to run the OT reconstruction with box geometry
// =================================================================

/*
  This example explains how to reconstruct the displacement field of a
  tracer catalogue with box geometry. The box geometry can be used with
  any catalogue in comoving Cartesian coordinates, including non-cubic
  ones, provided that the mean particle separation is approximately
  constant throughout the catalogue.

  The reconstruction matches the tracers to randoms by optimal
  transport: in each realization every tracer is paired with one
  random, and its displacement points from the tracer (the Eulerian,
  observed position) to that random (the Lagrangian position).
  Independent realizations are then averaged.

  The box is not periodic: nothing flows through its faces, so modes on
  the scale of the box itself are not reconstructed, and the
  displacements are least reliable near the faces. A box cut from a
  periodic simulation is treated as a plain sub-volume.

  Usage: box [data_dir]. The catalogues are read from data_dir, by
  default the examples' data folder; the output is written to the
  folder output/ next to this file, created if missing.
*/

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "otswap/OT.h"


// ==========================================================================


int main (int argc, char** argv)
{
  try {

    // -------------------------------------------------------------
    // ------------ Locate the input and output folders ------------
    // -------------------------------------------------------------

    // OTSWAP_EXAMPLES_DATA and OTSWAP_EXAMPLE_OUTPUT are set by this
    // example's CMakeLists.txt to the data folder and to output/ next to
    // this file. A first command-line argument replaces the data folder
    const std::string data = (argc > 1) ? argv[1] : OTSWAP_EXAMPLES_DATA;
    const std::string output = OTSWAP_EXAMPLE_OUTPUT;
    std::filesystem::create_directories(output);


    // ---------------------------------------------------------------
    // ------------ Read the tracer and random catalogues ------------
    // ---------------------------------------------------------------

    // io::read reads the listed columns of a table: 0-based indices for
    // ASCII files, column names for FITS files (.fits, .fit, .fits.gz).
    // Table::values holds the rows one after the other, x, y, z for each
    // object: the flat layout every otswap function takes. Coordinates
    // are comoving, in Mpc/h
    const std::vector<double> tracers = otswap::io::read(data + "/box_halos_real_space.dat", {"0", "1", "2"}).values;
    const std::vector<double> randoms = otswap::io::read(data + "/box_randoms.dat", {"0", "1", "2"}).values;

    const std::size_t nTracers = tracers.size() / 3;
    const std::size_t nRandoms = randoms.size() / 3;

    std::cout << "Tracers: " << nTracers << std::endl;
    std::cout << "Randoms: " << nRandoms << ", " << (double)nRandoms / (double)nTracers << " per tracer" << std::endl;


    // --------------------------------------------------------------
    // ------------ Compute the mean particle separation ------------
    // --------------------------------------------------------------

    // In box geometry the mean particle separation (mps) is a single
    // number given by the caller: here (V/N)^(1/3), with V the volume of
    // the tracers' bounding box and N their number. It sets the scale of
    // the search for swap partners, so it must describe the whole
    // catalogue: with a strongly varying density use the lightcone
    // geometry, which measures it as a function of redshift
    double volume = 1.;
    for (int c = 0; c < 3; ++c) {
      double lo = tracers[c], hi = tracers[c];
      for (std::size_t i = 0; i < nTracers; ++i) {
        lo = std::min(lo, tracers[3*i+c]);
        hi = std::max(hi, tracers[3*i+c]);
      }
      volume *= hi - lo;
    }
    const double mps = std::pow(volume / (double)nTracers, 1. / 3.);

    std::cout << "Mean particle separation: " << mps << " Mpc/h" << std::endl;


    // ------------------------------------------------------------------
    // ------------ Set the parameters of the reconstruction ------------
    // ------------------------------------------------------------------

    otswap::Config config;

    // Independent realizations. Each one consumes N randoms, disjoint
    // from those of the others, so the random catalogue must hold at least
    // nRealizations x N objects: the catalogue read here has exactly 8 per
    // tracer. An empty random array lets otswap draw the randoms itself,
    // uniformly in the tracers' bounding box
    config.nRealizations = 8;

    // The swap loop stops when a sweep changes fewer than this fraction of
    // the pairs. Larger values are faster and less accurate; the default
    // is 1e-3
    config.convergence = 1.e-2;

    // A fixed seed makes the run reproducible; 0 draws a new seed at each
    // run. The number of threads follows OMP_NUM_THREADS, and with a fixed
    // seed the result is the same whatever the number of threads
    config.seed = 12345;

    // The grid cell, in units of mps, affects the speed only, never the
    // result; the default is 4
    config.cellSize = 4.;


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    const otswap::Result result = otswap::reconstructBox(tracers, randoms, mps, config);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // In a box result every realization is valid for every tracer, and
    // meanDisplacement is the average over all of them. Its average over
    // the tracers is the offset between the centres of mass of the
    // randoms used and of the tracers
    std::size_t allValid = 0;
    for (std::size_t i = 0; i < nTracers; ++i)
      allValid += result.validRealizations[i] == result.nRealizations ? 1 : 0;
    std::cout << "Realizations: " << result.nRealizations << ", all valid for " << allValid << " tracers" << std::endl;

    double length = 0., average[3] = {0., 0., 0.};
    for (std::size_t i = 0; i < nTracers; ++i) {
      const double* d = &result.meanDisplacement[3*i];
      length += std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      for (int c = 0; c < 3; ++c) average[c] += d[c];
    }
    std::cout << "Mean length of the mean displacement: " << length / (double)nTracers << " Mpc/h" << std::endl;
    std::cout << "Average of the mean displacement over the tracers: (" << average[0] / (double)nTracers << ", "
              << average[1] / (double)nTracers << ", " << average[2] / (double)nTracers << ") Mpc/h" << std::endl;
    std::cout << "Displacement of the first tracer in the first realization: (" << result.displacement[0] << ", "
              << result.displacement[1] << ", " << result.displacement[2] << ") Mpc/h" << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // One row per tracer, in the order of the input: its position, its
    // mean Lagrangian position (the position plus the mean displacement),
    // and the mean displacement. io::write generates the rows on demand,
    // so nothing the size of the catalogue is held twice; the format
    // follows the extension, ASCII here
    const std::vector<otswap::io::Column> columns = {
      {"tracX", 'D', "tracer position, in Mpc/h", {}},
      {"tracY", 'D', "", {}},
      {"tracZ", 'D', "", {}},
      {"lagrX", 'D', "mean Lagrangian position, tracer + mean displacement, in Mpc/h", {}},
      {"lagrY", 'D', "", {}},
      {"lagrZ", 'D', "", {}},
      {"displX", 'D', "mean displacement over the valid realizations, in Mpc/h", {}},
      {"displY", 'D', "", {}},
      {"displZ", 'D', "", {}}};

    otswap::io::write(output + "/displacement_box.dat", columns, nTracers,
                      [&] (const std::size_t i, std::vector<double>& row) {
                        for (int c = 0; c < 3; ++c) {
                          row[c]   = tracers[3*i+c];
                          row[c+3] = tracers[3*i+c] + result.meanDisplacement[3*i+c];
                          row[c+6] = result.meanDisplacement[3*i+c];
                        }
                      });

    std::cout << "Written: output/displacement_box.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
