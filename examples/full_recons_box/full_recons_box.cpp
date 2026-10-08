// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// =====================================================================
// Example code: how to move a box catalogue from redshift to real space
// =====================================================================

/*
  This example explains how to correct a catalogue in box geometry for
  redshift-space distortions with the OT reconstruction. The tracers
  are given in redshift space, distorted along one Cartesian axis, the
  line of sight. The example runs the reconstruction of the box example
  on them, then moves each tracer along that axis to the position it
  would have in real space.

  The correction of tracer i is a shift along the axis,

    s_i = f(z) / (b + 3 f(z) / 5) * <Psi_axis>_i ,

  with f the linear growth rate at the box's redshift z, b the tracers'
  linear bias, and <Psi_axis>_i the component along the axis of the
  reconstructed mean displacement, averaged with a gaussian of width
  sigma over the tracers around i.

  The box is not periodic, here or in the reconstruction: nothing flows
  through its faces, so modes on the scale of the box itself are not
  reconstructed, and the corrected positions are not wrapped. A box cut
  from a periodic simulation is treated as a plain sub-volume.

  Usage: full_recons_box [data_dir]. The catalogues are read from
  data_dir, by default the examples' data folder; the output is written
  to the folder output/ next to this file, created if missing.
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

    // Comoving Cartesian coordinates in Mpc/h, x, y, z for each object,
    // one row after the other. The tracers are in redshift space, with
    // the line of sight along y; the randoms are uniform in the volume
    const std::vector<double> tracers = otswap::io::read(data + "/box_halos_redshift_space.dat", {"0", "1", "2"}).values;
    const std::vector<double> randoms = otswap::io::read(data + "/box_randoms.dat", {"0", "1", "2"}).values;

    const std::size_t nTracers = tracers.size() / 3;
    const std::size_t nRandoms = randoms.size() / 3;

    std::cout << "Tracers: " << nTracers << std::endl;
    std::cout << "Randoms: " << nRandoms << ", " << (double)nRandoms / (double)nTracers << " per tracer" << std::endl;


    // --------------------------------------------------------------
    // ------------ Compute the mean particle separation ------------
    // --------------------------------------------------------------

    // (V/N)^(1/3), with V the volume of the tracers' bounding box, as in
    // the box example
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


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // As in the box example: 8 realizations, a convergence threshold of
    // 1e-2 and a fixed seed
    otswap::Config config;
    config.nRealizations = 8;
    config.convergence = 1.e-2;
    config.seed = 12345;

    const otswap::Result result = otswap::reconstructBox(tracers, randoms, mps, config);

    std::cout << "Realizations: " << result.nRealizations << std::endl;


    // ----------------------------------------------------
    // ------------ Set the cosmological model ------------
    // ----------------------------------------------------

    // Only the linear growth rate f at the box's redshift is needed. The
    // table computes it from Omega_m, h, w0 and wa of a flat cosmology
    const otswap::DistanceTable distances(0.3186, 0.67, -1., 0.);


    // ----------------------------------------------------------------
    // ------------ Correct the redshift-space distortions ------------
    // ----------------------------------------------------------------

    // The box's redshift, the tracers' linear bias, and the line of sight:
    // axis 0, 1 or 2 for x, y or z
    const double redshift = 1.;
    const double bias = 1.91255;
    const unsigned axis = 1;

    // sigma, the width of the gaussian average in Mpc/h: 10 Mpc/h is a
    // reasonable starting value, but the best value depends on the sample
    // and should be checked in each analysis. With
    // weightByRealizations, each neighbour would be weighted by its number
    // of valid realizations.
    //
    // Pass the positions the reconstruction was run on, in the same order.
    // A tracer with no valid neighbour within 3 sigma, which in a box
    // happens only for an isolated one, is left uncorrected, with NaN
    // coordinates, and listed in uncorrected
    const double sigma = 10.;
    const otswap::RealSpaceCatalog catalogue = otswap::realSpaceBox(result, tracers, axis, redshift, distances,
                                                                    bias, sigma, false);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    std::cout << "RSD factor f/(b + 3f/5): " << otswap::rsdFactorBox(redshift, distances, bias) << std::endl;

    std::size_t corrected = 0;
    double sum = 0., squares = 0.;
    for (std::size_t i = 0; i < nTracers; ++i) {
      if (std::isnan(catalogue.positions[3*i+axis])) continue;
      ++corrected;
      const double shift = catalogue.positions[3*i+axis] - tracers[3*i+axis];
      sum += shift;
      squares += shift * shift;
    }
    std::cout << "Tracers corrected: " << corrected << std::endl;
    std::cout << "Tracers left uncorrected: " << catalogue.uncorrected.size() << std::endl;
    std::cout << "Shift along the line of sight: mean " << sum / (double)corrected << ", rms "
              << std::sqrt(squares / (double)corrected) << " Mpc/h" << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // One row per tracer, in the order of the input: the corrected
    // position, NaN for an uncorrected tracer, then the number of valid
    // realizations and the two diagnostics of the average
    const std::vector<otswap::io::Column> columns = {
      {"X", 'D', "corrected position, in Mpc/h", {}},
      {"Y", 'D', "corrected position, in Mpc/h", {}},
      {"Z", 'D', "corrected position, in Mpc/h", {}},
      {"nValidRec", 'J', "number of valid OT realizations of the tracer", {}},
      {"nNeighbours", 'J', "number of tracers with a valid OT realization averaged within 3 sigma, the tracer included if valid", {}},
      {"nRealizationsAveraged", 'J', "sum of nValidRec over those tracers", {}}};

    otswap::io::write(output + "/reconstructed_catalogue_box.dat", columns, nTracers,
                      [&] (const std::size_t i, std::vector<double>& row) {
                        for (int c = 0; c < 3; ++c) row[c] = catalogue.positions[3*i+c];
                        row[3] = result.validRealizations[i];
                        row[4] = catalogue.nNeighbours[i];
                        row[5] = catalogue.nRealizationsAveraged[i];
                      });

    std::cout << "Written: output/reconstructed_catalogue_box.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
