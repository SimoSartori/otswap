// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// ======================================================================
// Example code: how to run the OT reconstruction in a periodic box
// ======================================================================

/*
  This example explains how to reconstruct the displacement field of
  halos in a cube with periodic boundary conditions, such as a
  simulation box. The box is given by its limits; every distance is then
  the minimum image, so no face is a boundary: a tracer near a face may
  be matched to a random near the opposite one, as anywhere else in the
  box.

  Each side of a periodic box must be an integer multiple of the
  shortest one, so a cube, or a cuboid such as 500 x 500 x 1000, is
  accepted. Tracers outside the limits are wrapped into the box, by a
  whole number of periods per axis, and reported; the displacements
  point along the shortest way, within half a period on each axis.

  Usage: periodic_box [data_dir]. The catalogue is read from data_dir, by
  default the examples' data folder; the output is written to the folder
  output/ next to this file, created if missing.
*/

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


    // ---------------------------------------------------
    // ------------ Read the tracer catalogue ------------
    // ---------------------------------------------------

    // Halos of a simulation in a cube of side 1000 Mpc/h, [0, 1000) on each
    // axis, comoving coordinates in Mpc/h, x, y, z for each object
    const std::vector<double> tracers = otswap::io::read(data + "/periodic_halos.dat", {"0", "1", "2"}).values;
    const std::size_t nTracers = tracers.size() / 3;

    std::cout << "Tracers: " << nTracers << std::endl;


    // ------------------------------------------------------------------
    // ------------ Set the parameters of the reconstruction ------------
    // ------------------------------------------------------------------

    otswap::Config config;

    // Independent realizations, each on its own set of randoms. With no
    // randoms given, as here, otswap draws them itself, uniformly in the
    // box, a different set for each realization; 8 realizations keep the
    // example fast
    config.nRealizations = 8;

    // The sweeps stop once the fraction of successful swaps in a sweep no
    // longer exceeds this threshold; 1e-2 keeps the example fast. For an
    // analysis, start from the defaults, 30 realizations and 1e-3: the
    // guide's section "Choosing the parameters" advises on both
    config.convergence = 1.e-2;

    // A fixed seed makes the run reproducible, whatever the number of threads
    config.seed = 12345;

    // The grid cell, in units of mps, affects the speed only, never the
    // result. In a periodic box the cells divide the shortest side a whole
    // number of times, so the cell used can be smaller: Result::cellSide
    config.cellSize = 3.;

    config.verbosity = otswap::Verbosity::Normal;
    config.log = &std::cout;


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // The box is periodic: the limits are required, and the mean particle
    // separation, (V/N)^(1/3) with V the volume of the limits, comes from
    // them. A tracer outside the limits would be wrapped into the box and
    // counted in Result::selection; here none is
    otswap::BoxLimits box;
    box.lower = {0., 0., 0.};
    box.upper = {1000., 1000., 1000.};

    const otswap::Result result = otswap::reconstructBox(tracers, {}, config, box, true);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    std::cout << "Periodic: " << (result.periodic ? "yes" : "no") << ", tracers wrapped: "
              << result.selection.tracersWrapped << std::endl;
    std::cout << "Mean particle separation: " << result.mps << " Mpc/h, grid cell side: " << result.cellSide
              << " Mpc/h" << std::endl;

    // Each displacement component lies within half a period, [-500, 500)
    // Mpc/h here, however far apart the tracer and its random lie in the
    // plain coordinates
    double length = 0., largest = 0.;
    for (std::size_t i = 0; i < nTracers; ++i) {
      const double* d = &result.meanDisplacement[3*i];
      length += std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      for (int c = 0; c < 3; ++c) largest = std::fmax(largest, std::fabs(d[c]));
    }
    std::cout << "Mean length of the mean displacement: " << length / (double)nTracers << " Mpc/h" << std::endl;
    std::cout << "Largest component of a mean displacement: " << largest << " Mpc/h" << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // One row per tracer, in the order of the input: its position, its
    // mean Lagrangian position, wrapped into the box, the mean
    // displacement and the number of valid realizations, with a header
    // describing the columns (lengths in coordinate units), the box, the
    // periodicity and the parameters of the run
    otswap::io::writeDisplacementMean(output + "/displacement_periodic_box.dat", result);

    std::cout << "Written: output/displacement_periodic_box.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
