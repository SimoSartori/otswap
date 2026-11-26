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

  The example then keeps every realization in a file, with the
  catalogue's IDs, reads it back, and uses the realizations to check
  whether their number is enough.

  Without the periodic option, the box is not periodic: nothing flows
  through its faces, so modes on the scale of the box itself are not
  reconstructed, and the displacements are least reliable near the
  faces. A box cut from a periodic simulation is treated as a plain
  sub-volume; a whole simulation box with periodic boundary conditions
  is reconstructed with the periodic option: see the periodic_box
  example.

  Usage: box [data_dir]. The catalogues are read from data_dir, by
  default the examples' data folder; the output is written to the
  folder output/ next to this file, created if missing.
*/

#include <cmath>
#include <cstdint>
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

    // io::read reads the listed columns of a table: for ASCII files
    // 0-based indices, or the names of the "###" line of a file otswap
    // wrote; for FITS files (.fits, .fit, .fits.gz) the column names.
    // Table::values holds the rows one after the other, x, y, z for each
    // object: the flat layout every otswap function takes. Coordinates
    // are comoving, in Mpc/h
    const std::vector<double> tracers = otswap::io::read(data + "/box_halos_real_space.dat", {"0", "1", "2"}).values;
    const std::vector<double> randoms = otswap::io::read(data + "/box_randoms.dat", {"0", "1", "2"}).values;

    const std::size_t nTracers = tracers.size() / 3;
    const std::size_t nRandoms = randoms.size() / 3;

    std::cout << "Tracers: " << nTracers << std::endl;
    std::cout << "Randoms: " << nRandoms << ", " << (double)nRandoms / (double)nTracers << " per tracer" << std::endl;


    // ------------------------------------------------------------------
    // ------------ Set the parameters of the reconstruction ------------
    // ------------------------------------------------------------------

    otswap::Config config;

    // Independent realizations. Each one consumes N randoms, disjoint
    // from those of the others, so the random catalogue must hold at least
    // nRealizations x N objects: the catalogue read here has exactly 8 per
    // tracer (404,376 randoms for 50,547 halos), hence 8 realizations. An
    // empty random array lets otswap draw the randoms itself, uniformly in
    // the box limits, or without them in the bounding box of the tracers
    config.nRealizations = 8;

    // The sweeps stop once the fraction of successful swaps in a sweep,
    // swaps per tracer visited, no longer exceeds this threshold. Larger
    // values are faster and less accurate; 1e-2 keeps the example fast.
    // For an analysis, start from the defaults, 30 realizations and 1e-3:
    // the guide's section "Choosing the parameters" advises on both
    config.convergence = 1.e-2;

    // A fixed seed makes the run reproducible; 0 draws a new seed at each
    // run, and Result::config.seed records the one drawn. The number of
    // threads follows OMP_NUM_THREADS, and with a fixed seed the result is
    // the same whatever the number of threads
    config.seed = 12345;

    // The grid cell, in units of mps, affects the speed only, never the
    // result; the default is 3
    config.cellSize = 3.;

    // What the reconstruction reports, and where: with box limits, how many
    // tracers and randoms lie outside them, then the line of the call;
    // Detailed adds, before that line, the mean particle separation and its
    // source. The report goes to std::clog by default; here to std::cout,
    // in order with this program's own lines
    config.verbosity = otswap::Verbosity::Detailed;
    config.log = &std::cout;


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // The catalogues fill a known volume, the cube [269, 669] Mpc/h, so it
    // is given as the box limits, closed on every face. Tracers outside
    // them would keep their row, flagged in Result::outsideBox with NaN
    // displacements, and randoms outside them would be dropped; here none
    // is. Without limits the box is the bounding box of the tracers and the
    // randoms, as in the full_recons_box example.
    //
    // In box geometry the mean particle separation (mps) is a single
    // number. Without one, as here, otswap takes (V/N)^(1/3), with V the
    // volume of the limits and N the number of tracers inside them, and
    // records it in Result::mps; an overload takes it from the caller. It
    // sets the radius of the first guess, 8 mps, within which each tracer
    // is paired with a free random, and the cell of the neighbour grids, so
    // the density must be roughly uniform across the box
    otswap::BoxLimits box;
    box.lower = {269., 269., 269.};
    box.upper = {669., 669., 669.};

    otswap::Result result = otswap::reconstructBox(tracers, randoms, config, box);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // Every realization is valid for every tracer inside the limits, and
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
    // the mean displacement, the number of valid realizations and the flag
    // outsideBox, with a header describing the columns (lengths in
    // coordinate units), the box and the parameters of the run. The format
    // follows the extension, ASCII here; io::readDisplacementMean reads
    // this file back (see the lightcone example)
    otswap::io::writeDisplacementMean(output + "/displacement_box.dat", result);

    std::cout << "Written: output/displacement_box.dat" << std::endl;


    // --------------------------------------------------------------
    // ------------ Keep every realization, with the IDs ------------
    // --------------------------------------------------------------

    // A catalogue usually names its objects by an ID column, read exactly
    // by the overload of io::read that takes integer columns, and set as
    // Result::ids. These files have none, so here the position of each halo
    // among the data rows of its file, counted from 1, stands for it. The
    // reconstruction sets no IDs and no computation reads them: the writers
    // put them in a first column, id, and realSpaceBox and
    // realSpaceLightcone copy them into the corrected catalogue. An empty ids, or exactly 0..N-1, means none
    result.ids.resize(nTracers);
    for (std::size_t i = 0; i < nTracers; ++i) result.ids[i] = (std::int64_t)i + 1;

    // writeDisplacementComplete writes every realization of every tracer,
    // with its valid flag, the whole configuration and the box, exactly:
    // with 17 significant digits in ASCII, in binary in FITS, as here,
    // which is smaller and faster to read. readDisplacementComplete
    // rebuilds the Result, every double identical, the means recomputed
    // from the realizations
    otswap::io::writeDisplacementComplete(output + "/displacement_complete_box.fits", result);

    std::cout << "Written: output/displacement_complete_box.fits" << std::endl;

    const otswap::Result readBack = otswap::io::readDisplacementComplete(output + "/displacement_complete_box.fits");

    const bool identical = readBack.displacement == result.displacement && readBack.valid == result.valid &&
                           readBack.meanDisplacement == result.meanDisplacement && readBack.ids == result.ids;
    std::cout << "Read back: " << readBack.nRealizations << " realizations of " << readBack.nObjects
              << " tracers, IDs " << readBack.ids.front() << " to " << readBack.ids.back() << ", "
              << (identical ? "identical" : "NOT identical") << " to the result written" << std::endl;


    // ----------------------------------------------------------
    // ------------ Check the number of realizations ------------
    // ----------------------------------------------------------

    // recomputeMeans recomputes validRealizations, meanDisplacement and the
    // Lagrangian positions from displacement and valid, for instance after
    // a filter of one's own has cleared some entries of valid. Here it
    // keeps the first half of the realizations. The rms difference between
    // the mean over half of them and the mean over all of them shows how
    // much the mean still depends on their number: the larger it is
    // against the mean length printed above, the more realizations are
    // worth adding
    otswap::Result half = readBack;
    const std::size_t nHalf = half.nRealizations / 2;
    for (std::size_t rec = nHalf; rec < half.nRealizations; ++rec)
      for (std::size_t i = 0; i < nTracers; ++i) half.valid[rec * nTracers + i] = 0;
    otswap::recomputeMeans(half);

    double squares = 0.;
    for (std::size_t k = 0; k < 3 * nTracers; ++k) {
      const double difference = half.meanDisplacement[k] - readBack.meanDisplacement[k];
      squares += difference * difference;
    }
    std::cout << "rms difference between the mean displacements of " << nHalf << " and " << readBack.nRealizations
              << " realizations: " << std::sqrt(squares / (double)nTracers) << " Mpc/h" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
