// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// =======================================================================
// Example code: how to run the OT reconstruction with lightcone geometry
// =======================================================================

/*
  This example explains how to reconstruct the displacement field of a
  survey catalogue with lightcone geometry. The lightcone geometry takes
  observed coordinates (right ascension, declination, redshift), a mean
  particle separation that changes with redshift, which otswap measures
  from the tracers, and the survey's footprint, here a HEALPix mask.

  The reconstruction matches the tracers to randoms by optimal
  transport: in each realization every tracer is paired with one
  random, and its displacement points from the tracer (the Eulerian,
  observed position) to that random (the Lagrangian position).
  Independent realizations are then averaged.

  The randoms carry the survey's geometry and selection function, which
  otswap does not model: use the survey's own random catalogue.

  Usage: lightcone [data_dir]. The catalogues are read from data_dir, by
  default the examples' data folder; the output is written to the
  folder output/ next to this file, created if missing.
*/

#include <cmath>
#include <filesystem>
#include <initializer_list>
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
    // Table::values holds the rows one after the other: right ascension,
    // declination and redshift for each object, the flat layout every
    // otswap function takes. The catalogues give the angles in degrees,
    // and otswap takes radians: the factor below is pi/180 as a double,
    // the one numpy.deg2rad and otswap's Python interface multiply by. The
    // degrees as read are kept for the output
    const std::vector<double> tracersDeg = otswap::io::read(data + "/lightcone_tracers.dat", {"0", "1", "2"}).values;
    const std::vector<double> randomsDeg = otswap::io::read(data + "/lightcone_randoms.dat", {"0", "1", "2"}).values;

    const double degToRad = 3.14159265358979323846 / 180.;
    std::vector<double> tracersSky = tracersDeg, randomsSky = randomsDeg;
    for (std::vector<double>* sky : {&tracersSky, &randomsSky})
      for (std::size_t i = 0; i < sky->size(); i += 3) {
        (*sky)[i]   *= degToRad;
        (*sky)[i+1] *= degToRad;
      }

    const std::size_t nTracers = tracersSky.size() / 3;
    const std::size_t nRandoms = randomsSky.size() / 3;

    std::cout << "Tracers read: " << nTracers << std::endl;
    std::cout << "Randoms read: " << nRandoms << ", " << (double)nRandoms / (double)nTracers << " per tracer" << std::endl;


    // ---------------------------------------
    // ------------ Read the mask ------------
    // ---------------------------------------

    // A HEALPix map in RING or NESTED ordering, read from a FITS file: a
    // pixel is observed when its value is greater than 0. A map already in
    // memory can be passed to the constructor taking an array instead.
    // When a mask is available, give it to the reconstruction rather than
    // a sky area: otswap then takes the area from it, and leaves out the
    // objects on unobserved pixels
    const otswap::Mask mask(data + "/lightcone_mask.fits");

    std::cout << "Mask: NSIDE " << mask.nside() << ", " << mask.skyAreaDeg2() << " deg^2" << std::endl;


    // ----------------------------------------------------
    // ------------ Set the cosmological model ------------
    // ----------------------------------------------------

    // The comoving distances, from Omega_m, h, w0 and wa of a flat
    // cosmology, tabulated by default from z = 0 to 10. The table must
    // cover the redshifts of tracers and randoms
    const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);


    // ------------------------------------------------------------------
    // ------------ Set the parameters of the reconstruction ------------
    // ------------------------------------------------------------------

    otswap::Config config;

    // Independent realizations. Each one consumes as many randoms as there
    // are tracers, disjoint from those of the others, so the randoms kept
    // must number at least nRealizations times the tracers kept
    config.nRealizations = 8;

    // The swap loop stops when a sweep changes fewer than this fraction of
    // the pairs. Larger values are faster and less accurate; the default
    // is 1e-3
    config.convergence = 1.e-2;

    // A fixed seed makes the run reproducible; 0 draws a new seed at each
    // run. The number of threads follows OMP_NUM_THREADS, and with a fixed
    // seed the result is the same whatever the number of threads
    config.seed = 12345;

    // The grid cell, in units of the mean particle separation, affects the
    // speed only, never the result; the default is 4
    config.cellSize = 4.;

    // After the reconstruction, a displacement whose path on the sky (the
    // great-circle arc from the tracer to its random) crosses more than
    // maxUnobservedPixelsCrossed unobserved pixels of the mask is marked
    // invalid. The count is of pixels, so the same threshold means a
    // different angle at a different NSIDE. These are the defaults
    config.rejectCrossings = true;
    config.maxUnobservedPixelsCrossed = 0;

    // With verbose, the default, the reconstruction writes to std::clog
    // how many tracers and randoms the redshift cut and the mask left out,
    // and how many displacements the crossing filter rejected. false
    // silences it; the counts stay in Result::selection either way
    config.verbose = true;

    // The mean particle separation is measured in nBins redshift bins of
    // equal width. Each bin needs at least 278 tracers: if one is refused,
    // the message says how many bins the catalogue supports
    const unsigned nBins = 30;

    // Tracers and randoms outside this redshift range are left out of the
    // reconstruction. The range here spans the whole catalogue: narrow it
    // to reconstruct a slice. The default, {}, cuts nothing
    const otswap::RedshiftCut cut {0.885, 1.10};


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // Tracers left out by the cut or the mask keep their row in the
    // result, flagged in outsideRedshiftCut or outsideMask, with NaN
    // displacements
    const otswap::Result result = otswap::reconstructLightcone(tracersSky, randomsSky, mask, nBins,
                                                               distances, config, cut);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // valid holds, for each realization and tracer, whether the
    // displacement survived the crossing filter; validRealizations counts
    // them per tracer, and meanDisplacement averages the valid ones, NaN
    // where there are none
    for (unsigned rec = 0; rec < result.nRealizations; ++rec) {
      std::size_t n = 0;
      for (std::size_t i = 0; i < nTracers; ++i) n += result.valid[rec * nTracers + i];
      std::cout << "Realization " << rec << ": " << n << " valid displacements of " << nTracers << std::endl;
    }

    std::vector<std::size_t> byCount(result.nRealizations + 1, 0);
    std::size_t leftOut = 0, noMean = 0;
    double length = 0.;
    for (std::size_t i = 0; i < nTracers; ++i) {
      ++byCount[result.validRealizations[i]];
      leftOut += (result.outsideRedshiftCut[i] || result.outsideMask[i]) ? 1 : 0;
      const double* d = &result.meanDisplacement[3*i];
      if (std::isnan(d[0])) ++noMean;
      else length += std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    }
    std::cout << "Tracers by number of valid realizations:";
    for (std::size_t n = 0; n < byCount.size(); ++n)
      std::cout << (n ? ", " : " ") << n << ": " << byCount[n];
    std::cout << std::endl;
    std::cout << "Tracers left out (redshift cut or mask): " << leftOut << std::endl;
    std::cout << "Tracers with a NaN mean displacement: " << noMean << std::endl;
    std::cout << "Mean length of the mean displacement: " << length / (double)(nTracers - noMean) << " Mpc/h" << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // toCartesian is the conversion the reconstruction applies, so these
    // are the positions the displacements start from
    const std::vector<double> tracers = otswap::toCartesian(tracersSky, distances);

    // Result::lagrangianSky holds the sky coordinates of each tracer's
    // mean Lagrangian position, in radians; they are written in degrees,
    // with a right ascension that rounds up to 360 taken back to 0
    const double radToDeg = 180. / 3.14159265358979323846;
    auto degrees = [&] (const double angle, const bool rightAscension) {
      const double d = angle * radToDeg;
      return (rightAscension && d >= 360.) ? 0. : d;
    };

    // One row per tracer, in the order of the input. A tracer left out, or
    // without a valid realization, has NaN Lagrangian coordinates and
    // displacement
    const std::vector<otswap::io::Column> columns = {
      {"tracRA", 'D', "tracer right ascension and declination, in degrees", {}},
      {"tracDec", 'D', "", {}},
      {"tracRed", 'D', "tracer redshift", {}},
      {"lagrRA", 'D', "mean Lagrangian position on the sky, in degrees, and its redshift", {}},
      {"lagrDec", 'D', "", {}},
      {"lagrRed", 'D', "", {}},
      {"tracX", 'D', "tracer position, in Mpc/h", {}},
      {"tracY", 'D', "", {}},
      {"tracZ", 'D', "", {}},
      {"lagrX", 'D', "mean Lagrangian position, tracer + mean displacement, in Mpc/h", {}},
      {"lagrY", 'D', "", {}},
      {"lagrZ", 'D', "", {}},
      {"displX", 'D', "mean displacement over the valid realizations, in Mpc/h", {}},
      {"displY", 'D', "", {}},
      {"displZ", 'D', "", {}},
      {"nValidRec", 'J', "number of valid realizations of the tracer", {}},
      {"outsideRedshiftCut", 'J', "1 if the tracer was left out by the redshift cut", {}},
      {"outsideMask", 'J', "1 if the tracer was left out by the mask", {}}};

    otswap::io::write(output + "/displacement_lightcone.dat", columns, nTracers,
                      [&] (const std::size_t i, std::vector<double>& row) {
                        for (int c = 0; c < 3; ++c) {
                          row[c]    = tracersDeg[3*i+c];
                          row[c+3]  = c < 2 ? degrees(result.lagrangianSky[3*i+c], c == 0)
                                            : result.lagrangianSky[3*i+c];
                          row[c+6]  = tracers[3*i+c];
                          row[c+9]  = tracers[3*i+c] + result.meanDisplacement[3*i+c];
                          row[c+12] = result.meanDisplacement[3*i+c];
                        }
                        row[15] = result.validRealizations[i];
                        row[16] = result.outsideRedshiftCut[i];
                        row[17] = result.outsideMask[i];
                      });

    std::cout << "Written: output/displacement_lightcone.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
