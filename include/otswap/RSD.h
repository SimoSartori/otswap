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
 *  @file include/otswap/RSD.h
 *
 *  @brief The redshift-space correction: moving a catalogue from redshift
 *  space to real space with the displacement field of a reconstruction.
 *
 *  The correction of object i is a shift along its line of sight,
 *
 *    s_i = f(z_i) / (b(z_i) + 3 f(z_i) / 5)  <Psi . r_hat>_i ,
 *
 *  where Psi is the reconstructed displacement (Result::meanDisplacement,
 *  from the observed to the reconstructed position), r_hat the line of
 *  sight, f the linear growth rate, b the linear bias of the tracers, and
 *  < > a gaussian average over the neighbouring objects. realSpaceLightcone
 *  and realSpaceBox run the whole chain; its four steps are also available
 *  on their own.
 *
 *  Included by otswap/OT.h; either header may be included first.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_RSD_H
#define OTSWAP_RSD_H

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "otswap/OT.h"

namespace otswap {

  // ==========================================================================
  // Redshift-space correction: the four steps
  // ==========================================================================

  /**
   *  @brief Component of each displacement along the line of sight of its
   *  object, r_i . d_i / |r_i|: positive when the displacement points away
   *  from the observer.
   *
   *  @param displacement 3 * N, for example Result::meanDisplacement. NaN
   *  is allowed, and gives NaN for that object; infinities are not.
   *  @param positions Cartesian, 3 * N, for example Result::tracers; every
   *  position finite and away from the origin.
   *
   *  @exception Error if the sizes differ or are not multiples of three,
   *  if a position is not finite or is at the origin, or if a displacement
   *  is infinite; the message names the object.
   */
  std::vector<double> radialProjection (const std::vector<double>& displacement,
                                        const std::vector<double>& positions);

  /**
   *  @brief Component of each displacement along a Cartesian axis, the
   *  line of sight of a box.
   *
   *  @param displacement 3 * N; NaN allowed, infinities not.
   *  @param axis 0, 1 or 2.
   *
   *  @exception Error if the size is not a multiple of three, if the axis
   *  is not 0, 1 or 2, or if a displacement is infinite.
   */
  std::vector<double> axisProjection (const std::vector<double>& displacement,
                                      unsigned axis);

  /// The result of neighbourAverage: N entries each.
  struct NeighbourAverage {

    /// The average at each object; NaN with no valid object within 3 sigma.
    std::vector<double> values;

    /// The number of valid objects averaged, the object itself included
    /// when it is valid. With sigma = 0, 1 for a valid object and 0
    /// otherwise.
    std::vector<unsigned> nNeighbours;

    /// The sum of validRealizations over those objects.
    std::vector<unsigned> nRealizationsAveraged;
  };

  /**
   *  @brief Gaussian average of a value over the neighbours of each object.
   *
   *  An object is valid when its validRealizations is greater than zero.
   *  Each object, valid or not, receives the weighted mean of the values
   *  of the valid objects within 3 sigma of its own position, its own value
   *  included when it is valid. The weight of a neighbour at distance d is
   *  exp(-d^2 / (2 sigma^2)), multiplied by its validRealizations when
   *  weightByRealizations is set, and the weights are normalised on their
   *  sum, so a constant value is returned unchanged, at the edges of the
   *  catalogue too. An object with no valid object within 3 sigma receives
   *  NaN.
   *
   *  With sigma = 0 no average is taken and values is returned unchanged.
   *
   *  The neighbours are found on a grid built inside, over the bounding box
   *  of the positions given; nothing passed in is modified. The result does
   *  not depend on that grid within rounding, and is the same bit for bit
   *  for the same inputs, whatever the number of threads.
   *
   *  @param positions Cartesian, 3 * N, finite.
   *  @param values N entries; finite wherever validRealizations > 0, and
   *  not read elsewhere.
   *  @param validRealizations N entries.
   *  @param sigma width of the gaussian, in the unit of the positions;
   *  finite and non-negative. 10 Mpc/h is a reasonable starting value; the
   *  best value depends on the sample, and should be checked in each
   *  analysis.
   *  @param weightByRealizations weight each neighbour by its count of
   *  valid realizations as well.
   *
   *  @return the average and two diagnostics per object.
   *
   *  @exception Error if the sizes differ, if a position is not finite, if
   *  the value of a valid object is not finite, or if sigma is negative or
   *  not finite; the message names the object.
   *  @exception meshsearch::Error if sigma > 0 and the positions span no
   *  volume (a single object, or all of them in one plane or on one line):
   *  the neighbour grid cannot be built.
   */
  NeighbourAverage neighbourAverage (const std::vector<double>& positions,
                                     const std::vector<double>& values,
                                     const std::vector<unsigned>& validRealizations,
                                     double sigma,
                                     bool weightByRealizations = false);

  /**
   *  @brief A b(z) table: the linear bias of the tracers at redshift nodes,
   *  interpolated linearly between them and extrapolated linearly beyond
   *  them from the first or last segment. Read one from a file with
   *  io::readBiasTable.
   */
  struct BiasTable {
    std::vector<double> redshift;   ///< at least two, finite, strictly increasing
    std::vector<double> bias;       ///< same length, finite and positive
  };

  /**
   *  @brief The factor f / (b + 3 f / 5) that turns the averaged
   *  line-of-sight displacement into the shift to real space, at each
   *  redshift.
   *
   *  f is distances.growthRateAt(z), b(z) the bias table at z. Nothing is
   *  written: the redshifts at which b(z) is extrapolated beyond the table
   *  are counted in *nExtrapolated, when given.
   *
   *  @param redshift observed redshifts, inside the distance table's range.
   *  @param bias the b(z) table.
   *  @param nExtrapolated if not null, receives the number of redshifts
   *  outside the bias table.
   *
   *  @exception Error if the distance table carries no growth rate, if a
   *  redshift is not finite or lies outside the distance table, if the bias
   *  table is malformed, or if the extrapolated bias at a redshift is not
   *  positive; the message names the object.
   */
  std::vector<double> rsdFactor (const std::vector<double>& redshift,
                                 const DistanceTable& distances,
                                 const BiasTable& bias,
                                 std::size_t* nExtrapolated = nullptr);

  /**
   *  @brief The factor f / (b + 3 f / 5) of a box at a single redshift,
   *  f = distances.growthRateAt(redshift), with a constant bias.
   *
   *  @param bias finite and positive.
   *
   *  @exception Error if the redshift is not finite or lies outside the
   *  distance table, if the table carries no growth rate, or if the bias is
   *  not positive.
   */
  double rsdFactorBox (double redshift, const DistanceTable& distances, double bias);

  /**
   *  @brief Move each position by shift_i along its own line of sight:
   *  r_i (1 + shift_i / |r_i|). A positive shift moves the object away from
   *  the observer.
   *
   *  @param positions Cartesian, 3 * N; every position finite and away
   *  from the origin.
   *  @param shift N entries; NaN gives NaN in the whole row, infinities are
   *  refused.
   *
   *  @exception Error if the sizes differ, if a position is not finite or
   *  is at the origin, or if a shift is infinite; the message names the
   *  object.
   */
  std::vector<double> shiftRadially (const std::vector<double>& positions,
                                     const std::vector<double>& shift);

  /**
   *  @brief Move each position by shift_i along a Cartesian axis. Positions
   *  are not wrapped: a box is not treated as periodic.
   *
   *  @param positions Cartesian, 3 * N, finite.
   *  @param shift N entries; NaN gives NaN along the axis, infinities are
   *  refused.
   *  @param axis 0, 1 or 2.
   */
  std::vector<double> shiftAlongAxis (const std::vector<double>& positions,
                                      const std::vector<double>& shift,
                                      unsigned axis);

  // ==========================================================================
  // Redshift-space correction: the whole chain
  // ==========================================================================

  /// Options of realSpaceLightcone and realSpaceBox.
  struct CorrectionConfig {

    /// Weight each neighbour by its count of valid realizations as well,
    /// as in neighbourAverage.
    bool weightByRealizations = false;

    /// What the call reports; see Verbosity.
    Verbosity verbosity = Verbosity::Normal;

    /// Where the report is written; null writes nothing, as Silent.
    std::ostream* log = &std::clog;
  };

  /// What the correction did with each tracer.
  enum class CorrectionStatus : std::uint8_t {
    Corrected         = 0,  ///< had a valid realization; moved
    MovedByNeighbours = 1,  ///< had no valid realization; moved with its neighbours' average
    NoValidNeighbour  = 2,  ///< no valid tracer within 3 sigma (with sigma = 0: no valid realization); not moved
    LeftOut           = 3   ///< outside the reconstruction's redshift cut or mask; not moved
  };

  /**
   *  @brief A catalogue moved to real space. Row i describes input object
   *  i.
   */
  struct RealSpaceCatalog {

    std::size_t nObjects = 0;

    /// The geometry of the reconstruction corrected.
    Geometry geometry = Geometry::Box;

    /// Lightcone: right ascension, folded into [0, 2 pi), and declination of
    /// the tracer, and its corrected redshift; radians. Flat, [object][3]. NaN
    /// rows for the tracers not moved. Empty in a box.
    std::vector<double> sky;

    /// The corrected Cartesian position, Mpc/h. Box: the tracer with its
    /// axis coordinate moved. Lightcone: r (1 + shift / |r|), the tracer r
    /// moved along its line of sight. Flat, [object][xyz]. NaN rows for the
    /// tracers not moved.
    std::vector<double> cartesian;

    /// The shift applied along the line of sight, in Mpc/h: the factor times
    /// the averaged projection. NaN for a tracer not moved. nObjects entries.
    std::vector<double> shift;

    /// The factor f / (b + 3 f / 5) at each tracer's observed redshift; in a
    /// box, the single factor repeated. NaN for a tracer left out of the
    /// reconstruction. nObjects entries.
    std::vector<double> factor;

    /// What the correction did with each tracer. nObjects entries.
    std::vector<CorrectionStatus> status;

    /// The valid realizations of each tracer, copied from the result
    /// corrected. nObjects entries.
    std::vector<unsigned> validRealizations;

    /// The diagnostics of neighbourAverage: valid tracers averaged, and the
    /// sum of their valid realizations; 0 for a tracer left out. nObjects
    /// entries each.
    std::vector<unsigned> nNeighbours;
    std::vector<unsigned> nRealizationsAveraged;

    /// Tracers at whose redshift b(z) was extrapolated beyond its table; 0
    /// in a box.
    std::size_t nExtrapolated = 0;

    /// The width of the average and its weighting, as given.
    double sigma = std::numeric_limits<double>::quiet_NaN();
    bool weightByRealizations = false;

    /// Box: the line-of-sight axis, the redshift and the bias, as given. 0
    /// and NaN in a lightcone.
    unsigned axis = 0;
    double boxRedshift = std::numeric_limits<double>::quiet_NaN();
    double boxBias = std::numeric_limits<double>::quiet_NaN();

    /// Wall time of the call that made the catalogue, in seconds.
    double elapsedSeconds = 0.;
  };

  /**
   *  @brief Move a lightcone catalogue from redshift space to real space.
   *
   *  The chain, for the tracers of the result:
   *  1. the projection of result.meanDisplacement on the line of sight of
   *     each tracer, at result.tracers (radialProjection);
   *  2. its gaussian average of width sigma over the tracers with a valid
   *     realization (neighbourAverage); a tracer without one still receives
   *     the average around its own position;
   *  3. the factor f / (b + 3 f / 5) at each observed redshift, as
   *     rsdFactor computes it; the tracers at which b(z) is extrapolated are
   *     counted in RealSpaceCatalog::nExtrapolated, and reported;
   *  4. each tracer keeps its right ascension and declination and moves to
   *     the redshift of comoving distance d(z_i) + s_i, where s_i is the
   *     factor times the average; its Cartesian position moves by s_i along
   *     its line of sight.
   *
   *  Tracers flagged in result.outsideRedshiftCut or result.outsideMask are
   *  left out: they take no part in any average, are not moved, have
   *  diagnostics 0, and the status LeftOut.
   *
   *  @param result a lightcone reconstruction, with its tracers and
   *  tracersSky, as reconstructLightcone returns it.
   *  @param bias the b(z) table.
   *  @param sigma width of the average, in Mpc/h; 0 for none. 10 Mpc/h is a
   *  reasonable starting value; the best value depends on the sample, and
   *  should be checked in each analysis.
   *  @param config the weighting, and what the call reports.
   *
   *  @exception Error if result is not a lightcone result, is malformed or
   *  lacks its tracers or tracersSky, if a tracer kept is not finite or is
   *  at the origin, if the table carries no growth rate, if the bias table
   *  is malformed or extrapolates to a bias that is not positive, or if a
   *  corrected comoving distance is not positive or lies outside the
   *  distance table; the message names the object, and in the last case
   *  suggests a table covering a wider redshift range.
   *  @exception meshsearch::Error as neighbourAverage.
   */
  RealSpaceCatalog realSpaceLightcone (const Result& result,
                                       const DistanceTable& distances,
                                       const BiasTable& bias,
                                       double sigma,
                                       const CorrectionConfig& config = {});

  /**
   *  @brief Move a box catalogue from redshift space to real space, with
   *  the line of sight along a Cartesian axis.
   *
   *  The chain is that of realSpaceLightcone, on result.tracers, with the
   *  projection on the axis, a single factor f / (b + 3 f / 5) at the box's
   *  redshift, and the shift applied along the axis. Positions are not
   *  wrapped.
   *
   *  The box is not periodic, here or in reconstructBox: nothing flows
   *  through its faces, so modes on the scale of the box itself cannot be
   *  reconstructed.
   *
   *  @param result a box reconstruction, with its tracers.
   *  @param axis the line of sight: 0, 1 or 2.
   *  @param redshift the box's redshift, at which f is taken.
   *  @param bias the tracers' linear bias; finite and positive.
   *  @param sigma as in realSpaceLightcone.
   *  @param config the weighting, and what the call reports.
   *
   *  @exception Error if result is not a box result, is malformed or lacks
   *  its tracers, if the axis is not 0, 1 or 2, or as rsdFactorBox.
   *  @exception meshsearch::Error as neighbourAverage.
   */
  RealSpaceCatalog realSpaceBox (const Result& result,
                                 unsigned axis,
                                 double redshift,
                                 const DistanceTable& distances,
                                 double bias,
                                 double sigma,
                                 const CorrectionConfig& config = {});

}

#endif
