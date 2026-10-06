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
   *  @param positions Cartesian, 3 * N; every position finite and away
   *  from the origin.
   *  @param displacement 3 * N, for example Result::meanDisplacement. NaN
   *  is allowed, and gives NaN for that object; infinities are not.
   *
   *  @exception Error if the sizes differ or are not multiples of three,
   *  if a position is not finite or is at the origin, or if a displacement
   *  is infinite; the message names the object.
   */
  std::vector<double> lineOfSightProjection (const std::vector<double>& positions,
                                             const std::vector<double>& displacement);

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
  std::vector<double> lineOfSightProjection (const std::vector<double>& displacement,
                                             unsigned axis);

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
   *  @exception Error if the sizes differ, if a position is not finite, if
   *  the value of a valid object is not finite, or if sigma is negative or
   *  not finite; the message names the object.
   *  @exception meshsearch::Error if sigma > 0 and the positions span no
   *  volume (a single object, or all of them in one plane or on one line):
   *  the neighbour grid cannot be built.
   */
  std::vector<double> neighbourAverage (const std::vector<double>& positions,
                                        const std::vector<double>& values,
                                        const std::vector<unsigned>& validRealizations,
                                        double sigma,
                                        bool weightByRealizations = false);

  /**
   *  @brief As above, with two diagnostics per object.
   *
   *  @param[out] nNeighbours N entries: the number of valid objects
   *  averaged, the object itself included when it is valid. With
   *  sigma = 0, 1 for a valid object and 0 otherwise.
   *  @param[out] nRealizations N entries: the sum of validRealizations
   *  over those objects.
   */
  std::vector<double> neighbourAverage (const std::vector<double>& positions,
                                        const std::vector<double>& values,
                                        const std::vector<unsigned>& validRealizations,
                                        double sigma,
                                        bool weightByRealizations,
                                        std::vector<unsigned>& nNeighbours,
                                        std::vector<unsigned>& nRealizations);

  /**
   *  @brief The factor f / (b + 3 f / 5) that turns the averaged
   *  line-of-sight displacement into the shift to real space, at each
   *  redshift.
   *
   *  f is distances.growthRateAt(z). b(z) is interpolated linearly between
   *  the nodes of the bias table, and extrapolated linearly beyond its ends
   *  from the first or last segment. When any redshift lies outside the
   *  table, one line is written to std::clog, giving how many and the
   *  table's range.
   *
   *  @param redshift observed redshifts, inside the distance table's range.
   *  @param biasRedshift redshifts of the bias nodes: at least two, finite,
   *  strictly increasing.
   *  @param bias bias at each node: same length, finite and positive.
   *
   *  @exception Error if the distance table carries no growth rate, if a
   *  redshift is not finite or lies outside the distance table, if the bias
   *  table is malformed, or if the extrapolated bias at a redshift is not
   *  positive; the message names the object.
   */
  std::vector<double> rsdFactor (const std::vector<double>& redshift,
                                 const DistanceTable& distances,
                                 const std::vector<double>& biasRedshift,
                                 const std::vector<double>& bias);

  /**
   *  @brief As above, writing nothing: the number of redshifts at which
   *  b(z) was extrapolated is returned in nExtrapolated instead.
   */
  std::vector<double> rsdFactor (const std::vector<double>& redshift,
                                 const DistanceTable& distances,
                                 const std::vector<double>& biasRedshift,
                                 const std::vector<double>& bias,
                                 std::size_t& nExtrapolated);

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
  std::vector<double> shiftAlongLineOfSight (const std::vector<double>& positions,
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
  std::vector<double> shiftAlongLineOfSight (const std::vector<double>& positions,
                                             const std::vector<double>& shift,
                                             unsigned axis);

  // ==========================================================================
  // Redshift-space correction: the whole chain
  // ==========================================================================

  /**
   *  @brief A catalogue moved to real space. Row i describes input object
   *  i.
   */
  struct RealSpaceCatalog {

    std::size_t nObjects = 0;

    /// Lightcone: right ascension and declination, copied from the input,
    /// and the corrected redshift. Box: the Cartesian position, corrected
    /// along the line-of-sight axis. Flat, [object][3]. NaN in all three
    /// for an object in uncorrected.
    std::vector<double> positions;

    /// The diagnostics of neighbourAverage: valid objects averaged, and
    /// the sum of their valid realizations. nObjects entries each.
    std::vector<unsigned> nNeighbours;
    std::vector<unsigned> nRealizationsAveraged;

    /// Objects left without a correction, in increasing order: those with
    /// no valid object within 3 sigma, themselves included (with sigma = 0,
    /// those with no valid realization), and those outside the
    /// reconstruction's redshift cut.
    std::vector<std::size_t> uncorrected;
  };

  /**
   *  @brief Move a lightcone catalogue from redshift space to real space.
   *
   *  The chain, for the tracers that result was reconstructed from:
   *  1. the projection of result.meanDisplacement on the line of sight of
   *     each tracer, at toCartesian(tracersSky, distances);
   *  2. its gaussian average of width sigma over the tracers with a valid
   *     realization (neighbourAverage); a tracer without one still receives
   *     the average around its own position;
   *  3. the factor f / (b + 3 f / 5) at each observed redshift (rsdFactor),
   *     which writes to std::clog when b(z) is extrapolated;
   *  4. each tracer keeps its right ascension and declination and moves to
   *     the redshift of comoving distance d(z_i) + s_i, where s_i is the
   *     factor times the average.
   *
   *  Tracers flagged in result.outsideRedshiftCut are left out: they are
   *  not converted, take no part in any average, are not corrected, have
   *  diagnostics 0, and are listed in uncorrected.
   *
   *  Nothing checks that result belongs to these tracers beyond their
   *  number: pass the arrays the reconstruction was run on, in the same
   *  order.
   *
   *  @param result the reconstruction of these tracers.
   *  @param tracersSky sky coordinates, 3 * N, radians, as observed.
   *  @param biasRedshift,bias the bias table, as in rsdFactor.
   *  @param sigma width of the average, in Mpc/h; 0 for none. 10 Mpc/h is a
   *  reasonable starting value; the best value depends on the sample, and
   *  should be checked in each analysis.
   *  @param weightByRealizations as in neighbourAverage.
   *
   *  @exception Error if result does not hold N objects or is malformed, if
   *  a sky coordinate is invalid or a redshift lies outside the distance
   *  table, if the table carries no growth rate, if the bias table is
   *  malformed or extrapolates to a bias that is not positive, or if a
   *  corrected comoving distance is not positive or lies outside the
   *  distance table; the message names the object, and in the last case
   *  suggests a table covering a wider redshift range.
   *  @exception meshsearch::Error as neighbourAverage.
   */
  RealSpaceCatalog realSpaceLightcone (const Result& result,
                                       const std::vector<double>& tracersSky,
                                       const DistanceTable& distances,
                                       const std::vector<double>& biasRedshift,
                                       const std::vector<double>& bias,
                                       double sigma,
                                       bool weightByRealizations = false);

  /**
   *  @brief As above, writing nothing: the number of tracers at which b(z)
   *  was extrapolated is returned in nExtrapolated instead.
   */
  RealSpaceCatalog realSpaceLightcone (const Result& result,
                                       const std::vector<double>& tracersSky,
                                       const DistanceTable& distances,
                                       const std::vector<double>& biasRedshift,
                                       const std::vector<double>& bias,
                                       double sigma,
                                       bool weightByRealizations,
                                       std::size_t& nExtrapolated);

  /**
   *  @brief Move a box catalogue from redshift space to real space, with
   *  the line of sight along a Cartesian axis.
   *
   *  The chain is that of realSpaceLightcone, with the projection on the
   *  axis, a single factor f / (b + 3 f / 5) at the box's redshift, and the
   *  shift applied along the axis. Positions are not wrapped.
   *
   *  The box is not periodic, here or in reconstructBox: nothing flows
   *  through its faces, so modes on the scale of the box itself cannot be
   *  reconstructed.
   *
   *  @param result the reconstruction of these tracers.
   *  @param tracers Cartesian, 3 * N, as observed.
   *  @param axis the line of sight: 0, 1 or 2.
   *  @param redshift the box's redshift, at which f is taken.
   *  @param bias the tracers' linear bias; finite and positive.
   *  @param sigma as in realSpaceLightcone.
   *
   *  @exception Error if result does not hold N objects or is malformed, if
   *  a position is not finite, if the axis is not 0, 1 or 2, or as
   *  rsdFactorBox.
   *  @exception meshsearch::Error as neighbourAverage.
   */
  RealSpaceCatalog realSpaceBox (const Result& result,
                                 const std::vector<double>& tracers,
                                 unsigned axis,
                                 double redshift,
                                 const DistanceTable& distances,
                                 double bias,
                                 double sigma,
                                 bool weightByRealizations = false);

  // ==========================================================================
  // Bias table
  // ==========================================================================

  namespace io {

    /// A b(z) table, in the form rsdFactor and realSpaceLightcone take.
    struct BiasTable {
      std::vector<double> redshift;   ///< strictly increasing, at least two
      std::vector<double> bias;       ///< same length, finite and positive
    };

    /**
     *  @brief Read a b(z) table.
     *
     *  The format follows the extension, as in io::read: .fits, .fit and
     *  .fits.gz are read as FITS, from the columns named REDSHIFT and BIAS
     *  (case insensitively); anything else as ASCII, from columns 0 and 1,
     *  with the given delimiter and comment character.
     *
     *  The table is checked as rsdFactor checks its arrays, so a table read
     *  without error is one rsdFactor accepts.
     *
     *  @exception Error if the file cannot be read (as io::read), if it
     *  holds fewer than two rows, if a redshift is not finite or not greater
     *  than the one before, or if a bias is not finite or not positive; the
     *  message names the file and the row.
     */
    BiasTable readBiasTable (const std::string& file,
                             char delimiter = ' ', char comment = '#');

  }

}

#endif
