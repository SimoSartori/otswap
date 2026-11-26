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
 *  @file include/otswap/OT.h
 *
 *  @brief Optimal transport reconstruction of the displacement field
 *  of a tracer catalog, by local swapping.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_OT_H
#define OTSWAP_OT_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace otswap {

  // ==========================================================================
  // Errors
  // ==========================================================================

  /// Base class of every exception otswap throws on invalid input or
  /// failure. std::bad_alloc is not converted, and propagates unchanged,
  /// from inside a parallel region too. An error of the neighbour search,
  /// a meshsearch::Error or an exception derived from it, reaches the
  /// caller as an Error whose message is "meshsearch: " followed by its
  /// own, with the original nested in it (std::throw_with_nested), so
  /// std::rethrow_if_nested recovers it.
  ///
  /// An Error records where it was constructed: the source file, the
  /// function and the line. The location is not part of what(), which
  /// holds the message alone. An error raised by a shared helper, or
  /// rethrown with context after a parallel region or a lookup, carries
  /// the location of that helper or of the rethrow.
  class Error : public std::runtime_error {
  public:

    /// @param what the message, returned by what()
    /// @param file,function,line the location; left to their defaults,
    /// they are those of the expression that constructs the Error
    explicit Error (const std::string& what,
                    const char* file = __builtin_FILE(),
                    const char* function = __builtin_FUNCTION(),
                    int line = __builtin_LINE())
      : std::runtime_error(what), m_file(file), m_function(function), m_line(line) {}

    /// Source file in which the Error was constructed, as the compiler
    /// named it.
    const char* file () const noexcept { return m_file; }

    /// Name of the function in which the Error was constructed.
    const char* function () const noexcept { return m_function; }

    /// Line at which the Error was constructed.
    int line () const noexcept { return m_line; }

  private:

    const char* m_file;
    const char* m_function;
    int m_line;
  };

  // ==========================================================================
  // Coordinate layout
  // ==========================================================================

  /**
   *  @defgroup coordinate_layout Coordinate layout
   *
   *  @brief How coordinate arrays and the arrays of a result are laid out.
   *
   *  Every coordinate array is flat and row-major, one object after the
   *  other. Cartesian arrays hold 3 * nObjects entries, ordered x, y, z, in
   *  coordinate units. Sky arrays hold 3 * nObjects entries, ordered right
   *  ascension, declination, redshift; angles are in radians.
   *
   *  Every length otswap takes or returns, the mean particle separation,
   *  box limits, displacements and smoothing widths included, is in
   *  coordinate units: the unit of the Cartesian coordinates. In a box it is
   *  the unit of the tracers given; otswap assumes none. In a lightcone it is
   *  the unit of the distance table: Mpc/h when the table comes from the
   *  cosmology constructor of DistanceTable, the unit of the table given when
   *  it comes from the table constructor.
   *
   *  Nothing is copied into a nested container at any point: the layout is
   *  the one the algorithm walks, and the one a numpy array maps onto
   *  without a copy. The arrays of a Result and of a RealSpaceCatalog:
   *
   *      [realization][object][3]  Result::displacement
   *      [realization][object]     Result::valid (uint8, 0 or 1)
   *      [object][3]               Result::tracers, tracersSky, meanDisplacement,
   *                                lagrangian, lagrangianSky;
   *                                RealSpaceCatalog::sky, cartesian
   *      [object]                  Result::ids (int64), validRealizations,
   *                                outsideRedshiftCut, outsideDistanceCut,
   *                                outsideMask, outsideBox (uint8),
   *                                lagrangianDistance;
   *                                RealSpaceCatalog::ids (int64), shift,
   *                                factor, status, validRealizations,
   *                                nNeighbours, nRealizationsAveraged
   *
   *  Flags are uint8, 0 or 1; counts are unsigned; indices of objects are
   *  their rows.
   */

  // ==========================================================================
  // Configuration
  // ==========================================================================

  /**
   *  @brief What a reconstruction or a redshift-space correction reports.
   *
   *  Each line starts with "otswap: " and ends in '\n'; the counts behind
   *  every line are in the object returned, so nothing printed is
   *  unavailable to a program.
   *
   *  Each line is written, and the log flushed, as soon as its step
   *  finishes, by the calling thread, never inside a parallel region; the
   *  lines keep the order given here, so the mps(z) line of a lightcone
   *  that rejects mask crossings follows the line of the crossings. A call
   *  that throws leaves the lines of the steps it finished, and never
   *  writes its call line.
   */
  enum class Verbosity {
    /// Nothing.
    Silent,

    /// The selection report of a lightcone, when a selection applied, and of
    /// a box, when limits were given (SelectionCounts::message), and the
    /// extrapolation of b(z), when it happened; then one line per call with
    /// its time and what was lost:
    ///
    ///     otswap: reconstructBox: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement
    ///     otswap: reconstructLightcone: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement
    ///     otswap: measureMpsProfile: <B> bins over [<zmin>, <zmax>]
    ///     otswap: realSpaceLightcone: corrected <C> of <N> tracers in <T> s; <U> left uncorrected
    ///     otswap: realSpaceBox: corrected <C> of <N> tracers in <T> s; <U> left uncorrected
    Normal,

    /// As Normal, with the lines that explain it: the mean particle
    /// separation of a box and its source, the mps(z) profile of a
    /// lightcone, and the breakdown of the corrected (moved by their
    /// neighbours) and of the uncorrected (left out of the reconstruction,
    /// or with no valid neighbour).
    Detailed
  };

  /**
   *  @brief Free parameters of the reconstruction.
   *
   *  The multipliers that govern the seeding radius, the neighbourhood
   *  size and the initial pairing are absent: they are internal constants
   *  of the algorithm, not parameters.
   */
  struct Config {

    /// Independent reconstructions to run. Each consumes nObjects randoms.
    unsigned nRealizations = 30;

    /// The sweeps stop once the fraction of successful swaps in a sweep,
    /// swaps per tracer visited, no longer exceeds this threshold. Every
    /// tracer is examined in every sweep.
    double convergence = 1.e-3;

    /// Seed of the generator. 0 draws one from std::random_device.
    unsigned seed = 0;

    /// Grid cell size, in units of the representative mean particle
    /// separation. Affects speed only, never the result.
    double cellSize = 3.;

    /// With a Mask given to reconstructLightcone, apply rejectMaskCrossings
    /// to the result, with maxUnobservedPixelsCrossed, before returning it.
    /// Read by the mask overloads of reconstructLightcone only.
    bool rejectCrossings = true;

    /// The threshold of that call, the longest run of unobserved pixels a
    /// displacement's arc may cross; see rejectMaskCrossings. Read by the
    /// mask overloads of reconstructLightcone only.
    unsigned maxUnobservedPixelsCrossed = 0;

    /// What the call reports; see Verbosity.
    Verbosity verbosity = Verbosity::Normal;

    /// Where the report is written; null writes nothing, as Silent.
    std::ostream* log = &std::clog;
  };

  // ==========================================================================
  // Result
  // ==========================================================================

  /**
   *  @brief A closed redshift range for reconstructLightcone. The default,
   *  (-inf, +inf), cuts nothing.
   */
  struct RedshiftCut {
    double min = -std::numeric_limits<double>::infinity();
    double max =  std::numeric_limits<double>::infinity();
  };

  /**
   *  @brief A closed range of comoving distance from the origin, in
   *  coordinate units, for the Cartesian overloads of reconstructLightcone
   *  and measureMpsProfile: the counterpart of RedshiftCut. The default,
   *  (-inf, +inf), cuts nothing.
   */
  struct DistanceCut {
    double min = -std::numeric_limits<double>::infinity();
    double max =  std::numeric_limits<double>::infinity();
  };

  /**
   *  @brief The closed limits of a box for reconstructBox: the points with
   *  lower[k] <= x[k] <= upper[k] on each axis k, in coordinate units.
   *  Each bound must
   *  be finite, with lower[k] < upper[k]; the default, NaN, is refused.
   *
   *  For a periodic box, the period of axis k is upper[k] - lower[k], and
   *  each period must be an integer multiple of the shortest one, within a
   *  relative 1e-12: a cube, or a cuboid such as 500 x 500 x 1000.
   */
  struct BoxLimits {
    /// x, y, z of the lower corner.
    std::array<double, 3> lower {std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN()};
    /// x, y, z of the upper corner.
    std::array<double, 3> upper {std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN()};
  };

  /**
   *  @brief The objects a reconstruction left out, and why, and the
   *  displacements its mask filter rejected.
   *
   *  In a lightcone result, the tracers removed are
   *  tracersOutsideRedshiftCut + tracersOutsideMask - tracersOutsideBoth
   *  from a sky overload of reconstructLightcone, and
   *  tracersOutsideDistanceCut + tracersOutsideMask - tracersOutsideBoth
   *  from a Cartesian one, and the same for the randoms; the counts of the
   *  other cut are 0. In a box result, the tracers removed are
   *  tracersOutsideBox, and the randoms removed randomsOutsideBox; these
   *  and the counts of what was given are set only when limits were given
   *  to reconstructBox, with randoms 0 when the randoms were drawn. In a
   *  periodic box nothing is removed: the objects outside the limits are
   *  wrapped into them, counted in tracersWrapped and randomsWrapped, and
   *  the counts outside the box are 0. Every other count is 0, and nothing
   *  else is applied, in a box result; the box counts are 0 in a lightcone
   *  result.
   */
  struct SelectionCounts {

    std::optional<RedshiftCut> redshiftCut;  ///< the cut, when a bound was finite
    std::optional<DistanceCut> distanceCut;  ///< the distance cut, when a bound was finite
    bool maskApplied = false;                ///< a Mask was given
    std::optional<BoxLimits> boxLimits;      ///< the limits, when given to reconstructBox
    bool periodic = false;                   ///< the box of reconstructBox was periodic

    std::size_t tracers = 0;                    ///< tracers given
    std::size_t tracersOutsideRedshiftCut = 0;  ///< of them, outside the cut
    std::size_t tracersOutsideDistanceCut = 0;  ///< of them, outside the distance cut
    std::size_t tracersOutsideMask = 0;         ///< of them, on an unobserved pixel
    std::size_t tracersOutsideBoth = 0;         ///< outside the cut, redshift or distance, and the mask
    std::size_t tracersOutsideBox = 0;          ///< of them, outside the box limits
    std::size_t tracersWrapped = 0;             ///< of them, wrapped into the periodic box

    std::size_t randoms = 0;                    ///< randoms given
    std::size_t randomsOutsideRedshiftCut = 0;  ///< of them, outside the cut
    std::size_t randomsOutsideDistanceCut = 0;  ///< of them, outside the distance cut
    std::size_t randomsOutsideMask = 0;         ///< of them, on an unobserved pixel
    std::size_t randomsOutsideBoth = 0;         ///< outside the cut, redshift or distance, and the mask
    std::size_t randomsOutsideBox = 0;          ///< of them, outside the box limits
    std::size_t randomsWrapped = 0;             ///< of them, wrapped into the periodic box

    /// The threshold of the latest rejectMaskCrossings applied to the
    /// result, by reconstructLightcone or afterwards; empty when none was.
    std::optional<unsigned> maxUnobservedPixelsCrossed;
    std::size_t displacements = 0;              ///< valid displacements before the first filter
    std::size_t displacementsCrossingMask = 0;  ///< of them, rejected by every filter so far

    /**
     *  @brief The selection report of a reconstruction.
     *
     *  One line for the tracers and one for the randoms, each listing the
     *  selections applied, then one for the rejected crossings when a
     *  filter was applied; each line ends in '\n'. A line is written even
     *  when nothing was removed. In the form:
     *
     *      otswap: kept <n> of <N> tracers: <a> outside the redshift cut [<zmin>, <zmax>], <b> outside the mask (<c> outside both)
     *      otswap: kept <m> of <M> randoms: <a> outside the redshift cut [<zmin>, <zmax>], <b> outside the mask (<c> outside both)
     *      otswap: rejected <r> of <D> displacements crossing more than <L> unobserved pixels in a row
     *
     *  With a distance cut, from a Cartesian overload, the cut reads
     *  "<a> outside the distance cut [<dmin>, <dmax>] in coordinate units".
     *
     *  For a box with limits, the tracer line, and the random line when
     *  randoms were given:
     *
     *      otswap: kept <n> of <N> tracers: <a> outside the box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]
     *      otswap: kept <m> of <M> randoms: <b> outside the box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]
     *
     *  For a periodic box, in their place:
     *
     *      otswap: wrapped <a> of <N> tracers into the periodic box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]
     *      otswap: wrapped <b> of <M> randoms into the periodic box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]
     *
     *  @return the report, or an empty string when no selection was
     *  applied.
     */
    std::string message () const;
  };

  class DistanceTable;

  /// The geometry a Result was reconstructed in.
  enum class Geometry {
    Box,       ///< reconstructBox
    Lightcone  ///< reconstructLightcone
  };

  /**
   *  @brief The mean particle separation of a lightcone as a function of
   *  redshift or of comoving distance, measured by reconstructLightcone, or
   *  alone by measureMpsProfile.
   *
   *  The overloads that take sky coordinates bin the tracers kept in nBins
   *  uniform redshift bins over [redshiftMin, redshiftMax]: a profile in
   *  redshift, with redshift filled and distance empty. The overloads that
   *  take Cartesian coordinates bin them in nBins uniform bins of comoving
   *  distance from the origin over [distanceMin, distanceMax]: a profile in
   *  distance, with distance filled and redshift empty. In each bin
   *  mps = (N / V)^(-1/3), V the shell volume implied by the sky area. The
   *  bin values are the nodes of a piecewise linear profile, extrapolated
   *  linearly beyond the terminal ones; at() evaluates it. Empty in a box
   *  result.
   */
  struct MpsProfile {
    std::vector<double>      redshift;   ///< bin centres of a profile in redshift
    std::vector<double>      mps;        ///< in coordinate units, at each centre
    std::vector<std::size_t> count;      ///< tracers in each bin

    /// The binned range of a profile in redshift: the smallest and largest
    /// redshift of the tracers kept. NaN in a profile in distance.
    double redshiftMin = std::numeric_limits<double>::quiet_NaN();
    double redshiftMax = std::numeric_limits<double>::quiet_NaN();

    /// The count-weighted median of the nodes, which sizes the grids of the
    /// reconstruction.
    double representative = std::numeric_limits<double>::quiet_NaN();

    /// Bin centres of a profile in distance, in coordinate units.
    std::vector<double>      distance;

    /// The binned range of a profile in distance: the smallest and largest
    /// distance of the tracers kept, in coordinate units. NaN in a profile
    /// in redshift.
    double distanceMin = std::numeric_limits<double>::quiet_NaN();
    double distanceMax = std::numeric_limits<double>::quiet_NaN();

    /**
     *  @brief The mps at @p value, in coordinate units: a redshift for a profile
     *  in redshift, a comoving distance for a profile in distance. Linear
     *  between the nodes, extrapolated linearly beyond the terminal ones
     *  along the terminal segment, constant with a single node. This is the
     *  value the reconstruction used at each tracer.
     *
     *  @exception Error if the profile is empty or holds both redshift and
     *  distance nodes, or if the value is not positive, which only an
     *  extrapolation can give.
     */
    double at (double value) const;
  };

  /**
   *  @brief Displacement field produced by a reconstruction.
   *
   *  Displacements point from the Eulerian (observed) position to the
   *  Lagrangian one: the Lagrangian position of object i is its
   *  coordinate plus its displacement.
   */
  struct Result {

    std::size_t nObjects = 0;
    unsigned    nRealizations = 0;

    /// The number of realizations averaged into each displacement, when
    /// the result holds a mean as its only realization, as
    /// io::readDisplacementMean returns it; 1, the default, when it does
    /// not, as in every result of a reconstruction. A value above 1
    /// requires nRealizations to be 1, and validRealizations then holds
    /// counts up to nRealizationsMean. The displacement writers record a
    /// value above 1 as the keyword NRECMEAN, with NREC = 1. No computation
    /// reads it.
    unsigned    nRealizationsMean = 1;

    /// The geometry of the reconstruction, which the writers and the
    /// redshift-space correction follow.
    Geometry geometry = Geometry::Box;

    /// The ID of each object, as the caller sets it: any 64-bit integer,
    /// negative ones included, each object its own. Empty means no IDs, and
    /// the writers then write object i's row in the input, i, from 0. Flat,
    /// [object]; size nObjects or 0. The reconstruction functions leave it
    /// empty, and no computation reads it: realSpaceBox and
    /// realSpaceLightcone copy it into RealSpaceCatalog::ids, and the writers
    /// write it as the id column, among their default columns when it is not
    /// empty.
    std::vector<std::int64_t> ids;

    /// The configuration the reconstruction ran with. seed is the one used:
    /// the drawn one when 0 was given, so that the result can be
    /// reproduced.
    Config config;

    /// Cartesian position of each tracer, the start of its displacements:
    /// the array given to reconstructBox or to a Cartesian overload of
    /// reconstructLightcone, or toCartesian of the sky array given to a sky
    /// overload. In a periodic
    /// box, the tracers wrapped into the box: a coordinate outside the
    /// limits is shifted by a whole number of periods, one inside them is
    /// kept as given. A tracer left
    /// out by the redshift cut or the mask of a sky overload has its
    /// conversion when its redshift lies in the distance table, NaN
    /// otherwise. Flat,
    /// [object][xyz]; size 3 * nObjects.
    std::vector<double> tracers;

    /// Sky coordinates of each tracer: right ascension, declination,
    /// redshift; radians. As given to a sky overload of
    /// reconstructLightcone; for a Cartesian overload, the right ascension
    /// and declination of pointToSky, and the redshift of the distance
    /// table given, NaN without one or where the distance lies outside it.
    /// Flat, [object][3]; size 3 * nObjects in a lightcone result, empty in
    /// a box result.
    std::vector<double> tracersSky;

    /// Displacements: the random matched to each object in each
    /// realization, minus the tracer. Flat, [realization][object][xyz].
    /// Size 3 * nRealizations * nObjects. In a periodic box, the minimum
    /// image: on each axis, d - L floor(d / L + 1/2), d the difference and
    /// L the period, in [-L/2, L/2) within rounding.
    std::vector<double> displacement;

    /// Validity of each displacement: 1 valid, 0 rejected by a filter.
    /// Flat, [realization][object]. Size nRealizations * nObjects. Every
    /// entry is 1 until a filter is applied; a filter only ever clears
    /// entries.
    std::vector<std::uint8_t> valid;

    /// Mean displacement per object, over its valid realizations only,
    /// summed in realization order, then divided by their count. NaN in all
    /// three components for an object with no valid realization. Flat,
    /// [object][xyz]. Size 3 * nObjects.
    std::vector<double> meanDisplacement;

    /// Valid realizations per object: the count of 1 entries of valid for
    /// that object, and the number of terms in its meanDisplacement. In a
    /// result holding a mean (nRealizationsMean above 1), the count of the
    /// realizations averaged into it that were valid. Always at most
    /// nRealizations x nRealizationsMean. Size nObjects.
    std::vector<unsigned> validRealizations;

    /// NSIDE of the mask this result was last filtered against, or 0 when
    /// it has not been filtered. Set by rejectMaskCrossings, which refuses
    /// to filter the same result against a mask of a different NSIDE.
    int filteredNside = 0;

    /// Whether each tracer lay outside the redshift cut of
    /// reconstructLightcone and so took no part in the reconstruction:
    /// 1 for such a tracer, 0 otherwise. Flat, [object]. Size nObjects, all
    /// 0 when no cut was given and in a box result.
    ///
    /// A cut tracer has NaN in every component of its displacement rows,
    /// valid 0 in every realization, validRealizations 0 and
    /// meanDisplacement NaN. No filter changes the flag. A Result
    /// built by hand may leave the field empty, which every function reads
    /// as all 0.
    std::vector<std::uint8_t> outsideRedshiftCut;

    /// Whether each tracer lay outside the distance cut of a Cartesian
    /// overload of reconstructLightcone, and so took no part in the
    /// reconstruction: 1 for such a tracer, 0 otherwise. Flat, [object].
    /// Size nObjects, all 0 when no cut was given, in a result of a sky
    /// overload and in a box result. A flagged tracer's rows are as for
    /// outsideRedshiftCut, and no filter changes the flag. An empty field is
    /// read as all 0.
    std::vector<std::uint8_t> outsideDistanceCut;

    /// Whether each tracer fell on an unobserved pixel of the Mask given to
    /// reconstructLightcone, and so took no part in the reconstruction: 1
    /// for such a tracer, 0 otherwise. Flat, [object]. Size nObjects, all 0
    /// when no mask was given and in a box result. Independent of
    /// outsideRedshiftCut: a tracer may carry both flags. A flagged
    /// tracer's rows are as for outsideRedshiftCut, and no filter changes
    /// the flag. An empty field is read as all 0.
    std::vector<std::uint8_t> outsideMask;

    /// Whether each tracer lay outside the limits given to reconstructBox,
    /// and so took no part in the reconstruction: 1 for such a tracer, 0
    /// otherwise. Flat, [object]. Size nObjects, all 0 when no limits were
    /// given, in a periodic box and in a lightcone result. A flagged
    /// tracer's rows are as for
    /// outsideRedshiftCut, and no filter changes the flag. An empty field is
    /// read as all 0.
    std::vector<std::uint8_t> outsideBox;

    /// What the reconstruction left out, and why, and what the mask filters
    /// rejected; updated by every rejectMaskCrossings.
    SelectionCounts selection;

    /// Mean Lagrangian position of each tracer, tracers + meanDisplacement,
    /// in a periodic box wrapped into the box as the tracers are.
    /// NaN rows where meanDisplacement is NaN. Flat, [object][xyz]; size
    /// 3 * nObjects. Recomputed with meanDisplacement by rejectMaskCrossings
    /// and recomputeMeans.
    std::vector<double> lagrangian;

    /// Sky coordinates of each tracer's mean Lagrangian position, as toSky
    /// computes them from lagrangian with the reconstruction's distance
    /// table: right ascension in [0, 2 pi), declination, redshift; radians.
    /// Flat, [object][3]; size 3 * nObjects in a lightcone result with a
    /// distance table or a profile in distance, empty otherwise (a box
    /// result included). Recomputed with lagrangian.
    ///
    /// NaN rows where lagrangian is NaN: the tracers flagged in
    /// outsideRedshiftCut, outsideDistanceCut or outsideMask and those with no valid
    /// realization. A position whose distance falls outside the table,
    /// which a table that does not start at distance 0 allows, has a NaN
    /// redshift and keeps its right ascension and declination. Without a
    /// table, as a Cartesian overload of reconstructLightcone given none
    /// returns it, every redshift is NaN and the angles are filled.
    std::vector<double> lagrangianSky;

    /// Comoving distance from the origin of each tracer's mean Lagrangian
    /// position, as pointToSky computes it from lagrangian, in coordinate
    /// units. Flat, [object]; size nObjects when lagrangianSky is filled,
    /// empty otherwise; NaN where lagrangian is NaN. Recomputed with
    /// lagrangian.
    std::vector<double> lagrangianDistance;

    /// Box: the mean particle separation used, given or computed, in
    /// coordinate units.
    /// NaN in a lightcone result, whose separation is mpsProfile.
    double mps = std::numeric_limits<double>::quiet_NaN();

    /// Box: the box of the reconstruction, the limits given to
    /// reconstructBox or, without them, the bounding box of the tracers
    /// and the given randoms (of the tracers alone when the randoms were
    /// drawn), also when mps was given. Every random matched to a tracer
    /// lies inside it. Empty in a lightcone result.
    std::optional<BoxLimits> box;

    /// Box: whether the box was periodic, as given to reconstructBox. A
    /// periodic result has geometry Box and its limits in box; its
    /// tracers, displacements and lagrangian follow the periodic box.
    /// Always false in a lightcone result.
    bool periodic = false;

    /// Lightcone: the mean particle separation profile measured on the
    /// tracers kept: in redshift for a sky overload of
    /// reconstructLightcone, in distance for a Cartesian one. Empty in a
    /// box result.
    MpsProfile mpsProfile;

    /// Lightcone: the distance table of the reconstruction, shared by the
    /// copies of this result, with which rejectMaskCrossings and
    /// recomputeMeans recompute lagrangianSky. Null in a box result, and in
    /// the result of a Cartesian overload of reconstructLightcone given no
    /// table.
    std::shared_ptr<const DistanceTable> distances;

    /// Wall time of the call that made the result, in seconds. Two runs
    /// with the same seed give the same result but for this field.
    double elapsedSeconds = 0.;
  };

  // ==========================================================================
  // Cosmology
  // ==========================================================================

  /**
   *  @brief Sampled redshift to comoving distance relation, interpolated
   *  in both directions, and the linear growth rate on the same grid.
   *
   *  The object is immutable once built and holds no external state, so
   *  one instance can serve any number of reconstructions concurrently.
   */
  class DistanceTable {

  public:

    /**
     *  @brief Build from the parameters of a flat cosmology.
     *
     *  Samples nSamples points over [zMin, zMax] and integrates 1/E(z).
     *  The defaults of w0 and wa give flat LCDM. Non-flat cosmologies are
     *  served by the table constructor below.
     *
     *  With the default range and sampling, for flat LCDM and w0-wa
     *  cosmologies near it, the relative error of distanceAt, redshiftAt
     *  and growthRateAt is below 1e-6 for z >= 0.01, and that of
     *  growthRateAt over the whole range; the absolute error of distanceAt,
     *  and that of redshiftAt expressed as a distance along the line of
     *  sight, is below 2e-5 Mpc/h over [0, 10].
     *  Below z = 0.01 the relative error of the distances grows, because
     *  the distance itself vanishes at z = 0 while the absolute error does
     *  not. The table holds three arrays of nSamples doubles.
     *
     *  The growth rate f = dlnD/dlna is sampled on the same grid, by
     *  integrating the linear growth equation from deep in matter
     *  domination and taking f from the solution. It is not approximated
     *  by Omega_m(z) raised to a power: that form is fitted to LCDM and
     *  these parameters describe a w0-wa cosmology.
     *
     *  @param OmegaM matter density parameter; must lie in (0, 1]
     *  @param h Hubble constant in units of 100 km/s/Mpc; must be positive
     *  @param w0 dark energy equation of state, constant term; default -1;
     *  must be finite
     *  @param wa dark energy equation of state, evolution term; default 0;
     *  must be finite
     *  @param zMin lower end of the sampled range; must be non-negative
     *  @param zMax upper end; must exceed zMin
     *  @param nSamples sampling points, uniform in redshift; must be at
     *  least two
     *
     *  @exception Error if any parameter is out of range.
     */
    DistanceTable (double OmegaM, double h, double w0 = -1., double wa = 0.,
                   double zMin = 0., double zMax = 10., unsigned nSamples = 50000);

    /**
     *  @brief Build from a table the caller has computed.
     *
     *  @param redshift strictly increasing, at least two entries
     *  @param distance strictly increasing, same length, comoving, in the
     *  unit that Cartesian coordinates and every length will then have
     *
     *  @param growthRate f = dlnD/dlna on the same grid, or empty. When
     *  empty the table carries no growth rate and growthRateAt throws.
     *  Need not be monotonic; every entry must be finite and positive.
     *
     *  @exception Error if either of the first two arrays is not strictly
     *  increasing, if the lengths differ, if fewer than two points are
     *  given, if an entry is not finite, or if growthRate is neither empty
     *  nor of the same length, or holds a value that is not finite and
     *  positive.
     */
    DistanceTable (std::vector<double> redshift, std::vector<double> distance,
                   std::vector<double> growthRate = {});

    /// Comoving distance at z. Throws outside the sampled range. A value
    /// outside it by at most 8 rounding errors of the range's largest
    /// magnitude is taken to be at the nearer end, so that values computed
    /// from the table's own ends, such as the radius of a position
    /// converted at the largest redshift, are accepted. The same holds for
    /// redshiftAt and growthRateAt.
    double distanceAt (double z) const;

    /// Redshift at a comoving distance. Throws outside the sampled range,
    /// with the end tolerance of distanceAt.
    double redshiftAt (double distance) const;

    /// Growth rate f = dlnD/dlna at z. Throws outside the sampled range,
    /// with the end tolerance of distanceAt, and when the table carries
    /// no growth rate.
    double growthRateAt (double z) const;

    /// True when the table carries a growth rate.
    bool hasGrowthRate () const;

    double minRedshift () const;
    double maxRedshift () const;

    /// The sampled redshifts, strictly increasing, as given to the table
    /// constructor or as the cosmology constructor sampled them. A table
    /// built from redshift(), distance() and growthRate() with the table
    /// constructor gives the same value, bit for bit, at every lookup.
    const std::vector<double>& redshift () const;

    /// The comoving distance at each sampled redshift: Mpc/h from the
    /// cosmology constructor, the unit of the table given to the table
    /// constructor.
    const std::vector<double>& distance () const;

    /// The growth rate at each sampled redshift; empty when the table
    /// carries none.
    const std::vector<double>& growthRate () const;

  private:

    std::vector<double> m_redshift;
    std::vector<double> m_distance;
    std::vector<double> m_growthRate;
  };

  // ==========================================================================
  // Angular mask
  // ==========================================================================

  /// Pixel ordering of a full-sky Healpix map given in memory.
  enum class PixelOrdering {
    Ring,    ///< the RING scheme
    Nested   ///< the NESTED scheme
  };

  /**
   *  @brief Healpix veto mask, read from a FITS file or built from a
   *  full-sky map in memory.
   *
   *  A pixel is allowed (observed) when its value is greater than 0,
   *  whatever its magnitude: fractional values, the smallest positive ones
   *  and +inf are allowed. A pixel holding 0, a negative value, -inf, NaN
   *  or Healpix's UNSEEN (-1.6375e30) is unobserved. Every value is
   *  accepted, and each is tested once, in double, when the mask is built;
   *  none is kept or used as a weight, so skyAreaDeg2 is a count of
   *  allowed pixels, not a sum of their values. The rule is the same for
   *  both sources.
   *
   *  NSIDE may be any integer from 1 to 2^29, Healpix's limit, and must be
   *  a power of 2 in the NESTED scheme. The mask holds one byte per pixel,
   *  12 NSIDE^2 bytes: 12.6 MB at NSIDE 1024, 805 MB at NSIDE 8192, 3.2 GB
   *  at NSIDE 16384. A map too large for memory raises std::bad_alloc.
   *
   *  The object is immutable once built.
   */
  class Mask {

  public:

    /**
     *  @param fitsFile Healpix map: a binary table in HDU 2 with 1024
     *  values to a row, in any numeric column type; NSIDE and ORDERING are
     *  taken from the header, in either RING or NESTED scheme. The values
     *  are read a block of rows at a time, so reading needs memory for the
     *  mask's bytes only. Each value is tested as cfitsio converts it to
     *  double, after any TSCAL and TZERO; a null value of an integer
     *  column (TNULL) is unobserved.
     *
     *  @exception Error if the file cannot be opened or read, if it is not
     *  a full-sky Healpix map in that layout, if NSIDE exceeds 2^29, or if
     *  the ordering is NESTED and NSIDE is not a power of 2.
     */
    explicit Mask (const std::string& fitsFile);

    /**
     *  @brief Build from the pixel values of a full-sky Healpix map.
     *
     *  NSIDE follows from the number of values, which must be 12 NSIDE^2.
     *  The values are read once, during construction, and not kept.
     *
     *  @param values one value per pixel, in the order of @p ordering
     *  @param count the number of values
     *  @param ordering the pixel ordering of @p values; no default
     *
     *  @exception Error if @p count is 0 or not 12 NSIDE^2 for an integer
     *  NSIDE, if NSIDE exceeds 2^29, if the ordering is NESTED and NSIDE is
     *  not a power of 2, or if @p values is null; the message gives the
     *  number of values and the NSIDE it implies.
     */
    Mask (const double* values, std::size_t count, PixelOrdering ordering);

    /// As the constructor from a pointer and a count, over all of
    /// @p values.
    Mask (const std::vector<double>& values, PixelOrdering ordering);

    /// True when the direction falls in an allowed pixel. Radians.
    /// @exception Error if the right ascension is not finite, or the
    /// declination is not in [-pi/2, pi/2].
    bool allows (double rightAscension, double declination) const;

    /// Sky area covered by the allowed pixels, in square degrees.
    double skyAreaDeg2 () const;

    /// NSIDE of the map, from 1 to 2^29.
    int nside () const;

  private:

    class Impl;
    std::shared_ptr<const Impl> m_impl;

    friend void rejectMaskCrossings (Result&, const Mask&, unsigned);
  };

  // ==========================================================================
  // Reconstruction
  // ==========================================================================

  /**
   *  @brief Reconstruct in box geometry, with a constant mean particle
   *  separation.
   *
   *  @param tracers Cartesian coordinates, 3 * nObjects entries.
   *
   *  @param randoms Cartesian coordinates of the randoms, or empty. When
   *  empty, config.nRealizations * N randoms are drawn uniformly in the
   *  box limits, or without them in the bounding box of the tracers, N the
   *  tracers reconstructed. When given, the array must hold at least
   *  config.nRealizations * N of them, both counted inside the limits when
   *  they are given, since each realization consumes a disjoint subset;
   *  without limits they may lie outside the bounding box of the tracers.
   *
   *  @param mps mean particle separation; must be positive.
   *
   *  @param config the parameters of the run; rejectCrossings and
   *  maxUnobservedPixelsCrossed are not read.
   *
   *  @param limits the box of the reconstruction, or none. Given, the
   *  tracers outside it take no part in the reconstruction and are flagged
   *  in Result::outsideBox, the randoms given outside it are dropped, and
   *  the neighbour grids span it. A point on a face is inside.
   *
   *  @param periodic whether the box is periodic; requires @p limits, each
   *  of whose sides must then be an integer multiple of the shortest (see
   *  BoxLimits). Every distance is then the minimum image, and nothing is
   *  left out: each tracer and each given random outside the limits is
   *  wrapped into them, by a whole number of periods per axis, and counted
   *  in Result::selection. Result::tracers holds the wrapped tracers, and
   *  the randoms are counted and drawn as with limits.
   *
   *  @exception Error if any array has a size that is not a multiple of
   *  three, if the tracer array is empty, if the randoms are too few, if
   *  any coordinate is not finite, if a bound of the limits is not finite
   *  or a lower bound is not below the upper one, if the limits keep too
   *  few tracers or randoms, or if mps is not positive; with @p periodic,
   *  if no limits are given or a side is not an integer multiple of the
   *  shortest. Also, from
   *  the neighbour search (see Error), if the grid cannot be built: a cell
   *  of config.cellSize * mps so small against the extent of the tracers
   *  and randoms that the grid would need more cells than it can count.
   *
   *  @return the displacement field, one row per tracer given; every entry
   *  of valid is 1 and validRealizations is config.nRealizations, but for
   *  the tracers outside the limits, whose rows are as for
   *  Result::outsideBox, and meanDisplacement is the mean over all
   *  realizations. Result::mps is mps, and Result::box the limits or,
   *  without them, the bounding box of the tracers and the given randoms,
   *  or of the tracers alone when the randoms are drawn.
   *
   *  @note Unless @p periodic, nothing flows through the faces of the
   *  box, so modes on the scale of the box itself cannot be reconstructed.
   *
   *  @note In a periodic box the cells of the neighbour grids divide the
   *  shortest side, so they may be smaller than config.cellSize * mps.
   *
   *  @note Given a seed, the result does not depend on the number of
   *  threads. The assignment of given randoms to realizations, and in
   *  each realization the random draw, the seeding pass and the swap
   *  loop, draw from separate streams, each seeded from the run seed,
   *  the realization index and the stream alone.
   */
  Result reconstructBox (const std::vector<double>& tracers,
                         const std::vector<double>& randoms,
                         double mps,
                         const Config& config,
                         const std::optional<BoxLimits>& limits = std::nullopt,
                         bool periodic = false);

  /**
   *  @brief Reconstruct in box geometry, with the mean particle separation
   *  of the tracers: mps = (V / N)^(1/3), N the tracers inside the box and
   *  V its volume, the box being Result::box: the limits or, without them,
   *  the bounding box of the tracers and the given randoms, or of the
   *  tracers alone when the randoms are drawn (the box drawn randoms
   *  fill). The value used is Result::mps. In a periodic box, V is the
   *  volume of the limits and N every tracer.
   *
   *  @exception Error as above, from the neighbour search included, and,
   *  without limits, if the tracers, with the given randoms, span no
   *  volume (all of them in one plane, on one line or at one point).
   */
  Result reconstructBox (const std::vector<double>& tracers,
                         const std::vector<double>& randoms,
                         const Config& config,
                         const std::optional<BoxLimits>& limits = std::nullopt,
                         bool periodic = false);

  /**
   *  @brief Convert sky coordinates to Cartesian comoving coordinates.
   *
   *  Object i at right ascension ra, declination dec and redshift z maps
   *  to x = d cos(dec) cos(ra), y = d cos(dec) sin(ra), z = d sin(dec),
   *  with d = distances.distanceAt(z). reconstructLightcone converts its
   *  sky arrays with this function, so Cartesian arrays made with it are
   *  exactly the ones that overload uses.
   *
   *  @param sky sky coordinates, 3 * nObjects entries, ordered right
   *  ascension, declination, redshift; angles in radians. May be empty.
   *  @param distances the redshift to distance relation.
   *
   *  @return Cartesian coordinates in the unit of the distance table,
   *  3 * nObjects entries ordered x, y, z.
   *
   *  @exception Error if the array's size is not a multiple of three, if
   *  an entry is not finite, if a declination lies outside [-pi/2, pi/2],
   *  or if a redshift falls outside the distance table; the message names
   *  the object.
   */
  std::vector<double> toCartesian (const std::vector<double>& sky,
                                   const DistanceTable& distances);

  /**
   *  @brief Convert Cartesian comoving coordinates to sky coordinates: the
   *  inverse of toCartesian, with the same observer at the origin.
   *
   *  Object i at (x, y, z) maps to the right ascension atan2(y, x) folded
   *  into [0, 2 pi), the declination asin(z / d), and the redshift at which
   *  the table's comoving distance is d = |(x, y, z)|. A position at the
   *  origin gets right ascension and declination 0. toSky(toCartesian(sky))
   *  returns sky within a few rounding errors, with the right ascension
   *  folded into [0, 2 pi).
   *
   *  @param cartesian 3 * nObjects entries, ordered x, y, z, in the
   *  unit of the distance table. May be empty. A row holding a NaN gives a NaN row.
   *  @param distances the redshift to distance relation; its range bounds
   *  the distances accepted.
   *
   *  @return sky coordinates, 3 * nObjects entries, ordered right
   *  ascension, declination, redshift; angles in radians.
   *
   *  @exception Error if the array's size is not a multiple of three, if an
   *  entry is infinite, or if a distance falls outside the table (beyond
   *  its end tolerance); the message names the object and, for the
   *  distance, suggests a table covering a wider redshift range.
   */
  std::vector<double> toSky (const std::vector<double>& cartesian,
                             const DistanceTable& distances);

  /**
   *  @brief Comoving distance, right ascension and declination of one
   *  Cartesian point, the observer at the origin; no distance table.
   *
   *  The geometric part of toSky, bit for bit: distance = |(x, y, z)|, the
   *  right ascension atan2(y, x) folded into [0, 2 pi), the declination
   *  asin(z / distance), the ratio clamped to [-1, 1]. The origin gets
   *  right ascension and declination 0. A NaN coordinate gives NaN in all
   *  three, as a NaN row of toSky does; an infinite one gives an infinite
   *  distance. Never throws.
   *
   *  @param x, y, z the point, in coordinate units.
   *  @param[out] distance the distance from the origin, in coordinate units.
   *  @param[out] ra the right ascension, radians.
   *  @param[out] dec the declination, radians.
   */
  void pointToSky (double x, double y, double z, double& distance, double& ra, double& dec) noexcept;

  /**
   *  @brief Right ascension and declination of a sky array from degrees to
   *  radians, multiplied by pi/180; redshifts unchanged.
   *
   *  The factor is the double numpy.deg2rad multiplies by, so an array
   *  converted here and one converted by the Python package are the same.
   *
   *  @param skyDegrees 3 * nObjects entries, ordered right ascension,
   *  declination, redshift; angles in degrees. May be empty.
   *
   *  @exception Error if the size is not a multiple of three, or a finite
   *  declination lies outside [-90, 90] degrees; the message names the
   *  object.
   */
  std::vector<double> skyToRadians (std::vector<double> skyDegrees);

  /**
   *  @brief Right ascension and declination of a sky array from radians to
   *  degrees, multiplied by 180/pi; redshifts unchanged. A right ascension
   *  that rounds to 360 is folded to 0, so that one in [0, 2 pi) stays in
   *  [0, 360). NaN stays NaN.
   *
   *  @exception Error if the size is not a multiple of three.
   */
  std::vector<double> skyToDegrees (std::vector<double> skyRadians);

  /**
   *  @brief Reconstruct in lightcone geometry, from sky coordinates.
   *
   *  The mean particle separation is measured from the tracers
   *  themselves: nBins uniform bins in redshift over the observed range,
   *  mps = (N / V)^(-1/3) per bin, with V the shell volume implied by
   *  skyAreaDeg2 and the distance table. A bin holding fewer than
   *  1/(9 eps^2) tracers, eps = 0.02, is refused: below that count the
   *  Poisson error on its mps exceeds eps. The bin values are the nodes
   *  of a piecewise linear profile, extrapolated linearly beyond the
   *  terminal ones.
   *
   *  Randoms are required and are not generated: they carry the survey
   *  geometry, the selection function and the completeness, none of which
   *  this package models.
   *
   *  @param tracersSky sky coordinates, 3 * nObjects entries, radians.
   *
   *  @param randomsSky sky coordinates of the randoms. Their number is
   *  checked after the redshift cut and, in the overloads with a mask, the
   *  mask: the randoms kept must number at least config.nRealizations
   *  times the tracers kept, which is config.nRealizations * nObjects when
   *  nothing is left out. They may lie outside the Cartesian bounding box
   *  of the tracers.
   *
   *  @param skyAreaDeg2 effective area of the survey. When a mask is
   *  available, pass the mask instead, to the overload below, which takes
   *  the area from it; otherwise supply the published value: a bounding box
   *  in right ascension and declination overestimates it for any footprint
   *  that is not rectangular, and the error propagates into the mean
   *  particle separation.
   *
   *  @param nBins redshift bins used to measure mps(z).
   *
   *  @param distances the redshift to distance relation, with which the sky
   *  coordinates are converted and the shell volumes computed; it must
   *  cover the redshifts of the objects kept.
   *
   *  @param config the parameters of the run; rejectCrossings and
   *  maxUnobservedPixelsCrossed are read by the mask overloads only.
   *
   *  @param cut tracers and randoms whose redshift lies outside
   *  [cut.min, cut.max] are left out of the reconstruction, before
   *  anything else: the randoms are dropped, and the tracers keep their
   *  row, flagged in Result::outsideRedshiftCut. mps(z) is measured on the
   *  tracers kept, and only their redshifts need lie in the distance
   *  table. The default cuts nothing.
   *
   *  The counts are recorded in Result::selection and, with a cut that
   *  has a finite bound, reported as config.verbosity says.
   *
   *  @exception Error if a bin holds too few tracers for its density to
   *  be meaningful; the message reports how many bins the catalog
   *  supports. Also if the profile extrapolates to a mps that is not
   *  positive at a tracer's redshift, if any array is malformed, if a
   *  declination lies outside [-pi/2, pi/2], if the randoms are too few,
   *  or if a redshift falls outside the distance table. Also if a bound of
   *  the cut is NaN or min > max, or if too few tracers or randoms remain
   *  inside it; the message gives the counts kept and dropped. Also, from
   *  the neighbour search (see Error), if the grid cannot be built: a cell
   *  of config.cellSize mean particle separations so small against the
   *  extent of the tracers and randoms that the grid would need more cells
   *  than it can count.
   */
  Result reconstructLightcone (const std::vector<double>& tracersSky,
                               const std::vector<double>& randomsSky,
                               double skyAreaDeg2,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config,
                               const RedshiftCut& cut = {});

  /**
   *  @brief Reconstruct in lightcone geometry, from Cartesian comoving
   *  coordinates, the observer at the origin; no cosmology.
   *
   *  As the overload from sky coordinates, with the mean particle
   *  separation measured in distance rather than redshift: nBins uniform
   *  bins of comoving distance from the origin over the range of the
   *  tracers, mps = (N / V)^(-1/3) per bin, V the shell volume implied by
   *  skyAreaDeg2, with the same per-bin threshold. Result::mpsProfile is a
   *  profile in distance, and each tracer's mps is its value at the
   *  tracer's distance.
   *
   *  The result is filled as by the sky overload, with tracers the array
   *  given and tracersSky, lagrangianSky and lagrangianDistance computed
   *  with pointToSky; the redshifts come from @p distances when one is
   *  given and are NaN otherwise, or where a distance lies outside it.
   *
   *  @param tracers Cartesian coordinates, 3 * nObjects entries, in
   *  coordinate units.
   *  @param randoms Cartesian coordinates of the randoms, required, as in
   *  the sky overload, in the same units.
   *  @param skyAreaDeg2 as in the sky overload.
   *  @param nBins distance bins used to measure mps(d).
   *  @param config as in the sky overload.
   *  @param distances optional; the table whose redshifts fill those of
   *  tracersSky and lagrangianSky, and which the result keeps. It plays no
   *  part in the reconstruction: the result is the same with or without
   *  it but for those redshifts and Result::distances.
   *  @param cut the redshift cut of the sky overload, in comoving distance
   *  from the origin: tracers and randoms whose distance lies outside
   *  [cut.min, cut.max] are left out before anything else, the randoms
   *  dropped and the tracers kept in their row, flagged in
   *  Result::outsideDistanceCut; mps(d) is measured on the tracers kept.
   *  The counts are recorded in Result::selection and, with a cut that has
   *  a finite bound, reported as config.verbosity says. The default cuts
   *  nothing.
   *
   *  @exception Error if a bin holds too few tracers, if the profile
   *  extrapolates to a mps that is not positive at a tracer's distance, if
   *  every tracer lies at the same distance, if any array is malformed or
   *  holds a non-finite coordinate, or if the randoms are too few; also
   *  from the neighbour search, as in the sky overload. Also if a bound of
   *  the cut is NaN or min > max, or if too few tracers or randoms remain
   *  inside it; the message gives the counts kept and dropped.
   */
  Result reconstructLightcone (const std::vector<double>& tracers,
                               const std::vector<double>& randoms,
                               double skyAreaDeg2,
                               unsigned nBins,
                               const Config& config,
                               const DistanceTable* distances = nullptr,
                               const DistanceCut& cut = {});

  /**
   *  @brief Reconstruct in lightcone geometry, from sky coordinates, within
   *  the observed pixels of a mask.
   *
   *  As the overload with the sky area, with the area taken from the mask,
   *  mask.skyAreaDeg2(), and the mask applied twice:
   *  - before the reconstruction, tracers and randoms on unobserved pixels
   *    are left out, decided on their right ascension and declination as
   *    Mask::allows decides: the randoms are dropped, and the tracers keep
   *    their row, flagged in Result::outsideMask;
   *  - after it, when config.rejectCrossings is set, rejectMaskCrossings is
   *    applied with config.maxUnobservedPixelsCrossed, so the result
   *    carries filteredNside = mask.nside().
   *
   *  The mask and the redshift cut are both evaluated on every object, so a
   *  tracer may be flagged by both; the objects kept are those that pass
   *  both. The counts are recorded in Result::selection and reported as
   *  config.verbosity says.
   *
   *  To filter with another threshold afterwards, set
   *  config.maxUnobservedPixelsCrossed instead: a second rejectMaskCrossings
   *  with the same mask and a higher threshold changes nothing, since the
   *  filter only clears entries. Or set config.rejectCrossings to false and
   *  call rejectMaskCrossings on the result.
   *
   *  @param tracersSky as in the overload with the sky area.
   *  @param randomsSky as in the overload with the sky area.
   *  @param nBins as in the overload with the sky area.
   *  @param distances as in the overload with the sky area.
   *  @param mask the survey's veto mask; borrowed for the call only.
   *  @param config as in the overload with the sky area; rejectCrossings
   *  and maxUnobservedPixelsCrossed are read here.
   *  @param cut as in the overload with the sky area, evaluated with the
   *  mask on every object.
   *
   *  @exception Error as the overload with the sky area, from the
   *  neighbour search included, the counts in the message naming the mask
   *  as well as the cut.
   */
  Result reconstructLightcone (const std::vector<double>& tracersSky,
                               const std::vector<double>& randomsSky,
                               const Mask& mask,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config,
                               const RedshiftCut& cut = {});

  /**
   *  @brief Reconstruct in lightcone geometry, from Cartesian comoving
   *  coordinates, within the observed pixels of a mask; no cosmology.
   *
   *  The Cartesian overload with the sky area, with the area taken from the
   *  mask and the mask applied as by the sky overload with a mask: the
   *  objects on unobserved pixels, decided on the right ascension and
   *  declination pointToSky gives, are left out before the reconstruction,
   *  and, when config.rejectCrossings is set, the displacements crossing
   *  unobserved pixels are rejected after it.
   *
   *  @param tracers as in the Cartesian overload with the sky area.
   *  @param randoms as in the Cartesian overload with the sky area.
   *  @param mask the survey's veto mask; borrowed for the call only.
   *  @param nBins as in the Cartesian overload with the sky area.
   *  @param config as in the sky overload with a mask.
   *  @param distances as in the Cartesian overload with the sky area.
   *  @param cut as in the Cartesian overload with the sky area, evaluated
   *  with the mask on every object, as in the sky overload with a mask.
   *
   *  @exception Error as the Cartesian overload with the sky area, and if
   *  too few tracers or randoms remain inside the mask; the counts in the
   *  message name the mask as well as the cut.
   */
  Result reconstructLightcone (const std::vector<double>& tracers,
                               const std::vector<double>& randoms,
                               const Mask& mask,
                               unsigned nBins,
                               const Config& config,
                               const DistanceTable* distances = nullptr,
                               const DistanceCut& cut = {});

  /**
   *  @brief The mps(z) profile reconstructLightcone measures, without
   *  reconstructing: from sky coordinates, with the sky area given.
   *
   *  The steps of reconstructLightcone that concern the tracers, in the
   *  same order and with the same code: the checks of the sky array and of
   *  the cut, the redshift cut, the minimum number of tracers kept, and
   *  the profile with its per-bin threshold. The profile equals, bit for
   *  bit, Result::mpsProfile of a reconstructLightcone call with the same
   *  tracers, skyAreaDeg2, nBins, distances and cut, whatever its randoms
   *  and config. Use it to choose nBins, or to inspect the profile, before
   *  a reconstruction.
   *
   *  @param tracersSky as in reconstructLightcone.
   *  @param skyAreaDeg2 as in reconstructLightcone.
   *  @param nBins as in reconstructLightcone.
   *  @param distances as in reconstructLightcone; it must cover the
   *  redshifts of the tracers kept.
   *  @param cut as in reconstructLightcone, applied to the tracers.
   *
   *  @param verbosity what the call reports. At Normal, when the cut has a
   *  finite bound, the tracer line of the selection report
   *  (SelectionCounts::message), then
   *
   *      otswap: measureMpsProfile: <B> bins over [<zmin>, <zmax>]
   *
   *  At Detailed, the mps(z) line of reconstructLightcone before the last.
   *
   *  @param log where the report is written; null writes nothing.
   *
   *  @exception Error as reconstructLightcone, with the same messages, if
   *  the tracer array is malformed, if a declination lies outside
   *  [-pi/2, pi/2], if a bound of the cut is NaN or min > max, if a cut
   *  with a finite bound keeps too few tracers, if a bin holds too few
   *  tracers, or if a redshift kept falls outside the distance table.
   */
  MpsProfile measureMpsProfile (const std::vector<double>& tracersSky,
                                double skyAreaDeg2,
                                unsigned nBins,
                                const DistanceTable& distances,
                                const RedshiftCut& cut = {},
                                Verbosity verbosity = Verbosity::Normal,
                                std::ostream* log = &std::clog);

  /**
   *  @brief The mps(z) profile reconstructLightcone measures, without
   *  reconstructing: from sky coordinates, within the observed pixels of a
   *  mask.
   *
   *  As the overload with the sky area, with the area taken from the mask,
   *  mask.skyAreaDeg2(), and the tracers on unobserved pixels left out as
   *  reconstructLightcone leaves them out. The profile equals, bit for bit,
   *  Result::mpsProfile of the mask overloads of reconstructLightcone with
   *  the same tracers, mask, nBins, distances and cut.
   *
   *  @param tracersSky as in reconstructLightcone.
   *  @param mask the survey's veto mask; borrowed for the call only.
   *  @param nBins as in reconstructLightcone.
   *  @param distances as in the overload with the sky area.
   *  @param cut as in reconstructLightcone, evaluated with the mask on
   *  every tracer.
   *  @param verbosity as in the overload with the sky area; the tracer line
   *  of the selection report is always written, since the mask applies.
   *  @param log as in the overload with the sky area.
   *
   *  @exception Error as the overload with the sky area, the counts in the
   *  message naming the mask as well as the cut.
   */
  MpsProfile measureMpsProfile (const std::vector<double>& tracersSky,
                                const Mask& mask,
                                unsigned nBins,
                                const DistanceTable& distances,
                                const RedshiftCut& cut = {},
                                Verbosity verbosity = Verbosity::Normal,
                                std::ostream* log = &std::clog);

  /**
   *  @brief The mps(d) profile the Cartesian overloads of
   *  reconstructLightcone measure, without reconstructing: from Cartesian
   *  coordinates, with the sky area given.
   *
   *  The steps of the Cartesian reconstructLightcone that concern the
   *  tracers, in the same order and with the same code: the checks of the
   *  array and of the cut, the distance cut, the minimum number of tracers
   *  kept, and the profile in distance with its per-bin threshold. The
   *  profile equals, bit for bit, Result::mpsProfile of the Cartesian
   *  reconstructLightcone with the same tracers, skyAreaDeg2, nBins and
   *  cut, whatever its randoms, config and table.
   *
   *  @param tracers as in the Cartesian reconstructLightcone.
   *  @param skyAreaDeg2 as in reconstructLightcone.
   *  @param nBins as in the Cartesian reconstructLightcone.
   *  @param cut as in the Cartesian reconstructLightcone, applied to the
   *  tracers.
   *  @param verbosity what the call reports. At Normal, when the cut has a
   *  finite bound, the tracer line of the selection report
   *  (SelectionCounts::message), then
   *
   *      otswap: measureMpsProfile: <B> bins over [<dmin>, <dmax>] in coordinate units
   *
   *  At Detailed, the mps(d) line of reconstructLightcone before the last.
   *  @param log where the report is written; null writes nothing.
   *
   *  @exception Error as the Cartesian reconstructLightcone, with the same
   *  messages, if the array is malformed, if a bound of the cut is NaN or
   *  min > max, if a cut with a finite bound keeps too few tracers, if
   *  every tracer kept lies at the same distance, or if a bin holds too
   *  few tracers.
   */
  MpsProfile measureMpsProfile (const std::vector<double>& tracers,
                                double skyAreaDeg2,
                                unsigned nBins,
                                const DistanceCut& cut = {},
                                Verbosity verbosity = Verbosity::Normal,
                                std::ostream* log = &std::clog);

  /**
   *  @brief The mps(d) profile of the Cartesian reconstructLightcone with a
   *  mask, without reconstructing.
   *
   *  As the Cartesian overload with the sky area, with the area taken from
   *  the mask, mask.skyAreaDeg2(), and the tracers on unobserved pixels
   *  left out as the Cartesian reconstructLightcone with a mask leaves them
   *  out; at Normal the tracer line of the selection report comes first.
   *  The profile equals, bit for bit, Result::mpsProfile of that overload
   *  with the same tracers, mask, nBins and cut. The cut is evaluated with
   *  the mask on every tracer.
   *
   *  @exception Error as the Cartesian overload with the sky area, and if
   *  too few tracers remain inside the mask; the counts in the message
   *  name the mask as well as the cut.
   */
  MpsProfile measureMpsProfile (const std::vector<double>& tracers,
                                const Mask& mask,
                                unsigned nBins,
                                const DistanceCut& cut = {},
                                Verbosity verbosity = Verbosity::Normal,
                                std::ostream* log = &std::clog);

  // ==========================================================================
  // Filtering
  // ==========================================================================

  /**
   *  @brief Reject the displacements whose path crosses masked sky.
   *
   *  The segment from the Eulerian position, the tracer (Result::tracers),
   *  to the Lagrangian one, the tracer plus its displacement, lies in the
   *  plane through the origin and both endpoints, so its projection on
   *  the sky is exactly the great circle arc between them. The arc is
   *  bisected adaptively: a sub-arc whose ends fall in the same pixel is
   *  not split further, nor one shorter than a millionth of the pixel
   *  size, nor one shorter than a pixel whose end pixels share an edge
   *  and have no common neighbour that could change a run: with both end
   *  pixels unobserved, none at all; otherwise, no unobserved one. The
   *  pixels of the midpoints are examined in order along the arc, and the
   *  unobserved pixels met one after another, with no observed pixel
   *  between them, form a run; a run is counted in distinct pixels. The
   *  pixels holding the two endpoints are exempt, and end a run as an
   *  observed pixel does.
   *
   *  A displacement is rejected when the longest run along its arc holds
   *  more than maxUnobservedPixelsCrossed pixels; the search stops as soon
   *  as a run exceeds it. Unobserved pixels in separate runs are not added
   *  together: an arc that crosses two small holes is judged by the larger
   *  one. With 0, the first unobserved pixel rejects. A displacement with
   *  an endpoint at the origin is kept, and one whose endpoints lie in
   *  exactly opposite directions, which defines no arc, is rejected.
   *
   *  A rejected displacement has its entry of valid cleared; entries
   *  already cleared stay cleared, so the filter composes with any other
   *  and applying it twice with the same mask and threshold changes
   *  nothing. The fields that follow from valid are then recomputed, as
   *  recomputeMeans does: validRealizations, meanDisplacement, lagrangian
   *  and lagrangianSky. On a result holding a mean (nRealizationsMean above
   *  1, one realization), validRealizations keeps the counts of the objects
   *  kept, and is 0 for those rejected.
   *
   *  Result::selection records the filter: maxUnobservedPixelsCrossed the
   *  threshold of this call, displacements the valid displacements before
   *  the first filter applied to the result, and displacementsCrossingMask
   *  those rejected by every filter so far.
   *
   *  @param result a lightcone result, updated in place; a box result is
   *  refused, so a Result built by hand must set Result::geometry.
   *  @param mask the mask whose unobserved pixels are counted.
   *  @param maxUnobservedPixelsCrossed the longest run of unobserved pixels
   *  tolerated along the arc, the endpoint pixels excluded.
   *
   *  @warning The run is counted in pixels, not as an angular length, so it
   *  depends on NSIDE: the same arc crosses about twice as many pixels when
   *  NSIDE doubles, and the same threshold is a different criterion.
   *
   *  A tracer flagged in Result::outsideRedshiftCut or Result::outsideMask
   *  is left as it is: its tracers and displacement rows may be NaN, and
   *  its entries of valid must be 0. Every other tracer and displacement
   *  must be finite.
   *
   *  @exception Error if the result is a box result, which has no sky to
   *  filter against; if the result is malformed (including a NaN in a
   *  tracer outside the flagged ones or in a valid displacement, a flagged
   *  tracer with a valid entry, a
   *  tracers array without 3 * nObjects entries, or a flag, tracersSky,
   *  lagrangian or lagrangianSky array whose size is neither 0 nor the one
   *  of its layout), if an entry of valid is
   *  neither 0 nor 1, or if the result was already filtered against a
   *  mask of a different NSIDE.
   */
  void rejectMaskCrossings (Result& result,
                            const Mask& mask,
                            unsigned maxUnobservedPixelsCrossed = 0);

  /**
   *  @brief Recompute every field of a result that follows from valid:
   *  validRealizations and meanDisplacement, then lagrangian and
   *  lagrangianSky.
   *
   *  meanDisplacement is summed over the valid realizations in realization
   *  order and divided by their count, as the reconstruction does, so a
   *  result whose displacement and valid were read back from a file gives
   *  the reconstruction's bits. lagrangian is tracers + meanDisplacement,
   *  left empty when tracers is; lagrangianSky is computed with
   *  result.distances, and left empty without a distance table.
   *
   *  @exception Error if displacement or valid has the wrong size, an entry
   *  of valid is neither 0 nor 1, or a flag, tracers or tracersSky array has
   *  a size that is neither 0 nor the one of its layout.
   */
  void recomputeMeans (Result& result);

}

#include "otswap/RSD.h"
#include "otswap/io.h"

#endif
