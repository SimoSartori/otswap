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
  /// from inside a parallel region too. The neighbour search of the
  /// redshift-space correction raises meshsearch::Error, a
  /// std::runtime_error, for positions that span no volume.
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
   *  Mpc/h, the distance table's unit. Sky arrays hold 3 * nObjects entries,
   *  ordered right ascension, declination, redshift; angles are in radians.
   *
   *  Nothing is copied into a nested container at any point: the layout is
   *  the one the algorithm walks, and the one a numpy array maps onto
   *  without a copy. The arrays of a Result and of a RealSpaceCatalog:
   *
   *      [realization][object][3]  Result::displacement, Result::matchedRandom
   *      [realization][object]     Result::valid (uint8, 0 or 1)
   *      [object][3]               Result::tracers, tracersSky, meanDisplacement,
   *                                lagrangian, lagrangianSky;
   *                                RealSpaceCatalog::sky, cartesian
   *      [object]                  Result::validRealizations, outsideRedshiftCut,
   *                                outsideMask (uint8);
   *                                RealSpaceCatalog::shift, factor, status,
   *                                validRealizations, nNeighbours,
   *                                nRealizationsAveraged
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
   */
  enum class Verbosity {
    /// Nothing.
    Silent,

    /// The selection report of a lightcone, when a selection applied
    /// (SelectionCounts::message), and the extrapolation of b(z), when it
    /// happened; then one line per call with its time and what was lost:
    ///
    ///     otswap: reconstructLightcone: <R> realizations of <N> tracers in <T> s; <K> without a valid displacement
    ///     otswap: realSpaceLightcone: corrected <C> of <N> tracers in <T> s; <U> left uncorrected
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
    unsigned nRealizations = 1;

    /// The sweeps stop once the fraction of successful swaps in a sweep,
    /// swaps per tracer visited, no longer exceeds this threshold. Every
    /// tracer is examined in every sweep.
    double convergence = 1.e-3;

    /// Seed of the generator. 0 draws one from std::random_device.
    unsigned seed = 0;

    /// Grid cell size, in units of the representative mean particle
    /// separation. Affects speed only, never the result.
    double cellSize = 4.;

    /// With a Mask given to reconstructLightcone, apply rejectMaskCrossings
    /// to the result, with maxUnobservedPixelsCrossed, before returning it.
    /// Read by the mask overloads of reconstructLightcone only.
    bool rejectCrossings = true;

    /// The threshold of that call; see rejectMaskCrossings. Read by the
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
   *  @brief The objects a lightcone reconstruction left out, and why, and
   *  the displacements its mask filter rejected.
   *
   *  The tracers removed are tracersOutsideRedshiftCut + tracersOutsideMask
   *  - tracersOutsideBoth, and the same for the randoms. All counts are 0,
   *  and nothing is applied, in a box result.
   */
  struct SelectionCounts {

    std::optional<RedshiftCut> redshiftCut;  ///< the cut, when a bound was finite
    bool maskApplied = false;                ///< a Mask was given

    std::size_t tracers = 0;                    ///< tracers given
    std::size_t tracersOutsideRedshiftCut = 0;  ///< of them, outside the cut
    std::size_t tracersOutsideMask = 0;         ///< of them, on an unobserved pixel
    std::size_t tracersOutsideBoth = 0;         ///< counted in both of the above

    std::size_t randoms = 0;                    ///< randoms given
    std::size_t randomsOutsideRedshiftCut = 0;  ///< of them, outside the cut
    std::size_t randomsOutsideMask = 0;         ///< of them, on an unobserved pixel
    std::size_t randomsOutsideBoth = 0;         ///< counted in both of the above

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
     *      otswap: rejected <r> of <D> displacements crossing more than <L> unobserved pixels
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
   *  redshift, measured by reconstructLightcone.
   *
   *  The tracers kept are binned in nBins uniform redshift bins over
   *  [redshiftMin, redshiftMax], and mps = (N / V)^(-1/3) in each, V the
   *  shell volume implied by the sky area. The bin values are the nodes of
   *  a piecewise linear profile, extrapolated linearly beyond the terminal
   *  ones; at() evaluates it. Empty in a box result.
   */
  struct MpsProfile {
    std::vector<double>      redshift;   ///< bin centres
    std::vector<double>      mps;        ///< Mpc/h at each centre
    std::vector<std::size_t> count;      ///< tracers in each bin

    /// The binned range: the smallest and largest redshift of the tracers
    /// kept.
    double redshiftMin = std::numeric_limits<double>::quiet_NaN();
    double redshiftMax = std::numeric_limits<double>::quiet_NaN();

    /// The count-weighted median of the nodes, which sizes the grids of the
    /// reconstruction.
    double representative = std::numeric_limits<double>::quiet_NaN();

    /**
     *  @brief The mps at redshift z, in Mpc/h: linear between the nodes,
     *  extrapolated linearly beyond the terminal ones along the terminal
     *  segment, constant with a single node. This is the value the
     *  reconstruction used at each tracer.
     *
     *  @exception Error if the profile is empty, or if the value is not
     *  positive, which only an extrapolation can give.
     */
    double at (double z) const;
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

    /// The geometry of the reconstruction, which the writers and the
    /// redshift-space correction follow.
    Geometry geometry = Geometry::Box;

    /// The configuration the reconstruction ran with. seed is the one used:
    /// the drawn one when 0 was given, so that the result can be
    /// reproduced.
    Config config;

    /// Cartesian position of each tracer, the start of its displacements:
    /// the array given to reconstructBox or to a Cartesian overload of
    /// reconstructLightcone, or toCartesian of the sky array. A tracer left
    /// out by the redshift cut or the mask has its conversion when its
    /// redshift lies in the distance table, NaN otherwise. Flat,
    /// [object][xyz]; size 3 * nObjects.
    std::vector<double> tracers;

    /// Sky coordinates of each tracer, as given to reconstructLightcone:
    /// right ascension, declination, redshift; radians. Flat, [object][3];
    /// size 3 * nObjects in a lightcone result, empty in a box result.
    std::vector<double> tracersSky;

    /// Displacements. Flat, [realization][object][xyz].
    /// Size 3 * nRealizations * nObjects.
    std::vector<double> displacement;

    /// Cartesian coordinates of the random matched to each object.
    /// Same layout and size as displacement.
    std::vector<double> matchedRandom;

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
    /// that object, and the number of terms in its meanDisplacement.
    /// Size nObjects.
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
    /// A cut tracer has NaN in every component of its displacement and
    /// matchedRandom rows, valid 0 in every realization, validRealizations
    /// 0 and meanDisplacement NaN. No filter changes the flag. A Result
    /// built by hand may leave the field empty, which every function reads
    /// as all 0.
    std::vector<std::uint8_t> outsideRedshiftCut;

    /// Whether each tracer fell on an unobserved pixel of the Mask given to
    /// reconstructLightcone, and so took no part in the reconstruction: 1
    /// for such a tracer, 0 otherwise. Flat, [object]. Size nObjects, all 0
    /// when no mask was given and in a box result. Independent of
    /// outsideRedshiftCut: a tracer may carry both flags. A flagged
    /// tracer's rows are as for outsideRedshiftCut, and no filter changes
    /// the flag. An empty field is read as all 0.
    std::vector<std::uint8_t> outsideMask;

    /// What reconstructLightcone left out, and why, and what the mask
    /// filters rejected; updated by every rejectMaskCrossings.
    SelectionCounts selection;

    /// Mean Lagrangian position of each tracer, tracers + meanDisplacement.
    /// NaN rows where meanDisplacement is NaN. Flat, [object][xyz]; size
    /// 3 * nObjects. Recomputed with meanDisplacement by rejectMaskCrossings
    /// and recomputeMeans.
    std::vector<double> lagrangian;

    /// Sky coordinates of each tracer's mean Lagrangian position, as toSky
    /// computes them from lagrangian with the reconstruction's distance
    /// table: right ascension in [0, 2 pi), declination, redshift; radians.
    /// Flat, [object][3]; size 3 * nObjects in a lightcone result, empty in a
    /// box result. Recomputed with lagrangian.
    ///
    /// NaN rows where lagrangian is NaN: the tracers flagged in
    /// outsideRedshiftCut or outsideMask and those with no valid
    /// realization. A position whose distance falls outside the table,
    /// which a table that does not start at distance 0 allows, has a NaN
    /// redshift and keeps its right ascension and declination.
    std::vector<double> lagrangianSky;

    /// Box: the mean particle separation used, given or computed, Mpc/h.
    /// NaN in a lightcone result, whose separation is mpsProfile.
    double mps = std::numeric_limits<double>::quiet_NaN();

    /// Lightcone: the mean particle separation profile measured on the
    /// tracers kept. Empty in a box result.
    MpsProfile mpsProfile;

    /// Lightcone: the distance table of the reconstruction, shared by the
    /// copies of this result, with which rejectMaskCrossings and
    /// recomputeMeans recompute lagrangianSky. Null in a box result.
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
     *  Non-flat cosmologies are served by the table constructor below.
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
     *  @param w0 dark energy equation of state, constant term
     *  @param wa dark energy equation of state, evolution term
     *  @param zMin lower end of the sampled range; must be non-negative
     *  @param zMax upper end; must exceed zMin
     *  @param nSamples sampling points, uniform in redshift; must be at
     *  least two
     *
     *  @exception Error if any parameter is out of range.
     */
    DistanceTable (double OmegaM, double h, double w0, double wa,
                   double zMin = 0., double zMax = 10., unsigned nSamples = 50000);

    /**
     *  @brief Build from a table the caller has computed.
     *
     *  @param redshift strictly increasing, at least two entries
     *  @param distance strictly increasing, same length, comoving Mpc/h
     *
     *  @param growthRate f = dlnD/dlna on the same grid, or empty. When
     *  empty the table carries no growth rate and growthRateAt throws.
     *  Need not be monotonic; every entry must be finite and positive.
     *
     *  @exception Error if either of the first two arrays is not strictly
     *  increasing, if the lengths differ, or if fewer than two points are
     *  given.
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

    /// The filter reads the pixels directly, so that pixel identity and
    /// the allowed test have one source rather than two.
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
   *  empty, config.nRealizations * nObjects randoms are drawn uniformly
   *  in the bounding box of the tracers. When given, the array must hold
   *  at least config.nRealizations * nObjects of them, since each
   *  realization consumes a disjoint subset; they may lie outside the
   *  tracers' bounding box.
   *
   *  @param mps mean particle separation; must be positive.
   *
   *  @param config the parameters of the run; rejectCrossings and
   *  maxUnobservedPixelsCrossed are not read.
   *
   *  @exception Error if any array has a size that is not a multiple of
   *  three, if the tracer array is empty, if the randoms are too few, if
   *  any coordinate is not finite, or if mps is not positive.
   *
   *  @return the displacement field; every entry of valid is 1,
   *  validRealizations is uniformly config.nRealizations, and
   *  meanDisplacement is the mean over all realizations. Result::mps is
   *  mps.
   *
   *  @note The box is not periodic: nothing flows through its faces, so
   *  modes on the scale of the box itself cannot be reconstructed.
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
                         const Config& config);

  /**
   *  @brief Reconstruct in box geometry, with the mean particle separation
   *  of the tracers: mps = (V / N)^(1/3), V the volume of their bounding
   *  box, the box that drawn randoms fill, and N their number. The value
   *  used is Result::mps.
   *
   *  @exception Error as above, and if the tracers span no volume (all of
   *  them in one plane, on one line or at one point).
   */
  Result reconstructBox (const std::vector<double>& tracers,
                         const std::vector<double>& randoms,
                         const Config& config);

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
   *  @return Cartesian coordinates in the distance table's unit, Mpc/h,
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
   *  @param cartesian 3 * nObjects entries, ordered x, y, z, in Mpc/h. May
   *  be empty. A row holding a NaN gives a NaN row.
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
   *  @param randomsSky sky coordinates of the randoms; at least
   *  config.nRealizations * nObjects of them. They may lie outside the
   *  Cartesian bounding box of the tracers.
   *
   *  @param skyAreaDeg2 effective area of the survey. When a mask is
   *  available, pass the mask instead, to the overload below, which takes
   *  the area from it; otherwise supply the published value: a bounding box in right ascension and declination
   *  overestimates it for any footprint that is not rectangular, and the
   *  error propagates into the mean particle separation.
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
   *  inside it; the message gives the counts kept and dropped.
   */
  Result reconstructLightcone (const std::vector<double>& tracersSky,
                               const std::vector<double>& randomsSky,
                               double skyAreaDeg2,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config,
                               const RedshiftCut& cut = {});

  /**
   *  @brief Reconstruct in lightcone geometry, with the Cartesian
   *  coordinates already computed.
   *
   *  Identical to the overload above, except that the conversion is
   *  skipped. The sky coordinates are still required: the tracers'
   *  redshifts drive mps(z), and the sky random array must describe the
   *  same objects as the Cartesian one.
   *
   *  @param tracers Cartesian coordinates, 3 * nObjects entries.
   *  @param randoms Cartesian coordinates of the randoms.
   *  @param tracersSky as above.
   *  @param randomsSky as above.
   *  @param skyAreaDeg2 as above.
   *  @param nBins as above.
   *  @param distances as above.
   *  @param config as above.
   *
   *  @param cut as above, decided on the sky redshifts; the same rows are
   *  dropped from the Cartesian arrays.
   *
   *  @exception Error as above, and if the Cartesian and sky arrays
   *  describe a different number of objects.
   *
   *  @warning No check verifies that the Cartesian coordinates agree with
   *  the sky ones under the given cosmology; disagreement is silent.
   *  Arrays made with toCartesian and the same distance table agree
   *  exactly, and give the same result as the sky overload.
   */
  Result reconstructLightcone (const std::vector<double>& tracers,
                               const std::vector<double>& randoms,
                               const std::vector<double>& tracersSky,
                               const std::vector<double>& randomsSky,
                               double skyAreaDeg2,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config,
                               const RedshiftCut& cut = {});

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
   *  @exception Error as the overload with the sky area, the counts in the
   *  message naming the mask as well as the cut.
   */
  Result reconstructLightcone (const std::vector<double>& tracersSky,
                               const std::vector<double>& randomsSky,
                               const Mask& mask,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config,
                               const RedshiftCut& cut = {});

  /**
   *  @brief Reconstruct in lightcone geometry, with the Cartesian
   *  coordinates already computed, within the observed pixels of a mask.
   *
   *  The overload above with the conversion skipped, as for the overloads
   *  with the sky area. The mask is applied on the sky coordinates; the
   *  same rows are dropped from the Cartesian arrays.
   *
   *  @exception Error as above, and if the Cartesian and sky arrays
   *  describe a different number of objects.
   *
   *  @warning No check verifies that the Cartesian coordinates agree with
   *  the sky ones under the given cosmology; disagreement is silent.
   */
  Result reconstructLightcone (const std::vector<double>& tracers,
                               const std::vector<double>& randoms,
                               const std::vector<double>& tracersSky,
                               const std::vector<double>& randomsSky,
                               const Mask& mask,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config,
                               const RedshiftCut& cut = {});

  // ==========================================================================
  // Filtering
  // ==========================================================================

  /**
   *  @brief Reject the displacements whose path crosses masked sky.
   *
   *  The segment from the Eulerian to the Lagrangian position lies in the
   *  plane through the origin and both endpoints, so its projection on
   *  the sky is exactly the great circle arc between them. The arc is
   *  bisected adaptively: a sub-arc whose ends fall in the same pixel is
   *  not split further, nor one shorter than a millionth of the pixel
   *  size, nor one shorter than a pixel whose end pixels share an edge
   *  and have no unobserved pixel as a common neighbour. The pixel of
   *  every midpoint is examined, and the distinct unobserved pixels met
   *  are counted. The pixels holding the two endpoints are exempt.
   *
   *  A displacement is rejected when its arc crosses more than
   *  maxUnobservedPixelsCrossed distinct unobserved pixels; the search stops as
   *  soon as the count exceeds it. With 0, the first one rejects. A
   *  displacement with an endpoint at the origin is kept, and one whose
   *  endpoints lie in exactly opposite directions, which defines no arc,
   *  is rejected.
   *
   *  A rejected displacement has its entry of valid cleared; entries
   *  already cleared stay cleared, so the filter composes with any other
   *  and applying it twice with the same mask and threshold changes
   *  nothing. The fields that follow from valid are then recomputed, as
   *  recomputeMeans does: validRealizations, meanDisplacement, lagrangian
   *  and lagrangianSky.
   *
   *  Result::selection records the filter: maxUnobservedPixelsCrossed the
   *  threshold of this call, displacements the valid displacements before
   *  the first filter applied to the result, and displacementsCrossingMask
   *  those rejected by every filter so far.
   *
   *  @param result a lightcone result, updated in place.
   *  @param mask the mask whose unobserved pixels are counted.
   *  @param maxUnobservedPixelsCrossed distinct unobserved pixels tolerated
   *  along the arc, the endpoint pixels excluded.
   *
   *  @warning The count is of pixels, not an angular length, so it depends
   *  on NSIDE: the same arc crosses about twice as many pixels when NSIDE
   *  doubles, and the same threshold is a different criterion.
   *
   *  A tracer flagged in Result::outsideRedshiftCut or Result::outsideMask
   *  is left as it is: its displacement and matchedRandom rows may be NaN,
   *  and its entries of valid must be 0. Every other displacement must be
   *  finite.
   *
   *  @exception Error if the result is malformed (including a NaN outside
   *  the flagged tracers' rows, a flagged tracer with a valid entry, or a
   *  flag, tracers, tracersSky, lagrangian or lagrangianSky array whose
   *  size is neither 0 nor the one of its layout), if an entry of valid is
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
