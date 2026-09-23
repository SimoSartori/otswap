/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
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
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace otswap {

  // ==========================================================================
  // Errors
  // ==========================================================================

  /// Base class of every exception thrown by otswap.
  class Error : public std::runtime_error {
  public:
    explicit Error (const std::string& what) : std::runtime_error(what) {}
  };

  // ==========================================================================
  // Coordinate layout
  // ==========================================================================
  //
  // Every coordinate array is flat and row-major, one object after the
  // other. Cartesian arrays hold 3 * nObjects entries, ordered x, y, z.
  // Sky arrays hold 3 * nObjects entries, ordered right ascension,
  // declination, redshift; angles are in radians.
  //
  // Nothing is copied into a nested container at any point: the layout is
  // the one the algorithm walks, and the one a numpy array maps onto
  // without a copy.

  // ==========================================================================
  // Configuration
  // ==========================================================================

  /**
   *  @brief Free parameters of the reconstruction.
   *
   *  The multipliers that govern the search radius, the neighbourhood
   *  size and the initial pairing are deliberately absent: they are
   *  internal constants, calibrated once, not knobs.
   */
  struct Config {

    /// Independent reconstructions to run. Each consumes nObjects randoms.
    unsigned nRealizations = 1;

    /// A sweep that changes no more than this fraction of pairs ends the
    /// loop. Every tracer is examined in every sweep.
    double convergence = 1.e-3;

    /// Seed of the generator. 0 draws one from std::random_device.
    unsigned seed = 0;

    /// Grid cell size, in units of the representative mean particle
    /// separation. Affects speed only, never the result.
    double cellSize = 4.;
  };

  // ==========================================================================
  // Result
  // ==========================================================================

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
    /// summed in realization order. NaN in all three components for an
    /// object with no valid realization. Flat, [object][xyz]. Size
    /// 3 * nObjects.
    std::vector<double> meanDisplacement;

    /// Valid realizations per object: the count of 1 entries of valid for
    /// that object, and the number of terms in its meanDisplacement.
    /// Size nObjects.
    std::vector<unsigned> validRealizations;

    /// NSIDE of the mask this result was last filtered against, or 0 when
    /// it has not been filtered. Set by rejectMaskCrossings, which refuses
    /// to filter the same result against a mask of a different NSIDE.
    int filteredNside = 0;
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
     *  @param nSamples sampling points; must be at least two
     *
     *  @exception Error if any parameter is out of range.
     */
    DistanceTable (double OmegaM, double h, double w0, double wa,
                   double zMin, double zMax, unsigned nSamples = 1000);

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

    /// Comoving distance at z. Throws outside the sampled range.
    double distanceAt (double z) const;

    /// Redshift at a comoving distance. Throws outside the sampled range.
    double redshiftAt (double distance) const;

    /// Growth rate f = dlnD/dlna at z. Throws outside the sampled range,
    /// and when the table carries no growth rate.
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

  /**
   *  @brief Healpix veto mask, read from a FITS file.
   *
   *  The mask is strictly binary: every pixel must be exactly 0 or 1, and
   *  the reader refuses the file otherwise. A pixel is allowed when its
   *  value is greater than zero. With a binary mask there is only one
   *  threshold to choose and skyAreaDeg2 is exact by construction, being
   *  a count of pixels rather than a sum of fractions.
   *
   *  The object is immutable once built.
   */
  class Mask {

  public:

    /**
     *  @param fitsFile Healpix map; NSIDE and ORDERING are taken from the
     *  header, in either RING or NESTED scheme.
     *
     *  @exception Error if the file cannot be opened, is not a valid
     *  Healpix map, or holds a pixel that is neither 0 nor 1; the message
     *  names the first such pixel and its value.
     */
    explicit Mask (const std::string& fitsFile);

    /// True when the direction falls in an allowed pixel. Radians.
    bool allows (double rightAscension, double declination) const;

    /// Sky area covered by the allowed pixels, in square degrees.
    double skyAreaDeg2 () const;

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
   *  @exception Error if any array has a size that is not a multiple of
   *  three, if the tracer array is empty, if the randoms are too few, if
   *  any coordinate is not finite, or if mps is not positive.
   *
   *  @return the displacement field; every entry of valid is 1,
   *  validRealizations is uniformly config.nRealizations, and
   *  meanDisplacement is the mean over all realizations.
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
   *  @param skyAreaDeg2 effective area of the survey. Take it from
   *  Mask::skyAreaDeg2 when a mask is available, otherwise supply the
   *  published value: a bounding box in right ascension and declination
   *  overestimates it for any footprint that is not rectangular, and the
   *  error propagates into the mean particle separation.
   *
   *  @param nBins redshift bins used to measure mps(z).
   *
   *  @exception Error if a bin holds too few tracers for its density to
   *  be meaningful; the message reports how many bins the catalog
   *  supports. Also if the profile extrapolates to a mps that is not
   *  positive at a tracer's redshift, if any array is malformed, if the
   *  randoms are too few, or if a redshift falls outside the distance
   *  table.
   */
  Result reconstructLightcone (const std::vector<double>& tracersSky,
                               const std::vector<double>& randomsSky,
                               double skyAreaDeg2,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config);

  /**
   *  @brief Reconstruct in lightcone geometry, with the Cartesian
   *  coordinates already computed.
   *
   *  Identical to the overload above, except that the conversion is
   *  skipped. The sky coordinates are still required: their redshifts
   *  drive mps(z), and they are written to the output table.
   *
   *  @param tracers Cartesian coordinates, 3 * nObjects entries.
   *  @param randoms Cartesian coordinates of the randoms.
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
                               double skyAreaDeg2,
                               unsigned nBins,
                               const DistanceTable& distances,
                               const Config& config);

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
   *  maxForbiddenPixels distinct unobserved pixels; the search stops as
   *  soon as the count exceeds it. With 0, the first one rejects. A
   *  displacement with an endpoint at the origin is kept, and one whose
   *  endpoints lie in exactly opposite directions, which defines no arc,
   *  is rejected.
   *
   *  A rejected displacement has its entry of valid cleared; entries
   *  already cleared stay cleared, so the filter composes with any other
   *  and applying it twice with the same mask and threshold changes
   *  nothing. validRealizations and meanDisplacement are then recomputed
   *  from valid.
   *
   *  @param maxForbiddenPixels distinct unobserved pixels tolerated along
   *  the arc, the endpoint pixels excluded.
   *
   *  @warning The count is of pixels, not an angular length, so it depends
   *  on NSIDE: the same arc crosses about twice as many pixels when NSIDE
   *  doubles, and the same threshold is a different criterion.
   *
   *  @exception Error if the result is malformed, if an entry of valid is
   *  neither 0 nor 1, or if the result was already filtered against a
   *  mask of a different NSIDE.
   */
  void rejectMaskCrossings (Result& result,
                            const Mask& mask,
                            unsigned maxForbiddenPixels = 0);

  // ==========================================================================
  // Tables
  // ==========================================================================

  namespace io {

    /// One column of an output table.
    struct Column {
      std::string name;
      char        type = 'D';        ///< 'D' double, 'J' integer
      std::string description;       ///< written as a header comment
      std::vector<double> data;
    };

    /// A table read from file.
    struct Table {
      std::size_t nRows = 0;
      std::size_t nColumns = 0;
      std::vector<double> values;    ///< flat, row-major [row][column]
    };

    /**
     *  @brief Read selected columns of a table.
     *
     *  The format follows the extension: .fits, .fit and .fits.gz are
     *  read as FITS, anything else as ASCII. In FITS the columns are
     *  named, case insensitively; in ASCII they are 0-based indices
     *  given as strings. Vector-valued FITS columns are rejected.
     *
     *  @exception Error if the file is missing, a column is absent, a
     *  requested index is not an integer, a row is too short, or an
     *  undefined value is met.
     */
    Table read (const std::string& file,
                const std::vector<std::string>& columns,
                char delimiter = ' ', char comment = '#');

    /**
     *  @brief Write a table. An existing file is overwritten.
     *
     *  @exception Error if the columns differ in length, or the file
     *  cannot be written.
     */
    void write (const std::string& file, const std::vector<Column>& columns);

    /**
     *  @brief Write a table whose rows are produced on demand.
     *
     *  The data member of each column is ignored; fillRow is called once
     *  per row with a buffer of one entry per column. Nothing is held in
     *  memory, so this is the variant for catalog-sized output.
     */
    void write (const std::string& file, const std::vector<Column>& columns,
                std::size_t nRows,
                const std::function<void(std::size_t, std::vector<double>&)>& fillRow);

  }

}

#endif
