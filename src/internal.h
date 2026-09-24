/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * Distributed under the BSD 3-Clause License. See LICENSE.         *
 ********************************************************************/

/**
 *  @file src/internal.h
 *
 *  @brief Declarations shared by the otswap translation units. Nothing
 *  here is part of the public API.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef OTSWAP_INTERNAL_H
#define OTSWAP_INTERNAL_H

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "otswap/OT.h"

namespace otswap {

  namespace internal {

    /// A pixel is allowed when its value exceeds this. The mask is
    /// validated as strictly binary on the way in, so the only threshold
    /// that can separate 0 from 1 is zero, and the sky area is a pixel
    /// count rather than a sum of fractions.
    constexpr double kMaskAllowedAbove = 0.;

    /// Uniform integer in [0, bound), by Lemire's multiply-and-reject method
    /// on the engine's 32-bit output. Unlike std::uniform_int_distribution,
    /// whose algorithm the standard leaves open, it yields the same sequence
    /// on every toolchain.
    unsigned uniform_int (std::mt19937& rng, unsigned bound);

    /// Uniform double in [lo, hi), from 53 random bits of two engine
    /// outputs. The same sequence on every toolchain, which
    /// std::uniform_real_distribution does not guarantee.
    double uniform_real (std::mt19937& rng, double lo, double hi);

    /// Fisher-Yates shuffle over uniform_int, from the last element down.
    /// The same permutation on every toolchain, which std::shuffle does not
    /// guarantee.
    template <typename T>
    void shuffle (std::vector<T>& v, std::mt19937& rng)
    {
      for (std::size_t i = v.size(); i > 1; --i)
        std::swap(v[i-1], v[uniform_int(rng, (unsigned)i)]);
    }

    /// The independent random streams of a run.
    enum class Stream : std::uint32_t {
      Supply  = 1,  ///< the shuffle that assigns given randoms to realizations
      Randoms = 2,  ///< the uniform draw of randoms, when none are given
      Seeding = 3,  ///< the seeding pass: its opening shuffles and the pairs it draws
      Swap    = 4   ///< the swap loop: its visiting order and swap partners
    };

    /// Engine for stream @p id of realization @p realization.
    ///
    /// The key is splitmix64(splitmix64(seed << 32 | realization) ^ id),
    /// and its two 32-bit halves fill the engine's state through
    /// std::seed_seq. Keys differ for any two triples that share the
    /// stream, and for any two that share seed and realization. Every step
    /// is specified exactly by the standard, so an engine yields the same
    /// sequence on every toolchain, and it depends on nothing but the
    /// triple, so no draw depends on the thread count.
    std::mt19937 stream (unsigned seed, unsigned realization, Stream id);

    /// Number of objects in a flat coordinate array, after checking that
    /// its size is a multiple of three and every entry is finite.
    std::size_t check_coordinates (const std::vector<double>& a,
                                   const std::string& name,
                                   bool allowEmpty = false);

    /// Check the fields of Config, and resolve seed 0 to a drawn seed.
    unsigned check_config (const Config& config);

    /// Check that the randoms can serve nRealizations disjoint subsets.
    void check_random_supply (std::size_t nRandoms, std::size_t nObjects,
                              unsigned nRealizations);

    /// The smallest catalog accepted: one more than the neighbours held
    /// per tracer. Defined in Reconstruct.cpp.
    std::size_t min_objects ();

    /// One reconstruction, shared by both geometries.
    ///
    /// @param tracers flat Cartesian, 3*nObjects
    /// @param randoms flat Cartesian, at least 3*nRealizations*nObjects,
    ///   or empty to draw them in the tracers' bounding box
    /// @param mps per-object mean particle separation, nObjects entries
    /// @param representativeMps the value the grids are sized from
    /// @param seed already resolved; never 0
    /// @param sweepCosts if given, receives the mean squared cost after
    ///   each sweep, per realization; its length is the sweep count
    ///
    /// The grids span the union of the tracers' and the given randoms'
    /// extents, so a random may lie outside the tracers' bounding box.
    Result reconstruct (const std::vector<double>& tracers,
                        const std::vector<double>& randoms,
                        const std::vector<double>& mps,
                        double representativeMps,
                        const Config& config,
                        unsigned seed,
                        std::vector<std::vector<double>>* sweepCosts = nullptr);

    /// Set validRealizations and meanDisplacement from displacement and
    /// valid: the count of valid realizations per object, and the mean
    /// over them, summed in realization order; NaN where there are none.
    void summarize (Result& result);

    /// Sky to Cartesian: x = distance cos(dec) cos(ra),
    /// y = distance cos(dec) sin(ra), z = distance sin(dec).
    void to_cartesian (double ra, double dec, double distance,
                       double& x, double& y, double& z);

    /// Cartesian to sky. Right ascension is returned in [0, 2*pi).
    void to_sky (double x, double y, double z,
                 double& ra, double& dec, double& distance);

    /// Right ascension folded into [0, 2*pi).
    double normalize_ra (double ra);

    /// Nodes of the mean particle separation profile, measured from the
    /// tracers' redshifts in uniform bins over the observed range.
    struct MpsProfile {
      std::vector<double>   redshift;   ///< bin centres
      std::vector<double>   mps;        ///< mps at each bin centre
      std::vector<unsigned> count;      ///< tracers in each bin
    };

    /// Build the profile. Throws when a bin is too thinly populated for
    /// its density to be meaningful.
    MpsProfile mps_profile (const std::vector<double>& tracersSky,
                            double skyAreaDeg2, unsigned nBins,
                            const DistanceTable& distances);

    /// Piecewise linear function through the nodes (x, y), x strictly
    /// increasing, extrapolated linearly beyond the terminal nodes along
    /// the terminal segment. A single node gives a constant.
    double profile_at (const std::vector<double>& x,
                       const std::vector<double>& y, double at);

    /// Evaluate the profile at @p z: linear between nodes, extrapolated
    /// linearly beyond the terminal ones along the terminal segment, and
    /// constant when there is a single node. Throws when the value is not
    /// positive, which only an extrapolation can produce.
    double mps_at (const MpsProfile& profile, double z);

    /// Count-weighted median of the profile nodes, used to size the grids.
    double representative (const MpsProfile& profile);

  }

}

#endif
