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

    /// A pixel is allowed when its value exceeds this: any positive value,
    /// fractional ones included. 0, negative values, NaN and Healpix's
    /// UNSEEN are not. The sky area is a count of allowed pixels, not a sum
    /// of their values.
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

    /// check_coordinates for a sky array, which also refuses a
    /// declination outside [-pi/2, pi/2]; the right ascension may be any
    /// finite value.
    std::size_t check_sky (const std::vector<double>& sky,
                           const std::string& name,
                           bool allowEmpty = false);

    /// Checks that each flag array of @p result, outsideRedshiftCut and
    /// outsideMask, is empty or holds nObjects entries, and that tracers,
    /// tracersSky, lagrangian and lagrangianSky are each empty or hold
    /// 3 * nObjects.
    void check_flags (const Result& result);

    /// True when tracer @p i of @p result took no part in the
    /// reconstruction: flagged in outsideRedshiftCut or in outsideMask.
    /// The flags must have passed check_flags.
    inline bool excluded (const Result& result, const std::size_t i)
    {
      return (!result.outsideRedshiftCut.empty() && result.outsideRedshiftCut[i] != 0) ||
             (!result.outsideMask.empty() && result.outsideMask[i] != 0);
    }

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

    /// Set the fields that follow from valid: validRealizations and
    /// meanDisplacement (the count of valid realizations per object, and
    /// the mean over them, summed in realization order; NaN where there are
    /// none), lagrangian = tracers + meanDisplacement when tracers holds
    /// 3 * nObjects entries (empty otherwise), and lagrangianSky from
    /// lagrangian with result.distances in a lightcone result (empty
    /// otherwise). An array that already has its size is overwritten in
    /// place, so its data pointer does not change.
    void summarize (Result& result);

    /// The line a reconstruction ends its report with, '\n' included:
    /// "otswap: <name>: <R> realizations of <N> tracers in <T> s; <K>
    /// without a valid displacement".
    std::string reconstruction_line (const std::string& name, const Result& result);

    /// Sky coordinates of Cartesian positions, as toSky gives them, except
    /// that a distance outside the table gives a NaN redshift, the right
    /// ascension and declination kept, rather than an error. Written into
    /// @p sky, resized to the size of @p cartesian.
    void sky_of_positions (const std::vector<double>& cartesian,
                           const DistanceTable& distances,
                           std::vector<double>& sky);

    /// (V / N)^(1/3) of positions (3*N, finite, N > 0), V the volume of their
    /// bounding box: the mean separation of N objects filling it. 0 when the
    /// positions span no volume.
    double bounding_box_separation (const std::vector<double>& positions);

    /// Sky to Cartesian: x = distance cos(dec) cos(ra),
    /// y = distance cos(dec) sin(ra), z = distance sin(dec).
    void to_cartesian (double ra, double dec, double distance,
                       double& x, double& y, double& z);

    /// Cartesian to sky. Right ascension is returned in [0, 2*pi).
    void to_sky (double x, double y, double z,
                 double& ra, double& dec, double& distance);

    /// Right ascension folded into [0, 2*pi): normalize_ra may return 2 pi
    /// itself for a negative angle smaller in magnitude than half its
    /// rounding step, and -0 for -0, and fold_ra folds both to +0.
    double normalize_ra (double ra);
    double fold_ra (double ra);

    /// Build the profile of the tracers' redshifts, in uniform bins over the
    /// observed range, with its representative value. Throws when a bin is
    /// too thinly populated for its density to be meaningful.
    MpsProfile mps_profile (const std::vector<double>& tracersSky,
                            double skyAreaDeg2, unsigned nBins,
                            const DistanceTable& distances);

    /// Piecewise linear function through the nodes (x, y), x strictly
    /// increasing, extrapolated linearly beyond the terminal nodes along
    /// the terminal segment. A single node gives a constant.
    double profile_at (const std::vector<double>& x,
                       const std::vector<double>& y, double at);

    /// Count-weighted median of the profile nodes, used to size the grids.
    double representative (const MpsProfile& profile);

    /// Cell side of the grid neighbourAverage searches, for these positions
    /// (3*N, finite, N > 0): 4 bounding_box_separation(positions). A pure
    /// function of the positions. It is 0 when the positions span no
    /// volume, and the grid then refuses it.
    double neighbour_cell (const std::vector<double>& positions);

    /// neighbourAverage with the grid's cell side given; the public
    /// overloads call it with neighbour_cell(positions). The grid spans the
    /// bounding box of every position and holds the valid objects only.
    std::vector<double> neighbour_average (const std::vector<double>& positions,
                                           const std::vector<double>& values,
                                           const std::vector<unsigned>& validRealizations,
                                           double sigma,
                                           bool weightByRealizations,
                                           double cellSize,
                                           std::vector<unsigned>& nNeighbours,
                                           std::vector<unsigned>& nRealizations);

    /// Check a b(z) table: at least two nodes, equal lengths, redshifts
    /// finite and strictly increasing, bias finite and positive. The
    /// message starts with @p what and names the offending entry as
    /// @p entry followed by its index plus @p first.
    void check_bias_table (const std::vector<double>& redshift,
                           const std::vector<double>& bias,
                           const std::string& what,
                           const std::string& entry = "node",
                           std::size_t first = 0);

  }

}

#endif
