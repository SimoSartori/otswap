# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ============================================================================
# Example code: how to run the OT reconstruction on a lightcone in Cartesian
# coordinates, and correct it for redshift-space distortions
# ============================================================================

"""This example explains how to reconstruct a survey catalogue that is given
in comoving Cartesian coordinates, with the observer at the origin, instead
of right ascension, declination and redshift. Its tracers and randoms are
those of the lightcone and full_recons_lightcone examples, converted to
Cartesian coordinates with the same cosmology, and it does what both do: it
measures the mean particle separation alone, to choose the bins, then runs
the reconstruction, here with a cut on the distance, and the correction. It
first finds the distance and the direction of each tracer, with point_to_sky,
which is where the cut and the mask act.

Without redshifts the mean particle separation is measured as a function of
the comoving distance from the origin, in uniform distance bins, and the mask
and the optional distance cut act on the angles and the distance of each
object. The distance table is optional for the reconstruction, and only fills
the redshifts of the result; it is needed here for the correction, which
takes the growth rate f and the bias b(z) at the tracers' redshifts. The
example runs the reconstruction a second time without it, to show that
nothing else changes.

Usage: python lightcone_cartesian.py [data_dir]. The catalogues are read from
data_dir, by default the examples' data folder; the output is written to the
folder output/ next to this file, created if missing.
"""

import sys
from pathlib import Path

import numpy as np

import otswap


def main():

    # -------------------------------------------------------------
    # ------------ Locate the input and output folders ------------
    # -------------------------------------------------------------

    # The data folder is the examples' one, next to this example's folder,
    # unless a first command-line argument replaces it
    here = Path(__file__).resolve().parent
    data = Path(sys.argv[1]) if len(sys.argv) > 1 else here.parent / "data"
    output = here / "output"
    output.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------
    # ------------ Read the tracer and random catalogues ------------
    # ---------------------------------------------------------------

    # x, y, z for each object, in Mpc/h, the observer at the origin
    tracers = otswap.io.read(data / "lightcone_tracers_cartesian.dat", [0, 1, 2]).values
    randoms = otswap.io.read(data / "lightcone_randoms_cartesian.dat", [0, 1, 2]).values

    print(f"Tracers read: {len(tracers)}")
    print(f"Randoms read: {len(randoms)}, {len(randoms) / len(tracers):g} per tracer")

    # ---------------------------------------
    # ------------ Read the mask ------------
    # ---------------------------------------

    mask = otswap.Mask(str(data / "lightcone_mask.fits"))

    print(f"Mask: NSIDE {mask.nside}, {mask.sky_area_deg2:g} deg^2")

    # ---------------------------------------------------------------
    # ------------ Find the distances and the directions ------------
    # ---------------------------------------------------------------

    # point_to_sky gives, for each point, its distance from the origin, in
    # coordinate units, and its right ascension and declination, here in
    # degrees: the direction the mask is applied to, and the distance the
    # distance cut and the mps bins are taken on. It needs no distance
    # table. The range of the distances is where to choose the cut
    tracers_sky = otswap.point_to_sky(tracers, angle_unit="deg")

    print(f"Tracer distances: {tracers_sky[:, 0].min():g} to {tracers_sky[:, 0].max():g} Mpc/h")

    # One row per tracer, in input order: distance, RA, Dec
    otswap.io.write(output / "tracers_sky_lightcone_cartesian.dat",
                    [otswap.io.Column("distance", tracers_sky[:, 0], description="distance from the origin"),
                     otswap.io.Column("ra", tracers_sky[:, 1], description="right ascension", unit="deg"),
                     otswap.io.Column("dec", tracers_sky[:, 2], description="declination", unit="deg")])

    print("Written: output/tracers_sky_lightcone_cartesian.dat")

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # The table gives the linear growth rate f(z) of the correction, and the
    # redshifts of the result. Its distances are in Mpc/h, so the
    # coordinates must be in Mpc/h too: the unit of the coordinates is that
    # of any distance table used with them. It must cover the corrected
    # distances too: the default range, z = 0 to 10, leaves room for that
    distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)

    # --------------------------------------------
    # ------------ Set the parameters ------------
    # --------------------------------------------

    # A fixed seed, and the displacements crossing an unobserved pixel of the
    # mask rejected. 8 realizations, which the example randoms are enough
    # for, and a convergence threshold of 1e-2 keep the example fast; for an
    # analysis, start from the defaults, 30 and 1e-3 (the guide's section
    # "Choosing the parameters").
    #
    # The mean particle separation is measured in 30 distance bins of equal
    # width, each with enough tracers, over the distances of the tracers
    # kept.
    #
    # distance_cut: tracers and randoms whose comoving distance from the
    # origin lies outside this closed range, in coordinate units, are left
    # out of the reconstruction: the randoms are dropped, and the tracers
    # keep their row, flagged in Result.outside_distance_cut. The catalogue
    # spans 2126 to 2500 Mpc/h. The default, None, cuts nothing; a
    # (d_min, d_max) pair is accepted too
    config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, reject_crossings=True,
                           max_unobserved_pixels_crossed=0, verbosity="normal")
    n_bins = 30
    distance_cut = otswap.DistanceCut(2200.0, 2450.0)

    # --------------------------------------------------------------
    # ------------ Measure the mean particle separation ------------
    # --------------------------------------------------------------

    # measure_mps_profile_cartesian measures mps(d) from the tracers alone,
    # as the reconstruction does, from the same mask, bins and cut, and in a
    # fraction of its time: the profile is Result.mps_profile of the
    # reconstruction below, bit for bit. Use it to choose n_bins before a
    # reconstruction. It needs no distance table
    profile = otswap.measure_mps_profile_cartesian(tracers, mask=mask, n_bins=n_bins,
                                                   distance_cut=distance_cut)

    print(f"mps(d): {len(profile.mps)} bins, mps from {profile.mps.min():g} to {profile.mps.max():g} Mpc/h")

    # Without a mask, the keyword sky_area_deg2 gives the sky area, in square
    # degrees: no object is then left out on the sky and no crossing is
    # filtered, so the catalogues must lie within the footprint already.
    # These do, so with the mask's area the profile is the same
    from_area = otswap.measure_mps_profile_cartesian(tracers, sky_area_deg2=mask.sky_area_deg2, n_bins=n_bins,
                                                     distance_cut=distance_cut, verbosity="silent")

    print(f"The same profile from the sky area: {'yes' if np.array_equal(from_area.mps, profile.mps) else 'no'}")

    # One row per bin: its centre, in distance, its mps and its number of
    # tracers, with the range and the representative mps as keywords
    otswap.io.write_mps_profile(output / "mps_lightcone_cartesian.dat", profile)

    print("Written: output/mps_lightcone_cartesian.dat")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # distances fills the redshifts of the result. angle_unit is required:
    # it is the unit of the angles in the result's sky coordinates
    result = otswap.reconstruct_lightcone_cartesian(tracers, randoms, mask=mask, n_bins=n_bins,
                                                    distances=distances, angle_unit="deg",
                                                    distance_cut=distance_cut, config=config)

    # The profile is in distance, in coordinate units (Mpc/h here): the one
    # measured above
    same_profile = np.array_equal(result.mps_profile.mps, profile.mps)
    print(f"The profile of the result is the one measured: {'yes' if same_profile else 'no'}")

    # ----------------------------------------------------------------
    # ------------ Reconstruct without the distance table ------------
    # ----------------------------------------------------------------

    # Without a table, distances omitted, the reconstruction does the same
    # work, from the same seed, and gives the same displacements, valid flags
    # and profile, bit for bit; only the redshifts of the result are NaN, and
    # its files hold no distance table. NaN counts as equal to NaN here: the
    # tracers left out by the cut, and those without a valid realization,
    # have NaN means
    without_table = otswap.reconstruct_lightcone_cartesian(tracers, randoms, mask=mask, n_bins=n_bins,
                                                           angle_unit="deg", distance_cut=distance_cut,
                                                           config=config)
    same = (np.array_equal(without_table.displacement, result.displacement, equal_nan=True)
            and np.array_equal(without_table.valid, result.valid)
            and np.array_equal(without_table.mps_profile.mps, result.mps_profile.mps))
    print(f"Without the table: {'the same' if same else 'NOT the same'} displacements, valid flags and profile")

    # ---------------------------------------------
    # ------------ Read the bias table ------------
    # ---------------------------------------------

    bias = otswap.io.read_bias_table(data / "lightcone_bias.dat")

    # ----------------------------------------------------------------
    # ------------ Correct the redshift-space distortions ------------
    # ----------------------------------------------------------------

    # As in full_recons_lightcone: real_space_lightcone takes the result of
    # either function. The distance table is required here, even for a
    # result reconstructed without one: the growth rate f and the bias b(z)
    # are taken at the tracers' redshifts, which the table gives from their
    # distances.
    #
    # sigma, the width of the Gaussian average in coordinate units (Mpc/h
    # here): 10 Mpc/h is a starting value, and the value must be calibrated
    # on the sample (the guide's section "Choosing the parameters")
    correction = otswap.CorrectionConfig(weight_by_realizations=False, verbosity="normal")
    catalogue = otswap.real_space_lightcone(result, distances=distances, bias=bias, sigma=10.0,
                                            correction=correction)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    mean = result.mean_displacement
    length = np.linalg.norm(mean, axis=1)
    length = length[~np.isnan(length)]
    print(f"Mean length of the mean displacement over {len(length)} tracers: {length.mean():g} Mpc/h")
    shift = catalogue.shift[~np.isnan(catalogue.shift)]
    print(f"Shift along the line of sight over the {len(shift)} corrected tracers: mean {shift.mean():g}, "
          f"rms {np.sqrt(np.mean(shift ** 2)):g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # The mean displacement, with the mps(d) profile and the distance table
    # in twin files, and the corrected catalogue; the columns are those of
    # the lightcone and full_recons_lightcone examples, but the mean file
    # adds the column outsideDistanceCut, and the profile is binned in
    # distance: a column distance and the keywords DISTMIN and DISTMAX, in
    # coordinate units, instead of the redshift ones. The file of the result
    # without the table has NaN redshift columns (tracRed, lagrRed) and no
    # distance table twin. The writers are the C++ library's, so the files
    # are the C++ example's, byte for byte
    otswap.io.write_displacement_mean(output / "displacement_lightcone_cartesian.dat", result)
    otswap.io.write_displacement_mean(output / "displacement_lightcone_cartesian_no_table.dat", without_table)
    otswap.io.write_real_space_catalog(output / "reconstructed_catalogue_lightcone_cartesian.dat", catalogue)

    print("Written: output/displacement_lightcone_cartesian.dat")
    print("Written: output/displacement_lightcone_cartesian_no_table.dat")
    print("Written: output/reconstructed_catalogue_lightcone_cartesian.dat")


if __name__ == "__main__":
    main()
