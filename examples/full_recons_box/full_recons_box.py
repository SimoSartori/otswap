# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# =====================================================================
# Example code: how to move a box catalogue from redshift to real space
# =====================================================================

"""This example explains how to correct a catalogue in box geometry for
redshift-space distortions with the OT reconstruction. The tracers are given
in redshift space, distorted along the y axis (axis 1), the line of sight.
The example runs the reconstruction of the box example on them, then moves
each tracer along that axis to the position it would have in real space.

The correction of tracer i is a shift along the axis,

    s_i = f(z) / (b + 3 f(z) / 5) * <Psi_axis>_i ,

with f the linear growth rate at the box's redshift z, b the tracers' linear
bias, and <Psi_axis>_i the component along the axis of the reconstructed
mean displacement, averaged with a Gaussian of width sigma over the tracers
around i.

Without the periodic option, the box is not periodic, here or in the
reconstruction: nothing flows through its faces, so modes on the scale of
the box itself are not reconstructed, and the corrected positions are not
wrapped. A box cut from a periodic simulation is treated as a plain
sub-volume; a whole simulation box with periodic boundary conditions is
reconstructed, and corrected, with the periodic option: see the
periodic_box example.

Usage: python full_recons_box.py [data_dir]. The catalogues are read from
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

    # Comoving Cartesian coordinates in Mpc/h, shape (N, 3). The tracers are
    # in redshift space, with the line of sight along y; the randoms are
    # uniform in the volume
    tracers = otswap.io.read(data / "box_halos_redshift_space.dat", [0, 1, 2]).values
    randoms = otswap.io.read(data / "box_randoms.dat", [0, 1, 2]).values

    print(f"Tracers: {len(tracers)}")
    print(f"Randoms: {len(randoms)}, {len(randoms) / len(tracers):g} per tracer")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # As in the box example, but without box limits: a fixed seed, and the
    # mean particle separation of the bounding box of the tracers and
    # randoms, computed by otswap and printed with "detailed". 8
    # realizations, which the example randoms are enough for, and a
    # convergence threshold of 1e-2 keep the example fast; for an analysis,
    # start from the defaults, 30 and 1e-3 (the guide's section "Choosing
    # the parameters")
    config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, verbosity="detailed")
    result = otswap.reconstruct_box(tracers, randoms, config=config)

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # Only the linear growth rate f at the box's redshift is needed. The
    # table computes it from Omega_m, h, w0 and wa of a flat cosmology
    distances = otswap.DistanceTable.flat(0.3186, 0.67, w0=-1.0, wa=0.0)

    # ----------------------------------------------------------------
    # ------------ Correct the redshift-space distortions ------------
    # ----------------------------------------------------------------

    # The box's redshift, the tracers' linear bias, and the line of sight:
    # axis 0, 1 or 2 for x, y or z
    redshift, bias, axis = 1.0, 1.91255, 1

    # sigma, the width of the Gaussian average in Mpc/h: 10 Mpc/h is a
    # starting value, and the value must be calibrated on the sample (the
    # guide's section "Choosing the parameters"). The options are an
    # otswap.CorrectionConfig: with weight_by_realizations=True, each
    # neighbour would be weighted by its number of valid realizations.
    #
    # The correction reads the tracers from the result. In a box without
    # limits every tracer is valid and counts among its own neighbours, so
    # every one is corrected; catalogue.status records what was done with
    # each
    correction = otswap.CorrectionConfig(weight_by_realizations=False)
    catalogue = otswap.real_space_box(result, axis=axis, redshift=redshift, distances=distances,
                                      bias=bias, sigma=10.0, correction=correction)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # The single factor f/(b + 3f/5) of the box, and the shift applied along
    # the line of sight to each corrected tracer, both kept in the catalogue
    print(f"RSD factor f/(b + 3f/5): {catalogue.factor[0]:g}")
    shift = catalogue.shift[~np.isnan(catalogue.shift)]
    print(f"Shift along the line of sight over the {len(shift)} corrected tracers: mean {shift.mean():g}, "
          f"rms {np.sqrt(np.mean(shift ** 2)):g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: the corrected position,
    # NaN for an uncorrected tracer, then the number of valid realizations,
    # the two diagnostics of the average and the status of the tracer. The
    # writer is the C++ library's, so the file is the C++ example's, byte
    # for byte
    otswap.io.write_real_space_catalog(output / "reconstructed_catalogue_box.dat", catalogue)

    print("Written: output/reconstructed_catalogue_box.dat")

    # -----------------------------------------------------------
    # ------------ Compare with the real-space halos ------------
    # -----------------------------------------------------------

    # The file holds the corrected position in the columns tracX, tracY and
    # tracZ, which otswap.io.read finds by name in a file otswap wrote. The
    # example data give the same halos in real space too, row by row, so the
    # correction can be checked here, as survey data do not allow: along the
    # line of sight, the offset of each halo from its real-space position,
    # before and after the correction
    corrected = otswap.io.read(output / "reconstructed_catalogue_box.dat", ["tracX", "tracY", "tracZ"]).values
    real = otswap.io.read(data / "box_halos_real_space.dat", [0, 1, 2]).values

    kept = ~np.isnan(corrected[:, axis])
    before = tracers[kept, axis] - real[kept, axis]
    after = corrected[kept, axis] - real[kept, axis]
    print(f"rms offset from the real-space position along the line of sight, over {np.count_nonzero(kept)} "
          f"halos: {np.sqrt(np.mean(before ** 2)):g} Mpc/h before the correction, "
          f"{np.sqrt(np.mean(after ** 2)):g} Mpc/h after")


if __name__ == "__main__":
    main()
