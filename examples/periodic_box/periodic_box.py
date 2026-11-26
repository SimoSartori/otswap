# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ======================================================================
# Example code: how to run the OT reconstruction and the redshift-space
# correction in a periodic box
# ======================================================================

"""This example explains how to reconstruct the displacement field of halos
in redshift space in a cube with periodic boundary conditions, such as a
simulation box, and how to move them back to real space with it. The box is
given by its limits; every distance is then the minimum image, so no face is
a boundary: a tracer near a face may be matched to a random near the
opposite one, as anywhere else in the box, and the correction averages the
displacements of neighbours across the faces.

Each side of a periodic box must be an integer multiple of the shortest one,
so a cube, or a cuboid such as 500 x 500 x 1000, is accepted. Tracers
outside the limits are wrapped into the box, by a whole number of periods
per axis, and reported; the displacements point along the shortest way,
within half a period on each axis.

The halos are at redshift 0, with the line of sight along the z axis. The
example data give each halo in real space too, so the corrected positions
can be compared with the true ones, as survey data do not allow.

Usage: python periodic_box.py [data_dir]. The catalogue is read from
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

    # ---------------------------------------------------
    # ------------ Read the tracer catalogue ------------
    # ---------------------------------------------------

    # Halos of a simulation at redshift 0 in a cube of side 1000 Mpc/h,
    # [0, 1000] on each axis, comoving coordinates in Mpc/h. The file holds
    # x, y and z in real space, then z in redshift space, with the line of
    # sight along the z axis: the positions observed are x, y and the fourth
    # column, those in real space the first three. Table.values has shape
    # (N, 3)
    tracers = otswap.io.read(data / "periodic_halos.dat", [0, 1, 3]).values
    real = otswap.io.read(data / "periodic_halos.dat", [0, 1, 2]).values
    side = 1000.0

    print(f"Tracers: {len(tracers)}")

    # ------------------------------------------------------------------
    # ------------ Set the parameters of the reconstruction ------------
    # ------------------------------------------------------------------

    # n_realizations: independent realizations, each on its own set of
    # randoms. With no randoms given, as here, otswap draws them itself,
    # uniformly in the box, a different set for each realization; 8
    # realizations keep the example fast.
    #
    # convergence: the sweeps stop once the fraction of successful swaps in
    # a sweep no longer exceeds this threshold; 1e-2 keeps the example fast.
    # For an analysis, start from the defaults, 30 realizations and 1e-3:
    # the guide's section "Choosing the parameters" advises on both.
    #
    # seed: a fixed seed makes the run reproducible, whatever the number of
    # threads.
    #
    # cell_size: the grid cell, in units of mps, affects the speed only,
    # never the result; the default is 3. In a periodic box the cells divide
    # the shortest side a whole number of times, so the cell used can be
    # smaller: result.cell_side
    config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, cell_size=3.0,
                           verbosity="normal")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # The box is periodic: the limits, the keyword box, are required, and
    # the mean particle separation, (V/N)^(1/3) with V the volume of the
    # limits, comes from them. The limits are closed: a tracer on an upper
    # face stays there, and one outside the limits would be wrapped into the
    # box and counted in result.selection; here none is. With the randoms
    # omitted, otswap draws them
    result = otswap.reconstruct_box(tracers, box=((0.0, 0.0, 0.0), (side, side, side)),
                                    periodic=True, config=config)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    print(f"Periodic: {'yes' if result.periodic else 'no'}, tracers wrapped: {result.selection.tracers_wrapped}")
    print(f"Mean particle separation: {result.mps:g} Mpc/h, grid cell side: {result.cell_side:g} Mpc/h")

    # Each displacement component lies within half a period, [-500, 500)
    # Mpc/h here, however far apart the tracer and its random lie in the
    # plain coordinates
    mean = result.mean_displacement
    print(f"Mean length of the mean displacement: {np.linalg.norm(mean, axis=1).mean():g} Mpc/h")
    print(f"Largest component of a mean displacement: {np.abs(mean).max():g} Mpc/h")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One row per tracer, in the order of the input: its position and its
    # mean Lagrangian position, both wrapped into the box, the mean
    # displacement, the number of valid realizations and the flag
    # outsideBox, always 0 in a periodic box, with a header describing the
    # columns (lengths in coordinate units), the box, the periodicity and
    # the parameters of the run. The writer is the C++ library's, so the
    # file is the C++ example's, byte for byte
    otswap.io.write_displacement_mean(output / "displacement_periodic_box.dat", result)

    print("Written: output/displacement_periodic_box.dat")

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # Only the linear growth rate f at the box's redshift is needed. The
    # table computes it from Omega_m, h, w0 and wa of a flat cosmology; here
    # Omega_m comes from the physical densities of baryons and cold dark
    # matter, omega_b = 0.02244 and omega_c = 0.1206, and h = 0.67
    h = 0.67
    omega_m = (0.02244 + 0.1206) / (h * h)
    distances = otswap.DistanceTable.flat(omega_m, h, w0=-1.0, wa=0.0)

    # ----------------------------------------------------------------
    # ------------ Correct the redshift-space distortions ------------
    # ----------------------------------------------------------------

    # The box's redshift, the tracers' linear bias, and the line of sight:
    # axis 0, 1 or 2 for x, y or z
    redshift, bias, axis = 0.0, 1.14, 2

    # sigma, the width of the Gaussian average in Mpc/h: 10 Mpc/h is a
    # starting value, and the value must be calibrated on the sample (the
    # guide's section "Choosing the parameters"). The comparison below shows
    # how close the corrected positions come to real space.
    #
    # The result is periodic, so the correction is too: the neighbours of a
    # tracer near a face include those beyond it, and each corrected
    # coordinate is wrapped into the box. catalogue.status records what was
    # done with each tracer
    correction = otswap.CorrectionConfig(weight_by_realizations=False)
    catalogue = otswap.real_space_box(result, axis=axis, redshift=redshift, distances=distances,
                                      bias=bias, sigma=10.0, correction=correction)

    # The single factor f/(b + 3f/5) of the box
    print(f"RSD factor f/(b + 3f/5): {catalogue.factor[0]:g}")

    # One row per tracer, in the order of the input: the corrected position,
    # NaN for an uncorrected tracer, then the number of valid realizations,
    # the two diagnostics of the average and the status of the tracer. The
    # writer is the C++ library's, so the file is the C++ example's, byte
    # for byte
    otswap.io.write_real_space_catalog(output / "reconstructed_catalogue_periodic_box.dat", catalogue)

    print("Written: output/reconstructed_catalogue_periodic_box.dat")

    # -----------------------------------------------------------
    # ------------ Compare with the real-space halos ------------
    # -----------------------------------------------------------

    # Along the line of sight, the offset of each halo from its real-space
    # position, before and after the correction. In a periodic box a halo
    # moved across a face reappears at the opposite one, so the offset is
    # the minimum image: d - L floor(d/L + 1/2), within half a period
    def minimum_image(d):
        return d - side * np.floor(d / side + 0.5)

    moved = catalogue.cartesian[:, axis]
    kept = ~np.isnan(moved)
    before = minimum_image(tracers[kept, axis] - real[kept, axis])
    after = minimum_image(moved[kept] - real[kept, axis])
    print(f"rms offset from the real-space position along the line of sight, over {np.count_nonzero(kept)} "
          f"halos: {np.sqrt(np.mean(before ** 2)):g} Mpc/h before the correction, "
          f"{np.sqrt(np.mean(after ** 2)):g} Mpc/h after")


if __name__ == "__main__":
    main()
