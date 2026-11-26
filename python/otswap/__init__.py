# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Optimal transport reconstruction of the displacement field of a tracer
catalogue, by local swapping.

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates, and
every length, are in coordinate units: the unit of the tracers given in a box,
that of the distance table in a lightcone (Mpc/h from
``DistanceTable.flat``); sky coordinates are ordered right ascension, declination, redshift,
with the angular unit given explicitly by ``angle_unit``. Inputs are
converted to float64; the reconstruction works on its own copy of them.

The number of threads follows ``OMP_NUM_THREADS``. The reconstruction
functions release the GIL while they run.

The reconstructions, ``measure_mps_profile`` and the corrections print their
report to ``sys.stdout`` as they run, each line as its step finishes. Ctrl+C
stops such a call at its next report line with ``KeyboardInterrupt``; a
silent call, or a long step, runs on until its next line or its end, and a
call made from a thread other than the main one is not interrupted. An
exception raised by ``sys.stdout`` stops the printing, and the call raises it
when it ends.

Tables are read and written by the submodule ``otswap.io``.

The redshift-space correction of tracer i is a shift along its line of
sight::

    s_i = f(z_i) / (b(z_i) + 3 f(z_i) / 5) * <Psi . r_hat>_i

with Psi the reconstructed displacement (``Result.mean_displacement``, from
the observed to the reconstructed position), r_hat the line of sight, f the
linear growth rate, b the linear bias and < > a gaussian average over the
neighbouring tracers. ``real_space_lightcone`` and ``real_space_box`` run the
whole chain; the four steps are also available on their own.
"""

import os
from enum import IntEnum
from typing import Literal, NamedTuple

import numpy as np
from numpy.typing import NDArray

from ._otswap import (
    BiasTable,
    BoxLimits,
    Config,
    CorrectionConfig,
    DistanceCut,
    DistanceTable,
    Error,
    ExtrapolationWarning,
    Mask,
    MpsProfile,
    RealSpaceCatalog,
    RedshiftCut,
    Result,
    SelectionCounts,
    axis_projection,
    measure_mps_profile,
    measure_mps_profile_cartesian,
    neighbour_average,
    point_to_sky,
    real_space_box,
    real_space_lightcone,
    recompute_means,
    reconstruct_box,
    reconstruct_lightcone,
    reconstruct_lightcone_cartesian,
    radial_projection,
    reject_mask_crossings,
    rsd_factor,
    rsd_factor_box,
    shift_along_axis,
    shift_radially,
    to_cartesian,
    to_sky,
)

from . import io

AngleUnit = Literal["deg", "rad"]
Verbosity = Literal["silent", "normal", "detailed"]


class CorrectionStatus(IntEnum):
    """The values of ``RealSpaceCatalog.status``."""

    CORRECTED = 0
    MOVED_BY_NEIGHBOURS = 1
    NO_VALID_NEIGHBOUR = 2
    LEFT_OUT = 3


class NeighbourAverage(NamedTuple):
    """What ``neighbour_average`` returns, shape (N,) each."""

    values: NDArray[np.float64]
    n_neighbours: NDArray[np.uint32]
    n_realizations_averaged: NDArray[np.uint32]


def openmp_root():
    """The directory of the OpenMP runtime a macOS wheel of otswap ships,
    with ``lib/libomp.dylib`` and ``include/omp.h``; None when otswap ships
    no runtime."""
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.join(here, "openmp")
    if (os.path.isfile(os.path.join(root, "lib", "libomp.dylib"))
            and os.path.isfile(os.path.join(root, "include", "omp.h"))):
        return root
    return None


__all__ = [
    "AngleUnit",
    "BiasTable",
    "BoxLimits",
    "Config",
    "CorrectionConfig",
    "CorrectionStatus",
    "DistanceCut",
    "DistanceTable",
    "Error",
    "ExtrapolationWarning",
    "Mask",
    "MpsProfile",
    "RealSpaceCatalog",
    "RedshiftCut",
    "Result",
    "NeighbourAverage",
    "SelectionCounts",
    "axis_projection",
    "measure_mps_profile",
    "measure_mps_profile_cartesian",
    "neighbour_average",
    "openmp_root",
    "point_to_sky",
    "radial_projection",
    "real_space_box",
    "real_space_lightcone",
    "recompute_means",
    "reconstruct_box",
    "reconstruct_lightcone",
    "reconstruct_lightcone_cartesian",
    "reject_mask_crossings",
    "rsd_factor",
    "rsd_factor_box",
    "shift_along_axis",
    "shift_radially",
    "to_cartesian",
    "to_sky",
    "Verbosity",
]
