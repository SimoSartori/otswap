# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Optimal transport reconstruction of the displacement field of a tracer
catalogue, by local swapping.

Coordinates are NumPy arrays of shape (N, 3). Cartesian coordinates are in
Mpc/h; sky coordinates are ordered right ascension, declination, redshift,
with the angular unit given explicitly by ``angle_unit``. Inputs are
converted to float64; the reconstruction works on its own copy of them.

The number of threads follows ``OMP_NUM_THREADS``. The reconstruction
functions release the GIL while they run.

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

from enum import IntEnum
from typing import Literal, NamedTuple

import numpy as np
from numpy.typing import NDArray

from ._otswap import (
    BiasTable,
    DistanceTable,
    Error,
    ExtrapolationWarning,
    Mask,
    MpsProfile,
    RealSpaceCatalog,
    Result,
    SelectionCounts,
    axis_projection,
    neighbour_average,
    real_space_box,
    real_space_lightcone,
    recompute_means,
    reconstruct_box,
    reconstruct_lightcone,
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


__all__ = [
    "AngleUnit",
    "BiasTable",
    "CorrectionStatus",
    "DistanceTable",
    "Error",
    "ExtrapolationWarning",
    "Mask",
    "MpsProfile",
    "RealSpaceCatalog",
    "Result",
    "NeighbourAverage",
    "SelectionCounts",
    "axis_projection",
    "neighbour_average",
    "radial_projection",
    "real_space_box",
    "real_space_lightcone",
    "recompute_means",
    "reconstruct_box",
    "reconstruct_lightcone",
    "reject_mask_crossings",
    "rsd_factor",
    "rsd_factor_box",
    "shift_along_axis",
    "shift_radially",
    "to_cartesian",
    "to_sky",
    "Verbosity",
]
