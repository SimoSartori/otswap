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
 *  @file tests/data/write_fixtures_v021.cpp
 *
 *  @brief Writes the otswap 0.2.1 mean displacement files of tests/data,
 *  which the reader of mean displacement files must still read: a box, and
 *  a lightcone with a mask, a redshift cut and IDs, each in ASCII and in
 *  FITS, with the default columns and no profile or distance table.
 *
 *  Built against an installed otswap 0.2.1 (tag v0.2.1), with
 *  -ffp-contract=off, not by the build of this repository:
 *
 *      c++ -std=c++17 -O2 -ffp-contract=off -I<otswap 0.2.1>/include -Itests \
 *          tests/data/write_fixtures_v021.cpp <libotswap.a, libmeshsearch.a, cfitsio, OpenMP> \
 *          -o write_fixtures
 *      ./write_fixtures tests/data
 */

#include <iostream>
#include <string>

#include "otswap/OT.h"
#include "fixtures_v010.h"

int main (int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << "usage: write_fixtures <directory>" << std::endl;
    return 1;
  }
  const std::string dir = std::string(argv[1]) + "/";

  using namespace fixtures_v010;
  const otswap::Result box = otswap::reconstructBox(box_tracers(), {}, config());
  otswap::Result lightcone = otswap::reconstructLightcone(lightcone_tracers(), lightcone_randoms(),
                                                          mask(), 1, table(), config(), cut());
  lightcone.ids = lightcone_ids();
  for (const std::string ext : {".dat", ".fits"}) {
    otswap::io::writeDisplacementMean(dir + "v0.2.1_mean_box" + ext, box);
    otswap::io::writeDisplacementMean(dir + "v0.2.1_mean_lightcone" + ext, lightcone);
  }
  return 0;
}
