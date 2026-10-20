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
 *  @file tests/data/write_fixtures_v010.cpp
 *
 *  @brief Writes the otswap 0.1.0 field files of tests/data, which the
 *  reader of complete displacement files must still read: a box, and a
 *  lightcone with a mask and a redshift cut, each in ASCII and in FITS.
 *
 *  Built against an installed otswap 0.1.0 (tag v0.1.0), with
 *  -ffp-contract=off, not by the build of this repository:
 *
 *      c++ -std=c++17 -O2 -ffp-contract=off -I<otswap 0.1.0>/include -Itests \
 *          tests/data/write_fixtures_v010.cpp <libotswap.a, cfitsio, OpenMP> -o write_fixtures
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
  const otswap::Result lightcone = otswap::reconstructLightcone(lightcone_tracers(), lightcone_randoms(),
                                                                mask(), 1, table(), config(), cut());
  for (const std::string ext : {".dat", ".fits"}) {
    otswap::io::writeDisplacementField(dir + "v0.1.0_box" + ext, box);
    otswap::io::writeDisplacementField(dir + "v0.1.0_lightcone" + ext, lightcone);
  }
  return 0;
}
