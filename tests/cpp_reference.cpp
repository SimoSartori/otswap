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
 *  @file tests/cpp_reference.cpp
 *
 *  @brief Runs the C++ library on inputs written to raw files, and writes
 *  its outputs to raw files, so that the Python tests can check that the
 *  bindings return exactly what the library returns.
 *
 *  Not a test by itself. Usage, with key=value arguments:
 *
 *    cpp_reference box tracers=F [randoms=F] mps=X n_realizations=N
 *                      convergence=X seed=N cell_size=X out=P
 *                      [mask=FITS max_forbidden_pixels=N]
 *
 *    cpp_reference lightcone tracers_sky=F randoms_sky=F [tracers=F randoms=F]
 *                      sky_area_deg2=X n_bins=N omega_m=X h=X z_min=X z_max=X
 *                      n_samples=N n_realizations=N convergence=X seed=N
 *                      cell_size=X out=P [mask=FITS max_forbidden_pixels=N]
 *
 *    cpp_reference cartesian sky=F omega_m=X h=X z_min=X z_max=X n_samples=N out=P
 *
 *  Input files hold native float64 values, flat and row-major; sky angles
 *  are in radians. A reconstruction writes P.displacement, P.matched_random
 *  and P.mean_displacement (float64), P.valid (uint8) and
 *  P.valid_realizations (uint32); cartesian writes P.xyz (float64). The
 *  table of a lightcone is DistanceTable(omega_m, h, -1, 0, z_min, z_max,
 *  n_samples).
 */

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "otswap/OT.h"

namespace {

  std::map<std::string, std::string> g_args;

  bool has (const std::string& key) { return g_args.count(key) != 0; }

  const std::string& arg (const std::string& key)
  {
    const auto it = g_args.find(key);
    if (it == g_args.end()) throw otswap::Error("missing argument " + key);
    return it->second;
  }

  double real (const std::string& key) { return std::strtod(arg(key).c_str(), nullptr); }

  unsigned integer (const std::string& key)
  {
    return (unsigned)std::strtoul(arg(key).c_str(), nullptr, 10);
  }

  std::vector<double> read (const std::string& key)
  {
    std::ifstream in(arg(key), std::ios::binary | std::ios::ate);
    if (!in) throw otswap::Error("cannot open " + arg(key));
    std::vector<double> values((std::size_t)in.tellg() / sizeof(double));
    in.seekg(0);
    in.read(reinterpret_cast<char*>(values.data()), (std::streamsize)(values.size() * sizeof(double)));
    return values;
  }

  template <typename T>
  void write (const std::string& suffix, const std::vector<T>& values)
  {
    std::ofstream out(arg("out") + "." + suffix, std::ios::binary);
    out.write(reinterpret_cast<const char*>(values.data()), (std::streamsize)(values.size() * sizeof(T)));
    if (!out) throw otswap::Error("cannot write " + arg("out") + "." + suffix);
  }

  otswap::Config config ()
  {
    otswap::Config c;
    c.nRealizations = integer("n_realizations");
    c.convergence = real("convergence");
    c.seed = integer("seed");
    c.cellSize = real("cell_size");
    return c;
  }

  otswap::DistanceTable table ()
  {
    return otswap::DistanceTable(real("omega_m"), real("h"), -1., 0., real("z_min"), real("z_max"),
                                 integer("n_samples"));
  }

  void finish (otswap::Result& result)
  {
    if (has("mask"))
      otswap::rejectMaskCrossings(result, otswap::Mask(arg("mask")), integer("max_forbidden_pixels"));
    write("displacement", result.displacement);
    write("matched_random", result.matchedRandom);
    write("mean_displacement", result.meanDisplacement);
    write("valid", result.valid);
    write("valid_realizations", result.validRealizations);
  }

}

int main (int argc, char** argv)
{
  if (argc < 2) {
    std::fprintf(stderr, "usage: cpp_reference box|lightcone|cartesian key=value...\n");
    return 2;
  }
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    const std::size_t eq = a.find('=');
    if (eq == std::string::npos) {
      std::fprintf(stderr, "cpp_reference: argument %s is not key=value\n", argv[i]);
      return 2;
    }
    g_args[a.substr(0, eq)] = a.substr(eq + 1);
  }

  try {
    const std::string mode = argv[1];
    if (mode == "box") {
      const std::vector<double> randoms = has("randoms") ? read("randoms") : std::vector<double>();
      otswap::Result result = otswap::reconstructBox(read("tracers"), randoms, real("mps"), config());
      finish(result);
    }
    else if (mode == "lightcone") {
      const otswap::DistanceTable distances = table();
      otswap::Result result = has("tracers")
        ? otswap::reconstructLightcone(read("tracers"), read("randoms"), read("tracers_sky"),
                                       read("randoms_sky"), real("sky_area_deg2"), integer("n_bins"),
                                       distances, config())
        : otswap::reconstructLightcone(read("tracers_sky"), read("randoms_sky"), real("sky_area_deg2"),
                                       integer("n_bins"), distances, config());
      finish(result);
    }
    else if (mode == "cartesian") {
      write("xyz", otswap::toCartesian(read("sky"), table()));
    }
    else {
      std::fprintf(stderr, "cpp_reference: unknown mode %s\n", mode.c_str());
      return 2;
    }
  }
  catch (const std::exception& e) {
    std::fprintf(stderr, "cpp_reference: %s\n", e.what());
    return 1;
  }
  return 0;
}
