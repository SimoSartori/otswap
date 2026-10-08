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
 *  @file tests/test_writers.cpp
 *
 *  @brief The writers of a Result, a RealSpaceCatalog and an MpsProfile:
 *  their columns and keywords, and the lossless round trip of the complete
 *  displacement field.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <fitsio.h>

#include "otswap/OT.h"
#include "internal.h"
#include "check.h"

using namespace otswap;

namespace {

  const double kPi = 3.14159265358979323846;

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "otswap_test_writers_" + name;
  }

  std::string text_of (const std::string& file)
  {
    std::ifstream in(file);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  }

  // The names of the ### line of an ASCII file, joined by single spaces.
  std::string names_of (const std::string& file)
  {
    std::ifstream in(file);
    std::string line, names;
    while (std::getline(in, line))
      if (line.compare(0, 3, "###") == 0) {
        for (std::size_t i = 3; i < line.size(); ++i)
          if (line[i] != ' ' || (!names.empty() && names.back() != ' ')) names += line[i];
        break;
      }
    while (!names.empty() && names.back() == ' ') names.pop_back();
    return names;
  }

  // The value of a FITS keyword of HDU 2, as written.
  std::string fits_keyword (const std::string& file, const std::string& name)
  {
    fitsfile* f = nullptr;
    int status = 0, hdutype = 0;
    char value[FLEN_VALUE] = "", comment[FLEN_COMMENT] = "";
    fits_open_file(&f, file.c_str(), READONLY, &status);
    fits_movabs_hdu(f, 2, &hdutype, &status);
    fits_read_keyword(f, name.c_str(), value, comment, &status);
    int closing = 0;
    fits_close_file(f, &closing);
    return status == 0 ? std::string(value) : std::string("missing");
  }

  bool same_bits (const std::vector<double>& a, const std::vector<double>& b)
  {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
      if (!(std::isnan(a[i]) && std::isnan(b[i])) && std::memcmp(&a[i], &b[i], sizeof(double)) != 0)
        return false;
    return true;
  }

  std::vector<double> sky_points (const std::size_t n, std::mt19937& rng)
  {
    std::vector<double> sky(3 * n);
    for (std::size_t i = 0; i < n; ++i) {
      sky[3*i]   = internal::uniform_real(rng, 0.2, 0.7);
      sky[3*i+1] = internal::uniform_real(rng, -0.25, 0.25);
      sky[3*i+2] = internal::uniform_real(rng, 0.2, 0.7);
    }
    return sky;
  }

  template <typename F>
  bool throws_naming (F&& f, const std::string& text)
  {
    try {
      f();
    }
    catch (const Error& e) {
      return std::string(e.what()).find(text) != std::string::npos;
    }
    return false;
  }

  // Column c of a table read as double, as a vector.
  std::vector<double> column_of (const io::Table& t, const std::size_t c)
  {
    std::vector<double> out(t.nRows);
    for (std::size_t r = 0; r < t.nRows; ++r) out[r] = t.values[r*t.nColumns + c];
    return out;
  }

}

int main ()
{
  std::mt19937 rng(20261008);
  const std::size_t n = 900;
  const DistanceTable table(0.3, 0.7, -1., 0., 0., 1.5, 4000);
  const std::vector<double> sky = sky_points(n, rng);
  const std::vector<double> randomsSky = sky_points(3600, rng);
  std::vector<double> pixels(12 * 64 * 64, 1.);
  for (std::size_t p = 0; p < pixels.size(); p += 9) pixels[p] = 0.;
  const Mask mask(pixels, PixelOrdering::Ring);
  Config config;
  config.nRealizations = 3;
  config.seed = 4321;
  config.verbosity = Verbosity::Silent;
  const Result lc = reconstructLightcone(sky, randomsSky, mask, 1, table, config, RedshiftCut{0.3, 0.6});
  const Result box = reconstructBox(toCartesian(sky, table), {}, config);
  CorrectionConfig quiet;
  quiet.verbosity = Verbosity::Silent;
  const BiasTable bias {{0.3, 0.6}, {1.2, 1.6}};
  const RealSpaceCatalog lcCatalog = realSpaceLightcone(lc, table, bias, 10., quiet);
  const RealSpaceCatalog boxCatalog = realSpaceBox(box, 2, 0.5, table, 1.5, 8., quiet);

  group("writeDisplacements: the default columns of each geometry, CosmoBolognaLib's first");
  {
    const std::string file = temporary("displacements.dat");
    io::writeDisplacements(file, lc);
    check(names_of(file) == "tracRA tracDec tracRed lagrRA lagrDec lagrRed tracX tracY tracZ lagrX lagrY "
                            "lagrZ displX displY displZ nValidRec outsideRedshiftCut outsideMask",
          "the lightcone: " + names_of(file));
    const std::string text = text_of(file);
    check(text.find("## PRODUCT = displacements / what the file holds\n## OTSWAPV = 0.1.0 / otswap version\n"
                    "## GEOMETRY = lightcone / geometry of the reconstruction\n"
                    "## NOBJECTS = 900 / rows: one per input object\n## NREC = 3 / realizations\n"
                    "## SEED = 4321 / seed of the random streams\n") == 0,
          "the first keywords: " + text.substr(0, 400));
    check(text.find("## ZCUTMIN = 0.3 /") != std::string::npos && text.find("## ZCUTMAX = 0.6 /") != std::string::npos &&
          text.find("## MASK = 1 /") != std::string::npos && text.find("## FILTNSID = 64 /") != std::string::npos &&
          text.find("## MAXPIXCR = 0 /") != std::string::npos && text.find("## MPS") == std::string::npos,
          "the selection keywords, and no MPS in a lightcone");
    check(text.find("# tracRA               observed right ascension, declination and redshift of the tracer [deg]\n") !=
            std::string::npos, "descriptions with units");

    const io::Table t = io::read(file, {"tracRA", "tracDec", "tracRed", "lagrRA", "lagrX", "displX"},
                                 {"nValidRec", "outsideRedshiftCut", "outsideMask"});
    bool rows = t.nRows == n, angles = true, flags = true, nan = true;
    for (std::size_t i = 0; rows && i < n; ++i) {
      const double* v = &t.values[6*i];
      const std::int64_t* k = &t.integers[3*i];
      angles = angles && std::fabs(v[0] - sky[3*i] * 180. / kPi) < 1.e-6 && v[0] >= 0. && v[0] < 360. &&
               std::fabs(v[2] - sky[3*i+2]) < 1.e-8;
      flags = flags && k[0] == lc.validRealizations[i] && k[1] == lc.outsideRedshiftCut[i] &&
              k[2] == lc.outsideMask[i];
      if (lc.validRealizations[i] == 0) nan = nan && std::isnan(v[3]) && std::isnan(v[4]) && std::isnan(v[5]);
      else nan = nan && std::fabs(v[4] - lc.lagrangian[3*i]) < 1.e-5 * std::fabs(lc.lagrangian[3*i]);
    }
    check(rows && angles && flags && nan, "one row per tracer, in degrees, flags, NaN where nothing was computed");

    const std::string boxFile = temporary("box.dat");
    io::writeDisplacements(boxFile, box);
    check(names_of(boxFile) == "tracX tracY tracZ lagrX lagrY lagrZ displX displY displZ nValidRec",
          "the box: " + names_of(boxFile));
    check(text_of(boxFile).find("## MPS = ") != std::string::npos &&
          text_of(boxFile).find("## ZCUTMIN") == std::string::npos, "MPS in a box, no cut");

    io::writeDisplacements(boxFile, box, {io::DisplacementGroup::Index, io::DisplacementGroup::Displacement});
    const io::Table chosen = io::read(boxFile, {"displZ"}, {"index"});
    check(names_of(boxFile) == "index displX displY displZ" && chosen.integers[7] == 7,
          "groups in the order given; index the row");
    check(throws_naming([&] { io::writeDisplacements(boxFile, box, {io::DisplacementGroup::TracerSky}); },
                        "TracerSky describes a lightcone"), "a sky group on a box is refused");
    check(throws_naming([&] { io::writeDisplacements(boxFile, box, {io::DisplacementGroup::Tracer,
                                                                    io::DisplacementGroup::Tracer}); },
                        "given twice"), "a repeated group is refused");
    Result noLagrangian = box;
    noLagrangian.lagrangian.clear();
    check(throws_naming([&] { io::writeDisplacements(boxFile, noLagrangian); }, "cannot write lagrangian"),
          "a field of the wrong size is refused, before anything is written");
    std::remove(file.c_str());
    std::remove(boxFile.c_str());
  }

  group("writeDisplacements: FITS carries the same values, the keywords and the units");
  {
    const std::string file = temporary("displacements.fits");
    io::writeDisplacements(file, lc);
    check(fits_keyword(file, "NOBJECTS") == "900" && fits_keyword(file, "ZCUTMIN") == "0.3" &&
          fits_keyword(file, "PRODUCT").find("displacements") != std::string::npos &&
          fits_keyword(file, "TUNIT1").find("deg") != std::string::npos,
          "keywords and units in HDU 2");
    const io::Table t = io::read(file, {"lagrX", "displY"});
    bool exact = true;
    for (std::size_t i = 0; i < n; ++i)
      exact = exact && same_bits({t.values[2*i], t.values[2*i+1]}, {lc.lagrangian[3*i], lc.meanDisplacement[3*i+1]});
    check(exact, "the doubles themselves, NaN rows included");
    std::remove(file.c_str());
  }

  group("writeDisplacementField: every realization, losslessly, in ASCII and in FITS");
  {
    for (const std::string name : {"field.dat", "field.fits"}) {
      const std::string file = temporary(name);
      io::writeDisplacementField(file, lc);
      if (name == "field.dat")
        check(names_of(file) == "realization index tracX tracY tracZ lagrX lagrY lagrZ displX displY displZ "
                                "valid outsideRedshiftCut outsideMask", "the columns: " + names_of(file));
      const io::Table t = io::read(file, {"tracX", "tracY", "tracZ", "lagrX", "lagrY", "lagrZ",
                                          "displX", "displY", "displZ"},
                                   {"realization", "index", "valid", "outsideRedshiftCut", "outsideMask"});
      check(t.nRows == 3 * n, name + ": one row per realization and tracer");
      Result back;
      back.nObjects = n;
      back.nRealizations = 3;
      back.geometry = Geometry::Lightcone;
      back.distances = std::make_shared<const DistanceTable>(table);
      back.tracers.resize(3 * n);
      back.matchedRandom.resize(9 * n);
      back.displacement.resize(9 * n);
      back.valid.resize(3 * n);
      back.outsideRedshiftCut.resize(n);
      back.outsideMask.resize(n);
      bool order = true, difference = true;
      for (std::size_t row = 0; row < 3 * n; ++row) {
        const double* v = &t.values[9*row];
        const std::int64_t* k = &t.integers[5*row];
        order = order && k[0] == (std::int64_t)(row / n) && k[1] == (std::int64_t)(row % n);
        for (std::size_t c = 0; c < 3; ++c) {
          back.tracers[3*(row % n)+c] = v[c];
          back.matchedRandom[3*row+c] = v[3+c];
          back.displacement[3*row+c] = v[6+c];
          if (!std::isnan(v[6+c])) difference = difference && v[6+c] == v[3+c] - v[c];
        }
        back.valid[row] = (std::uint8_t)k[2];
        back.outsideRedshiftCut[row % n] = (std::uint8_t)k[3];
        back.outsideMask[row % n] = (std::uint8_t)k[4];
      }
      recomputeMeans(back);
      check(order, name + ": ordered [realization][object]");
      check(difference, name + ": displacement = lagr - trac, bit for bit");
      check(same_bits(back.displacement, lc.displacement) && same_bits(back.matchedRandom, lc.matchedRandom) &&
            same_bits(back.tracers, lc.tracers) && back.valid == lc.valid &&
            back.outsideMask == lc.outsideMask && back.outsideRedshiftCut == lc.outsideRedshiftCut,
            name + ": every array read back identical");
      check(same_bits(back.meanDisplacement, lc.meanDisplacement) &&
            back.validRealizations == lc.validRealizations && same_bits(back.lagrangian, lc.lagrangian) &&
            same_bits(back.lagrangianSky, lc.lagrangianSky),
            name + ": recomputeMeans gives the reconstruction's means, lagrangian and lagrangianSky");
      std::remove(file.c_str());
    }
    Result malformed = lc;
    malformed.valid.pop_back();
    check(throws_naming([&] { io::writeDisplacementField(temporary("x.dat"), malformed); }, "cannot write valid"),
          "a malformed result is refused");
  }

  group("writeRealSpaceCatalog: the default columns, the status, the keywords");
  {
    const std::string file = temporary("catalogue.dat");
    io::writeRealSpaceCatalog(file, lcCatalog);
    check(names_of(file) == "tracRA tracDec tracRed tracX tracY tracZ nValidRec nNeighbours "
                            "nRealizationsAveraged status", "the lightcone: " + names_of(file));
    const std::string text = text_of(file);
    check(text.find("## PRODUCT = real-space catalogue /") == 0 && text.find("## SIGMA = 10 /") != std::string::npos &&
          text.find("## WEIGHTED = 0 /") != std::string::npos && text.find("## NEXTRAP = ") != std::string::npos,
          "keywords: " + text.substr(0, 300));
    const io::Table t = io::read(file, {"tracRA", "tracRed", "tracX"}, {"status", "nNeighbours"});
    bool right = t.nRows == n;
    for (std::size_t i = 0; right && i < n; ++i)
      right = t.integers[2*i] == (std::int64_t)lcCatalog.status[i] &&
              t.integers[2*i+1] == lcCatalog.nNeighbours[i] &&
              (std::isnan(lcCatalog.sky[3*i]) ? std::isnan(t.values[3*i])
                                             : std::fabs(t.values[3*i] - lcCatalog.sky[3*i] * 180. / kPi) < 1.e-6);
    check(right, "status, neighbours and degrees per row");

    const std::string boxFile = temporary("catalogue_box.fits");
    io::writeRealSpaceCatalog(boxFile, boxCatalog);
    check(fits_keyword(boxFile, "AXIS") == "2" && fits_keyword(boxFile, "ZBOX") == "0.5" &&
          fits_keyword(boxFile, "BIAS") == "1.5" && fits_keyword(boxFile, "SIGMA") == "8",
          "box keywords in FITS");
    const io::Table b = io::read(boxFile, {"tracZ"}, {"nValidRec"});
    check(b.nRows == n && same_bits(column_of(b, 0), [&] {
            std::vector<double> z(n);
            for (std::size_t i = 0; i < n; ++i) z[i] = boxCatalog.cartesian[3*i+2];
            return z; }()), "the corrected box positions, exactly");

    io::writeRealSpaceCatalog(boxFile, boxCatalog, {io::CatalogGroup::Index, io::CatalogGroup::Shift});
    const io::Table s = io::read(boxFile, {"shift", "rsdFactor"}, {"index"});
    check(s.values[2*5] == boxCatalog.shift[5] && s.values[2*5+1] == boxCatalog.factor[5] && s.integers[5] == 5,
          "the Shift group");
    check(throws_naming([&] { io::writeRealSpaceCatalog(boxFile, boxCatalog, {io::CatalogGroup::Sky}); },
                        "Sky describes a lightcone"), "Sky on a box is refused");
    std::remove(file.c_str());
    std::remove(boxFile.c_str());
  }

  group("writeMpsProfile: one row per bin, and the profile's keywords");
  {
    const std::string file = temporary("mps.dat");
    io::writeMpsProfile(file, lc.mpsProfile);
    check(names_of(file) == "redshift MPS nTracers", "the columns");
    const io::Table t = io::read(file, {"redshift", "MPS"}, {"nTracers"});
    check(t.nRows == 1 && std::fabs(t.values[1] - lc.mpsProfile.mps[0]) < 1.e-6 * lc.mpsProfile.mps[0] &&
          t.integers[0] == (std::int64_t)lc.mpsProfile.count[0], "the bin");
    const std::string text = text_of(file);
    check(text.find("## PRODUCT = mps profile /") == 0 && text.find("## MPSREPR = ") != std::string::npos &&
          text.find("## ZMIN = ") != std::string::npos, "keywords: " + text.substr(0, 300));
    check(throws_naming([&] { io::writeMpsProfile(file, MpsProfile()); }, "empty"), "an empty profile is refused");
    std::remove(file.c_str());
  }

  return report("test_writers");
}
