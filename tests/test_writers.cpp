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
 *  their columns and keywords, the lossless round trip of the complete
 *  displacement file, its profile and distance table parts, and the
 *  refusals of the displacement writers.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
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

  // The value of a keyword read by internal::read_keywords, or "missing".
  std::string value_of (const std::vector<io::Keyword>& keywords, const std::string& name)
  {
    for (const io::Keyword& k : keywords)
      if (k.name == name) return k.value;
    return "missing";
  }

  // The EXTNAME of an HDU of a FITS file, or "missing".
  std::string extname_of (const std::string& file, const int hdu)
  {
    fitsfile* f = nullptr;
    int status = 0, hdutype = 0;
    char value[FLEN_VALUE] = "";
    fits_open_file(&f, file.c_str(), READONLY, &status);
    fits_movabs_hdu(f, hdu, &hdutype, &status);
    fits_read_key(f, TSTRING, "EXTNAME", value, nullptr, &status);
    int closing = 0;
    fits_close_file(f, &closing);
    return status == 0 ? std::string(value) : std::string("missing");
  }

  bool exists (const std::string& file)
  {
    return std::ifstream(file).good();
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
  BoxLimits half = *box.box;
  half.upper[0] = 0.5 * (half.lower[0] + half.upper[0]);
  const Result limited = reconstructBox(toCartesian(sky, table), {}, config, half);
  CorrectionConfig quiet;
  quiet.verbosity = Verbosity::Silent;
  const BiasTable bias {{0.3, 0.6}, {1.2, 1.6}};
  const RealSpaceCatalog lcCatalog = realSpaceLightcone(lc, table, bias, 10., quiet);
  const RealSpaceCatalog boxCatalog = realSpaceBox(box, 2, 0.5, table, 1.5, 8., quiet);

  group("writeDisplacementMean: the default columns of each geometry, CosmoBolognaLib's first");
  {
    const std::string file = temporary("displacements.dat");
    io::writeDisplacementMean(file, lc);
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
    io::writeDisplacementMean(boxFile, box);
    check(names_of(boxFile) == "tracX tracY tracZ lagrX lagrY lagrZ displX displY displZ nValidRec",
          "the box: " + names_of(boxFile));
    check(text_of(boxFile).find("## MPS = ") != std::string::npos &&
          text_of(boxFile).find("## ZCUTMIN") == std::string::npos, "MPS in a box, no cut");

    io::writeDisplacementMean(boxFile, box, {io::DisplacementGroup::Id, io::DisplacementGroup::Displacement});
    const io::Table chosen = io::read(boxFile, {"displZ"}, {"id"});
    check(names_of(boxFile) == "id displX displY displZ" && chosen.integers[7] == 7,
          "groups in the order given; without IDs, id the row");
    check(throws_naming([&] { io::writeDisplacementMean(boxFile, box, {io::DisplacementGroup::TracerSky}); },
                        "TracerSky describes a lightcone"), "a sky group on a box is refused");
    check(throws_naming([&] { io::writeDisplacementMean(boxFile, box, {io::DisplacementGroup::Tracer,
                                                                    io::DisplacementGroup::Tracer}); },
                        "given twice"), "a repeated group is refused");
    Result noLagrangian = box;
    noLagrangian.lagrangian.clear();
    check(throws_naming([&] { io::writeDisplacementMean(boxFile, noLagrangian); }, "cannot write lagrangian"),
          "a field of the wrong size is refused, before anything is written");

    io::writeDisplacementMean(boxFile, limited);
    const io::Table flagged = io::read(boxFile, {"tracX"}, {"outsideBox"});
    bool same = flagged.nRows == n;
    for (std::size_t i = 0; same && i < n; ++i) same = flagged.integers[i] == limited.outsideBox[i];
    check(names_of(boxFile) == "tracX tracY tracZ lagrX lagrY lagrZ displX displY displZ nValidRec outsideBox" &&
          same && std::count(limited.outsideBox.begin(), limited.outsideBox.end(), 1) > 0,
          "a box with limits: outsideBox appended to the default, as in the result");
    io::writeDisplacementMean(boxFile, box, {io::DisplacementGroup::Id, io::DisplacementGroup::OutsideBox});
    const io::Table zeros = io::read(boxFile, {"outsideBox"});
    bool allZero = zeros.nRows == n;
    for (std::size_t i = 0; allZero && i < n; ++i) allZero = zeros.values[i] == 0.;
    check(names_of(boxFile) == "id outsideBox" && allZero,
          "a box without limits: OutsideBox written on request, all 0");
    check(throws_naming([&] { io::writeDisplacementMean(file, lc, {io::DisplacementGroup::OutsideBox}); },
                        "OutsideBox describes a box"), "OutsideBox on a lightcone is refused");
    std::remove(file.c_str());
    std::remove(boxFile.c_str());
  }

  group("writeDisplacementMean: FITS carries the same values, the keywords and the units");
  {
    const std::string file = temporary("displacements.fits");
    io::writeDisplacementMean(file, lc);
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

  group("writeDisplacementComplete: every realization, losslessly, in ASCII and in FITS");
  {
    for (const std::string name : {"complete.dat", "complete.fits"}) {
      const std::string file = temporary(name);
      io::writeDisplacementComplete(file, lc);
      if (name == "complete.dat")
        check(names_of(file) == "realization id tracRA tracDec tracRed tracX tracY tracZ lagrX lagrY lagrZ "
                                "displX displY displZ valid outsideRedshiftCut outsideMask outsideBox",
              "the columns: " + names_of(file));
      const std::string product = name == "complete.dat" ? text_of(file).substr(0, 38)
                                                         : fits_keyword(file, "PRODUCT");
      check(product == (name == "complete.dat" ? "## PRODUCT = complete displacements / " : "'complete displacements'"),
            name + ": PRODUCT is complete displacements: " + product);
      const io::Table t = io::read(file, {"tracX", "tracY", "tracZ", "lagrX", "lagrY", "lagrZ",
                                          "displX", "displY", "displZ"},
                                   {"realization", "id", "valid", "outsideRedshiftCut", "outsideMask"});
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
      const io::Table extra = io::read(file, {"tracRA", "tracDec", "tracRed"}, {"outsideBox"});
      bool degrees = extra.nRows == 3 * n, inside = degrees;
      std::vector<double> skyBack(3 * n);
      for (std::size_t row = 0; degrees && row < 3 * n; ++row) {
        const std::size_t i = row % n;
        const double* v = &extra.values[3*row];
        degrees = v[0] == sky[3*i] * (180. / kPi) && v[1] == sky[3*i+1] * (180. / kPi) && v[2] == sky[3*i+2];
        inside = inside && extra.integers[row] == 0;
        if (row < n) { skyBack[3*i] = v[0]; skyBack[3*i+1] = v[1]; skyBack[3*i+2] = v[2]; }
      }
      skyBack = skyToRadians(skyBack);
      bool close = degrees;
      for (std::size_t j = 0; close && j < 3 * n; ++j)
        close = std::fabs(skyBack[j] - sky[j]) <= 4.e-16 * std::fabs(sky[j]);
      check(degrees, name + ": tracRA, tracDec in degrees, tracRed, in every realization");
      check(close, name + ": the sky positions back in radians, within the rounding of the conversion");
      check(inside, name + ": outsideBox 0 in a lightcone");
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
    check(throws_naming([&] { io::writeDisplacementComplete(temporary("x.dat"), malformed); }, "cannot write valid"),
          "a malformed result is refused");
  }

  group("Both displacement writers: the reconstruction's keywords, in their order, each read back as written");
  {
    const char* const order[] = {"PRODUCT", "OTSWAPV", "GEOMETRY", "NOBJECTS", "NREC", "SEED", "CONVERG",
                                 "CELLSIZE", "REJCROSS", "CFMAXPIX", "VERBOSE", "MPS", "BOXXLO", "BOXXHI",
                                 "BOXYLO", "BOXYHI", "BOXZLO", "BOXZHI", "BOXLIM", "ZCUTMIN", "ZCUTMAX",
                                 "MASK", "FILTNSID", "MAXPIXCR", "NTRAC", "NTRCUT", "NTRMASK", "NTRBOTH",
                                 "NTRBOX", "NRAND", "NRNCUT", "NRNMASK", "NRNBOTH", "NRNBOX", "NDISPL",
                                 "NDISPCR", "MPSPART", "DISTPART"};
    // Each count is written when it is not 0, as it is.
    auto counts_match = [] (const std::vector<io::Keyword>& k, const SelectionCounts& s) {
      const std::pair<const char*, std::size_t> counts[] = {
        {"NTRAC", s.tracers}, {"NTRCUT", s.tracersOutsideRedshiftCut}, {"NTRMASK", s.tracersOutsideMask},
        {"NTRBOTH", s.tracersOutsideBoth}, {"NTRBOX", s.tracersOutsideBox}, {"NRAND", s.randoms},
        {"NRNCUT", s.randomsOutsideRedshiftCut}, {"NRNMASK", s.randomsOutsideMask},
        {"NRNBOTH", s.randomsOutsideBoth}, {"NRNBOX", s.randomsOutsideBox}, {"NDISPL", s.displacements},
        {"NDISPCR", s.displacementsCrossingMask}};
      for (const auto& c : counts)
        if (value_of(k, c.first) != (c.second == 0 ? std::string("missing") : std::to_string(c.second)))
          return false;
      return true;
    };
    // The box keywords read back as the bounds, exactly.
    auto box_matches = [] (const std::vector<io::Keyword>& k, const std::optional<BoxLimits>& b) {
      const char* const names[3][2] = {{"BOXXLO", "BOXXHI"}, {"BOXYLO", "BOXYHI"}, {"BOXZLO", "BOXZHI"}};
      for (int c = 0; c < 3; ++c) {
        if (!b) {
          if (value_of(k, names[c][0]) != "missing" || value_of(k, names[c][1]) != "missing") return false;
          continue;
        }
        if (std::strtod(value_of(k, names[c][0]).c_str(), nullptr) != b->lower[c] ||
            std::strtod(value_of(k, names[c][1]).c_str(), nullptr) != b->upper[c]) return false;
      }
      return true;
    };
    Result hand = box;
    SelectionCounts& counts = hand.selection;
    counts.tracers = 1; counts.tracersOutsideRedshiftCut = 2; counts.tracersOutsideMask = 3;
    counts.tracersOutsideBoth = 4; counts.tracersOutsideBox = 5; counts.randoms = 6;
    counts.randomsOutsideRedshiftCut = 7; counts.randomsOutsideMask = 8; counts.randomsOutsideBoth = 9;
    counts.randomsOutsideBox = 10; counts.displacements = 11; counts.displacementsCrossingMask = 12;
    hand.config.cellSize = 2.5;
    hand.config.rejectCrossings = false;
    hand.config.maxUnobservedPixelsCrossed = 7;
    hand.config.verbosity = Verbosity::Detailed;
    hand.config.convergence = 0.1;

    for (const std::string extension : {".dat", ".fits"})
      for (const bool complete : {false, true}) {
        const std::string file = temporary("keywords" + extension);
        const std::string label = extension + (complete ? ", complete: " : ", mean: ");
        auto keywords_of = [&] (const Result& r) {
          if (complete) io::writeDisplacementComplete(file, r);
          else          io::writeDisplacementMean(file, r);
          return internal::read_keywords(file);
        };

        const std::vector<io::Keyword> k = keywords_of(lc);
        std::size_t next = 0;
        bool ordered = true;
        for (const io::Keyword& keyword : k) {
          std::size_t at = next;
          while (at < std::size(order) && keyword.name != order[at]) ++at;
          if (at == std::size(order)) {
            const bool structural = extension == ".fits" && (keyword.name == "XTENSION" ||
              keyword.name == "BITPIX" || keyword.name.compare(0, 5, "NAXIS") == 0 ||
              keyword.name == "PCOUNT" || keyword.name == "GCOUNT" || keyword.name == "TFIELDS" ||
              keyword.name.compare(0, 5, "TTYPE") == 0 || keyword.name.compare(0, 5, "TFORM") == 0 ||
              keyword.name.compare(0, 5, "TCOMM") == 0 || keyword.name.compare(0, 5, "TUNIT") == 0 ||
              keyword.name == "EXTNAME");
            ordered = ordered && structural;
            continue;
          }
          next = at + 1;
        }
        check(ordered, label + "the keywords in the order of the format");
        check(value_of(k, "CELLSIZE") == "4" && value_of(k, "REJCROSS") == "1" &&
              value_of(k, "CFMAXPIX") == "0" && value_of(k, "VERBOSE") == "silent" &&
              value_of(k, "NREC") == "3" && value_of(k, "SEED") == "4321",
              label + "the configuration of the lightcone");
        check(counts_match(k, lc.selection) && value_of(k, "NTRAC") == "900" &&
              value_of(k, "NRAND") == "3600" && lc.selection.displacementsCrossingMask > 0,
              label + "the lightcone's counts, those not 0");
        check(box_matches(k, lc.box) && value_of(k, "BOXLIM") == "missing", label + "no box in a lightcone");
        if (extension == ".fits")
          check(value_of(k, "EXTNAME") == "DISPLACEMENTS" && internal::has_part(file, "DISPLACEMENTS"),
                label + "the table named DISPLACEMENTS");
        if (complete)
          check(value_of(k, "MPSPART") == "1" && value_of(k, "DISTPART") == "1",
                label + "the parts declared");
        else
          check(value_of(k, "MPSPART") == "missing" && value_of(k, "DISTPART") == "missing",
                label + "no part in the mean file");

        const std::vector<io::Keyword> b = keywords_of(box);
        check(box_matches(b, box.box) && value_of(b, "BOXLIM") == "missing" && counts_match(b, box.selection) &&
              value_of(b, "NTRAC") == "missing" && value_of(b, "MPSPART") == "missing" &&
              value_of(b, "DISTPART") == "missing" && std::strtod(value_of(b, "MPS").c_str(), nullptr) == box.mps,
              label + "a box without limits: its bounding box, MPS, no count");
        const std::vector<io::Keyword> l = keywords_of(limited);
        check(box_matches(l, limited.box) && value_of(l, "BOXLIM") == "1" && counts_match(l, limited.selection) &&
              value_of(l, "NTRAC") == "900" && limited.selection.tracersOutsideBox > 0,
              label + "a box with limits: the limits, BOXLIM, the counts");
        const std::vector<io::Keyword> h = keywords_of(hand);
        check(counts_match(h, hand.selection) && value_of(h, "CELLSIZE") == "2.5" &&
              value_of(h, "REJCROSS") == "0" && value_of(h, "CFMAXPIX") == "7" &&
              value_of(h, "VERBOSE") == "detailed" && value_of(h, "CONVERG") == "0.1",
              label + "every count and the configuration of a hand-built result");
        std::remove(file.c_str());
      }
  }

  group("writeDisplacementComplete: the profile and the distance table, HDUs in FITS, twin files in ASCII");
  {
    const MpsProfile& p = lc.mpsProfile;
    std::vector<double> count(p.count.size());
    for (std::size_t i = 0; i < count.size(); ++i) count[i] = (double)p.count[i];

    const std::string fits = temporary("parts.fits");
    io::writeDisplacementComplete(fits, lc);
    check(extname_of(fits, 2) == "DISPLACEMENTS" && extname_of(fits, 3) == "MPSPROF" &&
          extname_of(fits, 4) == "DISTANCE" && extname_of(fits, 5) == "missing",
          "FITS: DISPLACEMENTS, MPSPROF, DISTANCE, in HDUs 2 to 4");
    const std::vector<io::Keyword> pk = internal::read_keywords(fits, "MPSPROF");
    const std::vector<io::Keyword> tk = internal::read_keywords(fits, "DISTANCE");
    check(value_of(pk, "PRODUCT") == "complete displacements mps profile" &&
          std::strtod(value_of(pk, "ZMIN").c_str(), nullptr) == p.redshiftMin &&
          std::strtod(value_of(pk, "ZMAX").c_str(), nullptr) == p.redshiftMax &&
          std::strtod(value_of(pk, "MPSREPR").c_str(), nullptr) == p.representative &&
          value_of(pk, "OTSWAPV") == value_of(internal::read_keywords(fits), "OTSWAPV"),
          "FITS: the profile's PRODUCT and keywords");
    check(value_of(tk, "PRODUCT") == "complete displacements distance table" &&
          std::strtod(value_of(tk, "ZMIN").c_str(), nullptr) == table.minRedshift() &&
          std::strtod(value_of(tk, "ZMAX").c_str(), nullptr) == table.maxRedshift(),
          "FITS: the table's PRODUCT and its ends");
    const io::Table fp = internal::read_part(fits, "MPSPROF", {"redshift", "MPS"}, {"nTracers"});
    const io::Table ft = internal::read_part(fits, "DISTANCE", {"redshift", "distance", "growthRate"}, {});
    check(fp.nRows == p.redshift.size() && same_bits(column_of(fp, 0), p.redshift) &&
          same_bits(column_of(fp, 1), p.mps) &&
          fp.integers == std::vector<std::int64_t>(p.count.begin(), p.count.end()),
          "FITS: the profile, exactly");
    check(same_bits(column_of(ft, 0), table.redshift()) && same_bits(column_of(ft, 1), table.distance()) &&
          same_bits(column_of(ft, 2), table.growthRate()), "FITS: the distance table, exactly");

    const std::string ascii = temporary("parts.complete.dat");
    const std::string mpsTwin = temporary("parts.complete.mps.dat");
    const std::string tableTwin = temporary("parts.complete.distances.dat");
    check(internal::twin_file(ascii, "mps") == mpsTwin && internal::twin_file(ascii, "distances") == tableTwin,
          "the twins of parts.complete.dat");
    check(internal::twin_file("run", "mps") == "run.mps" &&
          internal::twin_file("a.b/run", "mps") == "a.b/run.mps" &&
          internal::twin_file("a.b/run.txt", "distances") == "a.b/run.distances.txt" &&
          internal::twin_file("run.dat", "mps") == "run.mps.dat",
          "the twin rule: before the extension of the file name, or appended without one");
    io::writeDisplacementComplete(ascii, lc);
    check(value_of(internal::read_keywords(mpsTwin), "PRODUCT") == "complete displacements mps profile" &&
          value_of(internal::read_keywords(tableTwin), "PRODUCT") == "complete displacements distance table",
          "ASCII: the twins, each with its PRODUCT");
    const io::Table ap = io::read(mpsTwin, {"redshift", "MPS", "nTracers"});
    const io::Table at = io::read(tableTwin, {"redshift", "distance", "growthRate"});
    check(same_bits(column_of(ap, 0), p.redshift) && same_bits(column_of(ap, 1), p.mps) &&
          same_bits(column_of(ap, 2), count) && same_bits(column_of(at, 0), table.redshift()) &&
          same_bits(column_of(at, 1), table.distance()) && same_bits(column_of(at, 2), table.growthRate()),
          "ASCII: the profile and the table, exactly, at 17 digits");
    const std::string withTables = text_of(ascii);

    std::remove(mpsTwin.c_str());
    std::remove(tableTwin.c_str());
    io::writeDisplacementComplete(ascii, lc, false);
    io::writeDisplacementComplete(fits, lc, false);
    const std::vector<io::Keyword> without = internal::read_keywords(fits);
    check(!exists(mpsTwin) && !exists(tableTwin) && !internal::has_part(fits, "MPSPROF") &&
          !internal::has_part(fits, "DISTANCE") && value_of(without, "MPSPART") == "missing" &&
          value_of(without, "DISTPART") == "missing", "lightconeTables false: no part, none declared");
    std::string expected = withTables;
    bool declared = true;
    for (const std::string line : {"## MPSPART = 1 / the mps(z) profile is written\n",
                                   "## DISTPART = 1 / the distance table is written\n"}) {
      const std::size_t at = expected.find(line);
      declared = declared && at != std::string::npos;
      if (at != std::string::npos) expected.erase(at, line.size());
    }
    check(declared && text_of(ascii) == expected,
          "lightconeTables false: the same main table, without the two declarations");

    Result noProfile = lc;
    noProfile.mpsProfile = MpsProfile();
    noProfile.distances = std::make_shared<const DistanceTable>(table.redshift(), table.distance());
    io::writeDisplacementComplete(fits, noProfile);
    io::writeDisplacementComplete(ascii, noProfile);
    const std::vector<io::Keyword> partial = internal::read_keywords(fits);
    check(extname_of(fits, 3) == "DISTANCE" && extname_of(fits, 4) == "missing" &&
          value_of(partial, "MPSPART") == "missing" && value_of(partial, "DISTPART") == "1" &&
          !exists(mpsTwin) && exists(tableTwin),
          "a result without a profile: the table alone, in HDU 3");
    check(names_of(tableTwin) == "redshift distance" &&
          throws_naming([&] { internal::read_part(fits, "DISTANCE", {"growthRate"}, {}); }, "growthRate"),
          "a table without growth rate: no growthRate column");
    std::remove(tableTwin.c_str());

    const std::string boxFile = temporary("parts_box.dat");
    io::writeDisplacementComplete(boxFile, box);
    io::writeDisplacementComplete(fits, box);
    check(!exists(internal::twin_file(boxFile, "mps")) && !exists(internal::twin_file(boxFile, "distances")) &&
          extname_of(fits, 3) == "missing", "a box: no part");

    const std::string run = temporary("run.dat");
    const std::string runMps = temporary("run.mps.dat");
    io::writeMpsProfile(runMps, p);
    io::writeDisplacementComplete(run, lc);
    check(value_of(internal::read_keywords(runMps), "PRODUCT") == "complete displacements mps profile",
          "the twin run.mps.dat replaces a writeMpsProfile file of that name");
    for (const std::string& f : {fits, ascii, boxFile, run, runMps, temporary("run.distances.dat")})
      std::remove(f.c_str());
  }

  group("Both displacement writers refuse realization counts or box limits that disagree, before writing");
  {
    const std::string file = temporary("refused.dat");
    std::remove(file.c_str());
    Result counts = box;
    counts.config.nRealizations = 2;
    const std::string message = "the result holds nRealizations = 3 and config.nRealizations = 2; they must be equal";
    check(throws_naming([&] { io::writeDisplacementMean(file, counts); }, message) &&
          throws_naming([&] { io::writeDisplacementMean(file, counts, {io::DisplacementGroup::Id}); }, message) &&
          throws_naming([&] { io::writeDisplacementComplete(file, counts); }, message) && !exists(file),
          "differing realization counts, whatever the groups");
    Result noBox = limited;
    noBox.box.reset();
    Result moved = limited;
    moved.box->upper[1] += 1.;
    check(throws_naming([&] { io::writeDisplacementMean(file, noBox); }, "selection.boxLimits and no box") &&
          throws_naming([&] { io::writeDisplacementComplete(file, noBox); }, "selection.boxLimits and no box") &&
          throws_naming([&] { io::writeDisplacementMean(file, moved); },
                        "selection.boxLimits different from its box") &&
          throws_naming([&] { io::writeDisplacementComplete(file, moved); },
                        "selection.boxLimits different from its box") && !exists(file),
          "limits without the box, or differing from it");
    Result noSky = lc;
    noSky.tracersSky.clear();
    check(throws_naming([&] { io::writeDisplacementComplete(file, noSky); }, "cannot write tracersSky") &&
          !exists(file), "a lightcone without its sky positions");
    Result badProfile = lc;
    badProfile.mpsProfile.count.pop_back();
    check(throws_naming([&] { io::writeDisplacementComplete(file, badProfile); }, "the counts of the profile") &&
          !exists(file), "a malformed profile");
    io::writeDisplacementComplete(file, badProfile, false);
    check(exists(file), "a malformed profile left out is not read");
    std::remove(file.c_str());
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
    check(fits_keyword(boxFile, "TCOMM7") == "'0 corrected, 1 moved by neighbours, 2 no valid neighbour, 3 left out'",
          "the status description in TCOMMn, whole: " + fits_keyword(boxFile, "TCOMM7"));
    check(fits_keyword(boxFile, "AXIS") == "2" && fits_keyword(boxFile, "ZBOX") == "0.5" &&
          fits_keyword(boxFile, "BIAS") == "1.5" && fits_keyword(boxFile, "SIGMA") == "8",
          "box keywords in FITS");
    const io::Table b = io::read(boxFile, {"tracZ"}, {"nValidRec"});
    check(b.nRows == n && same_bits(column_of(b, 0), [&] {
            std::vector<double> z(n);
            for (std::size_t i = 0; i < n; ++i) z[i] = boxCatalog.cartesian[3*i+2];
            return z; }()), "the corrected box positions, exactly");

    io::writeRealSpaceCatalog(boxFile, boxCatalog, {io::CatalogGroup::Id, io::CatalogGroup::Shift});
    const io::Table s = io::read(boxFile, {"shift", "rsdFactor"}, {"id"});
    check(s.values[2*5] == boxCatalog.shift[5] && s.values[2*5+1] == boxCatalog.factor[5] && s.integers[5] == 5,
          "the Shift group");
    check(throws_naming([&] { io::writeRealSpaceCatalog(boxFile, boxCatalog, {io::CatalogGroup::Sky}); },
                        "Sky describes a lightcone"), "Sky on a box is refused");
    std::remove(file.c_str());
    std::remove(boxFile.c_str());
  }

  group("IDs: the id column holds the IDs, or the row without them, in the three writers and both formats");
  {
    std::vector<std::int64_t> ids(n);
    for (std::size_t i = 0; i < n; ++i)
      ids[i] = (i % 2 == 0 ? 1 : -1) * (std::int64_t)(9007199254740993ULL + 2 * i);
    ids[0] = std::numeric_limits<std::int64_t>::min();
    ids[1] = std::numeric_limits<std::int64_t>::max();
    Result lcIds = lc, boxIds = box;
    lcIds.ids = ids;
    boxIds.ids = ids;
    const RealSpaceCatalog lcCatalogIds = realSpaceLightcone(lcIds, table, bias, 10., quiet);
    const RealSpaceCatalog boxCatalogIds = realSpaceBox(boxIds, 2, 0.5, table, 1.5, 8., quiet);
    std::vector<std::int64_t> rows(n);
    for (std::size_t i = 0; i < n; ++i) rows[i] = (std::int64_t)i;

    // The id column of a file, read as 64-bit integers.
    auto read_ids = [] (const std::string& file) {
      const io::Table t = io::read(file, {"tracX"}, {"id"});
      return t.integers;
    };
    for (const std::string extension : {".dat", ".fits"}) {
      const std::string file = temporary("ids" + extension);
      io::writeDisplacementMean(file, lcIds);
      check(read_ids(file) == ids, extension + ": writeDisplacementMean, lightcone, the IDs exactly");
      io::writeDisplacementMean(file, boxIds);
      check(read_ids(file) == ids, extension + ": writeDisplacementMean, box, the IDs exactly");
      io::writeRealSpaceCatalog(file, lcCatalogIds);
      check(read_ids(file) == ids, extension + ": writeRealSpaceCatalog, lightcone, the IDs exactly");
      io::writeRealSpaceCatalog(file, boxCatalogIds);
      check(read_ids(file) == ids, extension + ": writeRealSpaceCatalog, box, the IDs exactly");
      io::writeRealSpaceCatalog(file, boxCatalog, {io::CatalogGroup::Id, io::CatalogGroup::Cartesian});
      check(read_ids(file) == rows, extension + ": writeRealSpaceCatalog without IDs, the rows");
      io::writeDisplacementComplete(file, lcIds);
      const std::vector<std::int64_t> field = read_ids(file);
      bool repeated = field.size() == 3 * n;
      for (std::size_t r = 0; repeated && r < field.size(); ++r) repeated = field[r] == ids[r % n];
      check(repeated, extension + ": writeDisplacementComplete, the IDs in every realization");
      io::writeDisplacementComplete(file, lc);
      const std::vector<std::int64_t> fieldRows = read_ids(file);
      bool counted = fieldRows.size() == 3 * n;
      for (std::size_t r = 0; counted && r < fieldRows.size(); ++r) counted = fieldRows[r] == (std::int64_t)(r % n);
      check(counted, extension + ": writeDisplacementComplete without IDs, the rows");
      if (extension == ".fits")
        check(fits_keyword(file, "TCOMM2") == "'tracer ID; its input row, from 0, when no IDs were set'",
              "the description in TCOMMn, whole: " + fits_keyword(file, "TCOMM2"));
      std::remove(file.c_str());
    }

    const std::string file = temporary("ids.dat");
    io::writeDisplacementMean(file, lcIds);
    check(names_of(file) == "id tracRA tracDec tracRed lagrRA lagrDec lagrRed tracX tracY tracZ lagrX lagrY "
                            "lagrZ displX displY displZ nValidRec outsideRedshiftCut outsideMask",
          "with IDs, id first in the lightcone default: " + names_of(file));
    const std::string header = text_of(file);
    const std::size_t idLine = header.find("\n# id ");
    const std::string line = idLine == std::string::npos ? ""
      : header.substr(idLine + 1, header.find('\n', idLine + 1) - idLine - 1);
    const std::string description = "tracer ID; its input row, from 0, when no IDs were set";
    check(line.size() > description.size() &&
          line.compare(line.size() - description.size(), description.size(), description) == 0,
          "the description in ASCII: " + line);
    io::writeDisplacementMean(file, boxIds);
    check(names_of(file) == "id tracX tracY tracZ lagrX lagrY lagrZ displX displY displZ nValidRec",
          "with IDs, id first in the box default: " + names_of(file));
    io::writeDisplacementMean(file, boxIds, {io::DisplacementGroup::Displacement, io::DisplacementGroup::Id});
    check(names_of(file) == "displX displY displZ id", "with IDs, id requested goes where it is given");
    io::writeDisplacementMean(file, boxIds, {io::DisplacementGroup::Displacement});
    check(names_of(file) == "displX displY displZ", "with IDs, groups given are written as given");
    io::writeRealSpaceCatalog(file, lcCatalogIds);
    check(names_of(file) == "id tracRA tracDec tracRed tracX tracY tracZ nValidRec nNeighbours "
                            "nRealizationsAveraged status",
          "with IDs, id first in the lightcone catalogue default: " + names_of(file));
    io::writeRealSpaceCatalog(file, boxCatalogIds);
    check(names_of(file) == "id tracX tracY tracZ nValidRec nNeighbours nRealizationsAveraged status",
          "with IDs, id first in the box catalogue default: " + names_of(file));

    const std::string length = "ids holds 3 entries for 900 objects; give one per object, or none";
    Result shortIds = box;
    shortIds.ids = {4, 5, 6};
    check(throws_naming([&] { io::writeDisplacementMean(file, shortIds); }, length) &&
          throws_naming([&] { io::writeDisplacementMean(file, shortIds, {io::DisplacementGroup::Displacement}); },
                        length) &&
          throws_naming([&] { io::writeDisplacementComplete(file, shortIds); }, length),
          "IDs of the wrong length are refused, whatever the groups");
    std::vector<std::int64_t> twice = rows;
    twice[300] = twice[10];
    twice[200] = twice[50];
    Result twiceIds = box;
    twiceIds.ids = twice;
    const std::string duplicate = "ids holds 50 twice, at entries 50 and 200; each object needs its own ID";
    check(throws_naming([&] { io::writeDisplacementMean(file, twiceIds); }, duplicate) &&
          throws_naming([&] { io::writeDisplacementMean(file, twiceIds, {io::DisplacementGroup::Displacement}); },
                        duplicate) &&
          throws_naming([&] { io::writeDisplacementComplete(file, twiceIds); }, duplicate),
          "a repeated ID is refused, the first one named, whatever the groups");
    RealSpaceCatalog shortCatalog = boxCatalog, twiceCatalog = boxCatalog;
    shortCatalog.ids = {4, 5, 6};
    twiceCatalog.ids = twice;
    check(throws_naming([&] { io::writeRealSpaceCatalog(file, shortCatalog); }, length) &&
          throws_naming([&] { io::writeRealSpaceCatalog(file, twiceCatalog, {io::CatalogGroup::Shift}); },
                        duplicate), "writeRealSpaceCatalog refuses the same IDs");
    std::remove(file.c_str());
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
