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
 *  @file app/otswap.cpp
 *
 *  @brief The otswap executable: a reconstruction, and optionally its
 *  redshift-space correction, run from a parameter file.
 *
 *  Usage: otswap <parameter-file>, otswap --help, otswap --version. Exit
 *  status 0 when every requested file is written, 1 on any error, with
 *  one line per error on stderr, "otswap: error: <message>". Every check
 *  that needs no data comes before any data is read. The keys are listed
 *  in docs/executable.md.
 *
 *  Uses the public headers only, as an outside program would.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

#include "otswap/OT.h"
#include "otswap/ParameterFile.h"
#include "otswap/RSD.h"
#include "otswap/io.h"

namespace {

  const char* const kUsage = "usage: otswap <parameter-file>";

  const char* const kHelp =
    "usage: otswap <parameter-file>\n"
    "  --help     print this help and exit\n"
    "  --version  print the version and exit\n"
    "The keys of the parameter file are listed in docs/executable.md.\n";

  // The defaults of the DistanceTable cosmology constructor, and the
  // redshift range of the table it builds with its other defaults.
  constexpr double kDefaultW0 = -1.;
  constexpr double kDefaultWa = 0.;
  constexpr double kTableMaxRedshift = 10.;

  // ==========================================================================
  // The keys, and the runs each applies to
  // ==========================================================================

  /// When a key applies to a run of one geometry: never, always, or only
  /// with a switch on.
  enum class Applies { Never, Always, PerformRecons, RandomFile, RejectCrossings };

  struct Key {
    const char* name;
    Applies box;
    Applies lightcone;
  };

  /// Every key, in the order of docs/executable.md, which is also the
  /// order of the typed reads.
  const Key kKeys[] = {
    {"geometry",                      Applies::Always,        Applies::Always},
    {"verbosity",                     Applies::Always,        Applies::Always},
    {"tracer_file",                   Applies::Always,        Applies::Always},
    {"column_tracer",                 Applies::Always,        Applies::Always},
    {"id_column",                     Applies::Always,        Applies::Always},
    {"random_file",                   Applies::Always,        Applies::Always},
    {"column_random",                 Applies::RandomFile,    Applies::Always},
    {"units",                         Applies::Never,         Applies::Always},
    {"delimiter",                     Applies::Always,        Applies::Always},
    {"comment",                       Applies::Always,        Applies::Always},
    {"veto_mask_filename",            Applies::Never,         Applies::Always},
    {"z_min",                         Applies::Never,         Applies::Always},
    {"z_max",                         Applies::Never,         Applies::Always},
    {"Omega_matter",                  Applies::PerformRecons, Applies::Always},
    {"little_h",                      Applies::PerformRecons, Applies::Always},
    {"w0",                            Applies::PerformRecons, Applies::Always},
    {"wa",                            Applies::PerformRecons, Applies::Always},
    {"Nrec",                          Applies::Always,        Applies::Always},
    {"threshold",                     Applies::Always,        Applies::Always},
    {"seed",                          Applies::Always,        Applies::Always},
    {"mesh_cellsize",                 Applies::Always,        Applies::Always},
    {"bins",                          Applies::Never,         Applies::Always},
    {"reject_mask_crossings",         Applies::Never,         Applies::Always},
    {"max_unobserved_pixels_crossed", Applies::Never,         Applies::RejectCrossings},
    {"box_min",                       Applies::Always,        Applies::Never},
    {"box_max",                       Applies::Always,        Applies::Never},
    {"mps",                           Applies::Always,        Applies::Never},
    {"perform_recons",                Applies::Always,        Applies::Always},
    {"rsd_smoothing_scale",           Applies::PerformRecons, Applies::PerformRecons},
    {"rsd_weight_by_realizations",    Applies::PerformRecons, Applies::PerformRecons},
    {"redshift",                      Applies::PerformRecons, Applies::Never},
    {"bias",                          Applies::PerformRecons, Applies::Never},
    {"los_axis",                      Applies::PerformRecons, Applies::Never},
    {"bias_file",                     Applies::Never,         Applies::PerformRecons},
    {"output_dir",                    Applies::Always,        Applies::Always},
    {"output_file",                   Applies::Always,        Applies::Always},
    {"output_format",                 Applies::Always,        Applies::Always},
    {"printDisplMean",                Applies::Always,        Applies::Always},
    {"printDisplComplete",            Applies::Always,        Applies::Always},
    {"printMps",                      Applies::Never,         Applies::Always},
    {"printRecCat",                   Applies::PerformRecons, Applies::PerformRecons},
  };

  /// The switches that decide which keys apply, read before the others.
  struct Switches {
    bool lightcone = false;
    bool performRecons = false;
    bool randomFile = false;
    bool rejectCrossings = false;
  };

  /// Every key of the file that does not apply to this run, and every
  /// unknown key, in file order, each as its message.
  std::vector<std::string> irrelevant_keys (const otswap::io::ParameterFile& parameters,
                                            const std::string& file, const Switches& run)
  {
    std::vector<std::string> messages;
    for (const otswap::io::ParameterFile::Entry& entry : parameters.unread()) {
      const std::string at = file + ":" + std::to_string(entry.line) + ": key '" + entry.key + "': ";
      const Key* key = nullptr;
      for (const Key& k : kKeys)
        if (entry.key == k.name) key = &k;
      if (key == nullptr) {
        messages.push_back(at + "unknown key");
        continue;
      }
      switch (run.lightcone ? key->lightcone : key->box) {
        case Applies::Never:
          messages.push_back(at + "does not apply to a " + (run.lightcone ? "lightcone" : "box"));
          break;
        case Applies::Always:
          break;
        case Applies::PerformRecons:
          if (!run.performRecons) messages.push_back(at + "does not apply to a run without perform_recons");
          break;
        case Applies::RandomFile:
          if (!run.randomFile) messages.push_back(at + "does not apply to a run without random_file");
          break;
        case Applies::RejectCrossings:
          if (!run.rejectCrossings)
            messages.push_back(at + "does not apply to a run without reject_mask_crossings");
          break;
      }
    }
    return messages;
  }

  // ==========================================================================
  // Typed reads, with the executable's own checks
  // ==========================================================================

  /// The parameter file with the helpers of the executable's messages,
  /// which name the file, the line and the key, and quote a value as
  /// written in the file.
  class Parameters {
  public:

    Parameters (const otswap::io::ParameterFile& parameters, const std::string& file)
      : m_parameters(parameters), m_file(file) {}

    const otswap::io::ParameterFile& file () const { return m_parameters; }

    std::string at (const std::string& key) const
    {
      return m_file + ":" + std::to_string(m_parameters.line(key)) + ": key '" + key + "': ";
    }

    std::string line (const std::string& key) const { return std::to_string(m_parameters.line(key)); }

    [[noreturn]] void fail (const std::string& key, const std::string& message) const
    {
      throw otswap::Error(at(key) + message);
    }

    /// The value of a scalar key as written, read after its typed read.
    std::string written (const std::string& key) const { return m_parameters.get<std::string>(key); }

    /// Element i, from 0, of a list key as written.
    std::string written (const std::string& key, const std::size_t i) const
    {
      return m_parameters.get<std::vector<std::string>>(key)[i];
    }

    /// X6: a value outside the range the library accepts.
    void require (const std::string& key, const bool inRange, const std::string& rule) const
    {
      if (!inRange) fail(key, written(key) + " is out of range; it must be " + rule);
    }

    /// X5: a list key of three elements.
    template <typename T>
    std::vector<T> three (const std::string& key) const
    {
      const std::vector<T> values = m_parameters.get<std::vector<T>>(key);
      if (values.size() != 3)
        fail(key, "expected 3 values, found " + std::to_string(values.size()));
      return values;
    }

    /// A key whose value is one of a few words (X1-X4); without a default
    /// the key is required.
    std::string word (const std::string& key, const std::vector<std::string>& words,
                      const char* defaultWord = nullptr) const
    {
      const std::string value = defaultWord == nullptr ? m_parameters.get<std::string>(key)
                                                       : m_parameters.get<std::string>(key, defaultWord);
      if (std::find(words.begin(), words.end(), value) != words.end()) return value;
      std::string list;
      for (std::size_t i = 0; i < words.size(); ++i)
        list += (i == 0 ? "" : i+1 == words.size() ? " or " : ", ") + words[i];
      fail(key, value + " is not " + list);
    }

  private:

    const otswap::io::ParameterFile& m_parameters;
    const std::string& m_file;
  };

  /// True if io::read reads the file as FITS: the rule of its extension.
  bool is_fits (const std::string& file)
  {
    auto ends_with = [&file] (const std::string& tail) {
      if (file.size() < tail.size()) return false;
      for (std::size_t i = 0; i < tail.size(); ++i) {
        char c = file[file.size()-tail.size()+i];
        if (c >= 'A' && c <= 'Z') c = char(c-'A'+'a');
        if (c != tail[i]) return false;
      }
      return true;
    };
    return ends_with(".fits") || ends_with(".fit") || ends_with(".fits.gz");
  }

  bool all_digits (const std::string& s)
  {
    return !s.empty() && std::all_of(s.begin(), s.end(), [] (const char c) { return c >= '0' && c <= '9'; });
  }

  /// X9: true if two column selectors of io::read name the same column,
  /// as far as that can be told without reading the file. In ASCII two
  /// indices are equal when their values are, two names when they are
  /// identical; an index and a name cannot be compared here. In FITS every
  /// selector is a name, compared ignoring case.
  bool same_column (const std::string& a, const std::string& b, const bool fits)
  {
    auto lower = [] (std::string s) {
      for (char& c : s) if (c >= 'A' && c <= 'Z') c = char(c-'A'+'a');
      return s;
    };
    if (fits) return lower(a) == lower(b);
    if (all_digits(a) && all_digits(b)) {
      auto strip = [] (const std::string& s) {
        const std::size_t first = s.find_first_not_of('0');
        return first == std::string::npos ? std::string("0") : s.substr(first);
      };
      return strip(a) == strip(b);
    }
    return a == b;
  }

  /// Everything the parameter file sets, typed and checked.
  struct Run {
    Switches switches;
    otswap::Config config;
    otswap::CorrectionConfig correction;

    std::string tracerFile;
    std::vector<std::string> tracerColumns;
    std::optional<std::string> idColumn;
    std::optional<std::string> randomFile;
    std::vector<std::string> randomColumns;
    bool degrees = false;
    char delimiter = ' ';
    char comment = '#';

    std::string maskFile;
    otswap::RedshiftCut cut;

    double omegaM = 0., h = 0., w0 = kDefaultW0, wa = kDefaultWa;

    unsigned bins = 0;
    std::optional<otswap::BoxLimits> limits;
    std::optional<double> mps;

    double sigma = 0.;
    double redshift = 0., bias = 0.;
    unsigned axis = 0;
    std::string biasFile;

    std::string outputDir, outputName, extension;
    bool printMean = false, printComplete = false, printMps = false, printCatalog = false;
  };

  /// The typed reads and the value checks, in the order of the keys, the
  /// first error stopping the run; then the checks of the outputs as a
  /// whole. Reads only the keys that apply to the run.
  Run read_run (const Parameters& p, const std::string& file, const Switches& switches)
  {
    const otswap::io::ParameterFile& f = p.file();
    const bool lightcone = switches.lightcone;
    const bool recons = switches.performRecons;

    Run run;
    run.switches = switches;

    if (f.has("verbosity")) {
      const std::string v = p.word("verbosity", {"silent", "normal", "detailed"});
      run.config.verbosity = v == "silent" ? otswap::Verbosity::Silent
                           : v == "normal" ? otswap::Verbosity::Normal : otswap::Verbosity::Detailed;
    }
    run.correction.verbosity = run.config.verbosity;

    // Input
    run.tracerFile = f.get<std::string>("tracer_file");
    run.tracerColumns = p.three<std::string>("column_tracer");
    if (f.has("id_column")) {
      run.idColumn = f.get<std::string>("id_column");
      for (const std::string& c : run.tracerColumns)
        if (same_column(*run.idColumn, c, is_fits(run.tracerFile)))
          p.fail("id_column", *run.idColumn + " is also a column of column_tracer, at line " +
                 p.line("column_tracer"));
    }
    if (lightcone || switches.randomFile) {
      run.randomFile = f.get<std::string>("random_file");
      run.randomColumns = p.three<std::string>("column_random");
    }
    if (lightcone)
      run.degrees = p.word("units", {"degrees", "radians"}) == "degrees";
    run.delimiter = f.get<char>("delimiter", run.delimiter);
    run.comment = f.get<char>("comment", run.comment);

    // Lightcone selection
    if (lightcone) {
      run.maskFile = f.get<std::string>("veto_mask_filename");
      run.cut.min = f.get<double>("z_min", run.cut.min);
      run.cut.max = f.get<double>("z_max", run.cut.max);
      if (f.has("z_min") && f.has("z_max") && run.cut.max < run.cut.min)
        p.fail("z_max", p.written("z_max") + " is below z_min, " + p.written("z_min") + ", at line " +
               p.line("z_min"));
    }

    // Cosmology
    if (lightcone || recons) {
      run.omegaM = f.get<double>("Omega_matter");
      p.require("Omega_matter", run.omegaM > 0. && run.omegaM <= 1., "in (0, 1]");
      run.h = f.get<double>("little_h");
      p.require("little_h", run.h > 0., "positive");
      run.w0 = f.get<double>("w0", run.w0);
      run.wa = f.get<double>("wa", run.wa);
    }

    // Reconstruction
    otswap::Config& c = run.config;
    c.nRealizations = f.get<unsigned>("Nrec", c.nRealizations);
    p.require("Nrec", c.nRealizations >= 1, "at least 1");
    c.convergence = f.get<double>("threshold", c.convergence);
    p.require("threshold", c.convergence > 0. && c.convergence < 1., "strictly between 0 and 1");
    c.seed = f.get<unsigned>("seed", c.seed);
    c.cellSize = f.get<double>("mesh_cellsize", c.cellSize);
    p.require("mesh_cellsize", c.cellSize > 0., "positive");
    if (lightcone) {
      run.bins = f.get<unsigned>("bins");
      p.require("bins", run.bins >= 1, "at least 1");
      c.rejectCrossings = switches.rejectCrossings;
      if (c.rejectCrossings)
        c.maxUnobservedPixelsCrossed = f.get<unsigned>("max_unobserved_pixels_crossed",
                                                       c.maxUnobservedPixelsCrossed);
    }
    else {
      std::vector<double> lower, upper;
      if (f.has("box_min")) lower = p.three<double>("box_min");
      if (f.has("box_max")) upper = p.three<double>("box_max");
      if (f.has("box_min") != f.has("box_max")) {
        const bool min = f.has("box_min");
        p.fail(min ? "box_min" : "box_max", std::string(min ? "box_max" : "box_min") +
               " is missing; give both or neither");
      }
      if (f.has("box_min")) {
        otswap::BoxLimits limits;
        for (std::size_t k = 0; k < 3; ++k) {
          if (!(upper[k] > lower[k]))
            p.fail("box_max", "element " + std::to_string(k+1) + ", " + p.written("box_max", k) +
                   ", is not above element " + std::to_string(k+1) + " of box_min, " +
                   p.written("box_min", k) + ", at line " + p.line("box_min"));
          limits.lower[k] = lower[k];
          limits.upper[k] = upper[k];
        }
        run.limits = limits;
      }
      if (f.has("mps")) {
        run.mps = f.get<double>("mps");
        p.require("mps", *run.mps > 0., "positive");
      }
    }

    // RSD correction
    if (recons) {
      run.sigma = f.get<double>("rsd_smoothing_scale");
      p.require("rsd_smoothing_scale", run.sigma >= 0., "non-negative");
      run.correction.weightByRealizations = f.get<bool>("rsd_weight_by_realizations",
                                                        run.correction.weightByRealizations);
      if (lightcone)
        run.biasFile = f.get<std::string>("bias_file");
      else {
        run.redshift = f.get<double>("redshift");
        p.require("redshift", run.redshift >= 0. && run.redshift <= kTableMaxRedshift,
                  "in [0, 10], the range of the distance table");
        run.bias = f.get<double>("bias");
        p.require("bias", run.bias > 0., "positive");
        run.axis = f.get<unsigned>("los_axis");
        p.require("los_axis", run.axis <= 2, "0, 1 or 2");
      }
    }

    // Output
    run.outputDir = f.get<std::string>("output_dir");
    run.outputName = f.get<std::string>("output_file");
    if (run.outputName.empty() || run.outputName.find('/') != std::string::npos)
      p.fail("output_file", "\"" + run.outputName + "\" is not a file name; it must be non-empty and "
             "hold no '/'");
    run.extension = p.word("output_format", {"dat", "fits"}, "dat");
    run.printMean = f.get<bool>("printDisplMean", false);
    run.printComplete = f.get<bool>("printDisplComplete", false);
    if (lightcone) run.printMps = f.get<bool>("printMps", false);
    if (recons) run.printCatalog = f.get<bool>("printRecCat", false);

    if (recons && !run.printCatalog)
      p.fail("perform_recons", "true needs printRecCat = true; the corrected catalogue is the "
             "correction's only output");
    if (!run.printMean && !run.printComplete && !run.printMps && !run.printCatalog)
      throw otswap::Error(file + ": no output requested; set " +
                          (lightcone ? "printDisplMean, printDisplComplete or printMps"
                                     : "printDisplMean or printDisplComplete") + " to true");

    return run;
  }

  // ==========================================================================
  // The run
  // ==========================================================================

  /// X14: the first ID met twice, in row order, with both rows from 1.
  void check_ids (const std::vector<std::int64_t>& ids, const Run& run)
  {
    std::unordered_map<std::int64_t, std::size_t> first;
    for (std::size_t row = 0; row < ids.size(); ++row) {
      const auto found = first.emplace(ids[row], row);
      if (!found.second)
        throw otswap::Error(run.tracerFile + ": the column " + *run.idColumn + " holds the ID " +
                            std::to_string(ids[row]) + " twice, at data rows " +
                            std::to_string(found.first->second+1) + " and " + std::to_string(row+1));
    }
  }

  int run_file (const std::string& file)
  {
    const otswap::io::ParameterFile parameters(file);
    const Parameters p(parameters, file);

    // The switches that decide which keys apply
    Switches switches;
    switches.lightcone = p.word("geometry", {"box", "lightcone"}) == "lightcone";
    switches.performRecons = parameters.get<bool>("perform_recons", false);
    switches.randomFile = parameters.has("random_file");
    if (switches.lightcone)
      switches.rejectCrossings = parameters.get<bool>("reject_mask_crossings",
                                                      otswap::Config{}.rejectCrossings);

    // Every key that does not apply, reported at once
    std::vector<std::string> irrelevant = irrelevant_keys(parameters, file, switches);
    if (!irrelevant.empty()) {
      for (const std::string& message : irrelevant) std::cerr << "otswap: error: " << message << '\n';
      return 1;
    }

    const Run run = read_run(p, file, switches);

    // Every key that applies has been read; one left is a bug here, reported
    // as an unknown key
    irrelevant.clear();
    for (const otswap::io::ParameterFile::Entry& entry : parameters.unread())
      irrelevant.push_back(file + ":" + std::to_string(entry.line) + ": key '" + entry.key + "': unknown key");
    if (!irrelevant.empty()) {
      for (const std::string& message : irrelevant) std::cerr << "otswap: error: " << message << '\n';
      return 1;
    }

    // The output directory, created before the long steps
    const std::filesystem::path directory(run.outputDir);
    std::error_code error;
    const std::filesystem::file_status status = std::filesystem::status(directory, error);
    if (std::filesystem::exists(status) && !std::filesystem::is_directory(status))
      p.fail("output_dir", "cannot create the directory " + run.outputDir + ": not a directory");
    std::filesystem::create_directories(directory, error);
    if (error)
      p.fail("output_dir", "cannot create the directory " + run.outputDir + ": " + error.message());

    // The objects built from the parameters, before the data
    std::optional<otswap::DistanceTable> distances;
    if (run.switches.lightcone || run.switches.performRecons)
      distances.emplace(run.omegaM, run.h, run.w0, run.wa);
    std::optional<otswap::Mask> mask;
    if (run.switches.lightcone) mask.emplace(run.maskFile);
    std::optional<otswap::BiasTable> biasTable;
    if (run.switches.lightcone && run.switches.performRecons)
      biasTable = otswap::io::readBiasTable(run.biasFile, run.delimiter, run.comment);

    // The data
    std::vector<double> tracers, randoms;
    std::vector<std::int64_t> ids;
    if (run.idColumn) {
      otswap::io::Table table = otswap::io::read(run.tracerFile, run.tracerColumns, {*run.idColumn},
                                                 run.delimiter, run.comment);
      tracers = std::move(table.values);
      ids = std::move(table.integers);
      check_ids(ids, run);
    }
    else
      tracers = otswap::io::read(run.tracerFile, run.tracerColumns, run.delimiter, run.comment).values;
    if (run.randomFile)
      randoms = otswap::io::read(*run.randomFile, run.randomColumns, run.delimiter, run.comment).values;
    if (run.switches.lightcone && run.degrees) {
      tracers = otswap::skyToRadians(std::move(tracers));
      randoms = otswap::skyToRadians(std::move(randoms));
    }

    // The reconstruction
    otswap::Result result =
      run.switches.lightcone
        ? otswap::reconstructLightcone(tracers, randoms, *mask, run.bins, *distances, run.config, run.cut)
        : run.mps ? otswap::reconstructBox(tracers, randoms, *run.mps, run.config, run.limits)
                  : otswap::reconstructBox(tracers, randoms, run.config, run.limits);
    if (run.idColumn) result.ids = std::move(ids);

    // Its files, before the correction, so that a correction that fails
    // leaves them on disk
    const bool report = run.config.verbosity != otswap::Verbosity::Silent;
    const std::string zspace = run.switches.performRecons ? "zspace_" : "";
    auto path = [&] (const std::string& prefix, const std::string& inserted = "") {
      return (directory / (prefix + run.outputName + inserted + "." + run.extension)).string();
    };
    auto wrote = [report] (const std::string& written) {
      if (report) std::clog << "otswap: wrote " << written << std::endl;
    };

    if (run.printMean) {
      otswap::io::writeDisplacementMean(path("displacement_mean_" + zspace), result);
      wrote(path("displacement_mean_" + zspace));
    }
    if (run.printComplete) {
      otswap::io::writeDisplacementComplete(path("displacement_complete_" + zspace), result, true);
      wrote(path("displacement_complete_" + zspace));
      if (run.switches.lightcone && run.extension == "dat") {
        wrote(path("displacement_complete_" + zspace, ".mps"));
        wrote(path("displacement_complete_" + zspace, ".distances"));
      }
    }
    if (run.printMps) {
      otswap::io::writeMpsProfile(path("MPS_" + zspace), result.mpsProfile);
      wrote(path("MPS_" + zspace));
    }

    // The correction
    if (run.switches.performRecons) {
      const otswap::RealSpaceCatalog catalog =
        run.switches.lightcone
          ? otswap::realSpaceLightcone(result, *distances, *biasTable, run.sigma, run.correction)
          : otswap::realSpaceBox(result, run.axis, run.redshift, *distances, run.bias, run.sigma,
                                 run.correction);
      otswap::io::writeRealSpaceCatalog(path("reconstructed_catalogue_"), catalog);
      wrote(path("reconstructed_catalogue_"));
    }

    return 0;
  }

}


// ============================================================================


int main (int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << kUsage << '\n';
    return 1;
  }
  const std::string argument = argv[1];
  if (argument == "--help") {
    std::cout << kHelp;
    return 0;
  }
  if (argument == "--version") {
    std::cout << "otswap " << OTSWAP_VERSION << '\n';
    return 0;
  }
  if (argument.compare(0, 2, "--") == 0) {
    std::cerr << kUsage << '\n';
    return 1;
  }

  try {
    return run_file(argument);
  }
  catch (const std::exception& e) {
    std::cerr << "otswap: error: " << e.what() << '\n';
  }
  catch (...) {
    std::cerr << "otswap: error: unknown exception\n";
  }
  return 1;
}
