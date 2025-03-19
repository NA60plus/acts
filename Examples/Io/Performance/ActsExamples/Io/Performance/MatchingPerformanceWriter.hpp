// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/MatchingPlotTool.hpp"

#include <cstddef>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {
struct AlgorithmContext;

class MatchingPerformanceWriter final : public WriterT<ConstTrackContainer> {
 public:
  struct Config {

    std::string inputTracks;
    std::string inputRecPrimaryVertex = "OutputRecPrimaryVertex";
    std::string inputGenPrimaryVertex = "OutputGenPrimaryVertex";
    std::string inputFitFunction = "OutputFitFunction";
    std::string inputZTracklets = "OutputZTracklets";
    std::string inputZTrackletsPeak = "OutputZTrackletsPeak";

    /// Output filename.
    std::string filePath = "performance_matching.root";
    /// Output file mode
    std::string fileMode = "RECREATE";
    bool verbose=false;
    MatchingPlotTool::Config matchingPlotToolConfig;
  };

  /// Construct from configuration and log level.
  /// @param config The configuration
  /// @param level
  MatchingPerformanceWriter(Config config, Acts::Logging::Level level);

  ~MatchingPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx, const ConstTrackContainer& track);

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};

  TDirectory *subdirHist{nullptr};
  TDirectory *subdirPeak{nullptr};

  /// Plot tool for efficiency
  MatchingPlotTool m_matchingPlotTool;
  MatchingPlotTool::MatchingPlotCache m_matchingPlotCache;
  int ev_counter;

  ReadDataHandle<double> m_inputRecPrimaryVertex{this, "OutputRecPrimaryVertex"};
  ReadDataHandle<double> m_inputGenPrimaryVertex{this, "OutputGenPrimaryVertex"};
  ReadDataHandle<std::vector<double>> m_inputFitFunction{this, "OutputFitFuncVtx"};
  ReadDataHandle<std::vector<double>> m_inputZTracklets{this, "OutputZTracklets"};
  ReadDataHandle<std::vector<double>> m_inputZTrackletsPeak{this, "OutputZTrackletsPeak"};

};

}  // namespace ActsExamples
