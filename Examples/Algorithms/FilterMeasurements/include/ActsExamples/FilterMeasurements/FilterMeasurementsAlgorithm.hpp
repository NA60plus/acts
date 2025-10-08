// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ActsExamples {
struct AlgorithmContext;

class FilterMeasurementsAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input trajectories collection.
    std::string inputTracks;
    /// Output measurements collection.
    std::string outputMeasurements;
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap ;
    /// Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap;
    /// Output collection to map particles to measurements.
    std::string outputParticleMeasurementsMap ;
    /// Output collection to map particles to simulated hits.
    std::string outputSimHitMeasurementsMap;
    /// Output collection to map particles to simulated hits.
    std::string inputSimHitMeasurementsMap;
  };

  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  FilterMeasurementsAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the ambiguity resolution algorithm.
  ///
  /// @param cxt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  WriteDataHandle<MeasurementContainer> m_outputMeasurements{this, "OutputMeasurements"};
  ReadDataHandle<SimHitContainer> m_inputHits{this, "InputHits"};

  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementSimHitsMap{
    this, "InputMeasurementSimHitsMap"};

  WriteDataHandle<IndexMultimap<SimBarcode>> m_outputMeasurementParticlesMap{
      this, "OutputMeasurementParticlesMap"};
  WriteDataHandle<IndexMultimap<Index>> m_outputMeasurementSimHitsMap{
      this, "OutputMeasurementSimHitsMap"};

  WriteDataHandle<InverseMultimap<SimBarcode>> m_outputParticleMeasurementsMap{
      this, "OutputParticleMeasurementsMap"};
  WriteDataHandle<InverseMultimap<Index>> m_outputSimHitMeasurementsMap{
      this, "OutputSimHitMeasurementsMap"};


};

}  // namespace ActsExamples
