// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using TrackProxyType = Acts::TrackContainer<Acts::ConstVectorTrackContainer,
                                            Acts::ConstVectorMultiTrajectory,
                                            std::shared_ptr>::ConstTrackProxy;

namespace ActsExamples {
struct AlgorithmContext;

class MatchingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Optional. Input track parameters collection
    std::string inputTracksMS;
    /// Optional. Input track parameters collection
    std::string inputTracksVT;
    /// Optional. Input track parameters collection
    std::string inputMeasurements = "measurements";

    std::string outputTracksVT = "matchVT";
    std::string outputTracksMS = "matchMS";
    std::string outputTracksRefit = "refit";
    std::string outputMatchedTracks = "matched";
    std::string inputMeasurementParticlesMapVT;
    std::string inputMeasurementParticlesMapMS;
    double chi2max = 1e12;
    std::uint64_t geoIdForPropagation;
    std::shared_ptr<ActsExamples::TrackFitterFunction> fit;

    /// The tracking geometry that should be used.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// The magnetic field that should be used.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
  };

  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  MatchingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the ambiguity resolution algorithm.
  ///
  /// @param cxt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

  inline std::vector<Acts::InputTrack> makeInputTracks(
      const TrackParametersContainer& trackParameters) const {
    std::vector<Acts::InputTrack> inputTracks;
    inputTracks.reserve(trackParameters.size());

    for (const auto& trackParam : trackParameters) {
      inputTracks.emplace_back(&trackParam);
    }
    return inputTracks;
  }

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracksVT{this, "InputTracksVT"};
  ReadDataHandle<ConstTrackContainer> m_inputTracksMS{this, "InputTracksMS"};

  ReadDataHandle<IndexMultimap<ActsFatras::Barcode>>
      m_inputMeasurementParticlesMapVT{this, "InputMeasurementParticlesMapVT"};
  ReadDataHandle<IndexMultimap<ActsFatras::Barcode>>
      m_inputMeasurementParticlesMapMS{this, "InputMeasurementParticlesMapMS"};

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  WriteDataHandle<std::vector<std::pair<TrackProxyType, TrackProxyType>>>
      m_outputMatchedTracks{this, "OutputMatchedTracks"};

  WriteDataHandle<ConstTrackContainer> m_outputTracksMS{this, "OutputRefitVT"};
  WriteDataHandle<ConstTrackContainer> m_outputTracksVT{this, "OutputRefitMS"};
  WriteDataHandle<ConstTrackContainer> m_outputTracksRefit{this, "OutputRefit"};
};

}  // namespace ActsExamples
