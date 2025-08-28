// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Matching/MatchingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>


#include "TRandom.h"

ActsExamples::MatchingAlgorithm::MatchingAlgorithm(
    ActsExamples::MatchingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("MatchingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
        

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_inputTracksVT.initialize(m_cfg.inputTracksVT);
  m_inputTracksMS.initialize(m_cfg.inputTracksMS);
  m_outputTracksVT.initialize(m_cfg.outputTracksVT);
  m_outputTracksMS.initialize(m_cfg.outputTracksMS);
  m_outputTracksRefit.initialize(m_cfg.outputTracksRefit);
  m_outputMatchedTracks.initialize(m_cfg.outputMatchedTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMapVT.initialize(
      m_cfg.inputMeasurementParticlesMapVT);
  m_inputMeasurementParticlesMapMS.initialize(
      m_cfg.inputMeasurementParticlesMapMS);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
}

ActsExamples::ProcessCode ActsExamples::MatchingAlgorithm::execute(
    const AlgorithmContext& ctx) const {

  const auto& hitParticlesMapMS = m_inputMeasurementParticlesMapMS(ctx);
  const auto& hitParticlesMapVT = m_inputMeasurementParticlesMapVT(ctx);

  const auto& tracksVT = m_inputTracksVT(ctx);
  const auto& tracksMS = m_inputTracksMS(ctx);

  const auto& measurements = m_inputMeasurements(ctx);

  auto pSurfaceMatching = m_cfg.trackingGeometry->findSurface(m_cfg.geoIdForPropagation);

  Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> extrapolator(
      Acts::EigenStepper<>(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  // VOID PROPAGATOR
  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.magneticField);

  Acts::PropagatorPlainOptions pOptions(ctx.geoContext, ctx.magFieldContext);
  
  //////////////
  // MATCHING
  //////////////
  std::vector<std::pair<TrackProxyType, TrackProxyType>> trackPairs;
  int indexMS = 0;
  int indexVT = 0;

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCountsVT;
  std::vector<ParticleHitCount> particleHitCountsMS;

  TrackContainer matchedTracksVT{
      std::make_shared<Acts::VectorTrackContainer>(),
      std::make_shared<Acts::VectorMultiTrajectory>()};

  TrackContainer matchedTracksMS{
      std::make_shared<Acts::VectorTrackContainer>(),
      std::make_shared<Acts::VectorMultiTrajectory>()};

  bool first = true;

  std::vector<int> ntracksVTs;

  for (const auto& trackMS : tracksMS) {
    indexVT = -1;
    std::pair<TrackProxyType, TrackProxyType> trackPairTmp(
        tracksMS.getTrack(indexMS),
        tracksMS.getTrack(indexMS));
    std::pair<int, int> trackPairIndex(indexMS, indexMS);

    identifyContributingParticles(hitParticlesMapMS,
                                  tracksMS.getTrack(indexMS),
                                  particleHitCountsMS);


    if (particleHitCountsMS.empty()) {
      continue;
    }
    /*
    ActsFatras::Barcode majorityParticleIdMS =
        particleHitCountsMS.front().particleId;
    ActsFatras::Barcode majorityParticleIdVT;
    */
    indexMS += 1;

    Acts::BoundTrackParameters params1(
      trackMS.referenceSurface().getSharedPtr(), trackMS.parameters(),
      trackMS.covariance(), trackMS.particleHypothesis());
    
    const auto resMS =
        extrapolator.propagateToSurface(params1, *pSurfaceMatching, pOptions);

    if (!resMS.ok()) {
      continue;
    }
    float chi2 = 1e12;
    const auto& endParamsMS = *resMS;

    Acts::BoundVector paramsMS = endParamsMS.parameters();
    const auto covMS = *endParamsMS.covariance();
    
    int ntracksVT = 0;
    for (const auto& trackVT : tracksVT) {
      ntracksVT++;
      indexVT++;

      Acts::BoundTrackParameters params2(
        trackVT.referenceSurface().getSharedPtr(), trackVT.parameters(),
        trackVT.covariance(), trackVT.particleHypothesis());
      const auto resVT =
          extrapolator.propagateToSurface(params2, *pSurfaceMatching, pOptions);

      identifyContributingParticles(hitParticlesMapVT,
                                    tracksVT.getTrack(indexVT),
                                    particleHitCountsVT);

      if (particleHitCountsVT.empty()) {
        continue;
      }

      if (!resVT.ok()) {
        continue;
      }

      const auto& endParamsVT = *resVT;
      Acts::BoundVector paramsMatch = endParamsVT.parameters() - paramsMS;
      const auto covMatch = *endParamsVT.covariance() + covMS;

      float chi2tmp  = 0;

      for(int i = 0; i < 5; ++i){
        for(int j = 0; j < 5; ++j){
          chi2tmp += paramsMatch(i)*paramsMatch(j)/covMatch(i,j);
        }
      }

      if (chi2tmp < chi2) {
        chi2 = chi2tmp;
        trackPairTmp.second = tracksVT.getTrack(indexVT);
        trackPairIndex.second = indexVT;
        //majorityParticleIdVT = particleHitCountsVT.front().particleId;
      }

    }
    ntracksVTs.push_back(ntracksVT);
    

    if (chi2 != m_cfg.chi2max && trackPairTmp.first != trackPairTmp.second) {
      //if (majorityParticleIdVT == majorityParticleIdMS) {
      trackPairs.push_back(trackPairTmp);
      if (first) {
        matchedTracksVT.ensureDynamicColumns(tracksVT);
        matchedTracksMS.ensureDynamicColumns(tracksMS);
        first = false;
      }

      auto destProxyMS = matchedTracksMS.makeTrack();
      auto srcProxyMS = trackPairTmp.first;
      destProxyMS.copyFrom(srcProxyMS, true);
      destProxyMS.tipIndex() = srcProxyMS.tipIndex();
      destProxyMS.parameters() = endParamsMS.parameters();
      if (endParamsMS.covariance().has_value()) {
        destProxyMS.covariance() = endParamsMS.covariance().value();
      }
    
      auto trackVTMatch = tracksVT.getTrack(trackPairIndex.second);
      Acts::BoundTrackParameters params2Tmp(
        trackVTMatch.referenceSurface().getSharedPtr(), trackVTMatch.parameters(),
        trackVTMatch.covariance(), trackVTMatch.particleHypothesis());
      const auto resVTki =
          extrapolator.propagateToSurface(params2Tmp, *pSurfaceMatching, pOptions);
      auto& matchParamsVT = *resVTki;

      auto destProxyVT = matchedTracksVT.makeTrack();
      auto srcProxyVT = trackPairTmp.second;
      destProxyVT.copyFrom(srcProxyVT, true);
      destProxyVT.tipIndex() = srcProxyVT.tipIndex();
      destProxyVT.parameters() = matchParamsVT.parameters();
      if (matchParamsVT.covariance().has_value()) {
        destProxyVT.covariance() = matchParamsVT.covariance().value();
      }

    }
        
  }
  
  std::vector<Acts::SourceLink> trackSourceLinks;

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  ActsExamples::PassThroughCalibrator pcalibrator;
  ActsExamples::MeasurementCalibratorAdapter calibrator(pcalibrator,
                                                        measurements);

  TrackFitterFunction::GeneralFitterOptions options{
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, pSurfaceMatching,
      pOptions};

  for (auto& pair : trackPairs) {
    auto trackMS = pair.first;
    auto trackVT = pair.second;

    trackSourceLinks.clear();

    const Acts::BoundTrackParameters initialParams(
        trackVT.referenceSurface().getSharedPtr(), trackVT.parameters(),
        trackVT.covariance(), trackMS.particleHypothesis());

    for (auto state : trackMS.trackStatesReversed()) {
      if (!state.hasCalibrated()) {
        continue;
      }
      auto source_link =
          state.getUncalibratedSourceLink().template get<IndexSourceLink>();
      trackSourceLinks.push_back(Acts::SourceLink{source_link});
    }

    for (auto state : trackVT.trackStatesReversed()) {
      if (!state.hasCalibrated()) {
        continue;
      }
      auto source_link =
          state.getUncalibratedSourceLink().template get<IndexSourceLink>();
      trackSourceLinks.push_back(Acts::SourceLink{source_link});
    }
    
    if (trackSourceLinks.empty()) {
      ACTS_WARNING("Empty track found.");
      continue;
    }
    std::reverse(trackSourceLinks.begin(), trackSourceLinks.end());

    auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
                               calibrator, tracks);

    if (result.ok()) {
      // Get the fit output object
      const auto& refittedTrack = result.value();

      if (refittedTrack.hasReferenceSurface()) {
        ACTS_VERBOSE("Refitted parameters for track ");
      } else {
        ACTS_DEBUG("No refitted parameters for track ");
      }
    } else {
      ACTS_WARNING("Fit failed for track with error: "
                   << result.error() << ", " << result.error().message());
    }
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  m_outputTracksRefit(ctx, std::move(constTracks));

  m_outputMatchedTracks(ctx, std::move(trackPairs));

  ActsExamples::ConstTrackContainer outputTracksVT{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(matchedTracksVT.container())),
      tracksVT.trackStateContainerHolder()};

  m_outputTracksVT(ctx, std::move(outputTracksVT));

  ActsExamples::ConstTrackContainer outputTracksMS{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(matchedTracksMS.container())),
      tracksMS.trackStateContainerHolder()};

  m_outputTracksMS(ctx, std::move(outputTracksMS));

  return ActsExamples::ProcessCode::SUCCESS;
}