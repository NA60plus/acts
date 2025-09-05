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
        

  m_inputTracksVT.initialize(m_cfg.inputTracksVT);
  m_inputTracksMS.initialize(m_cfg.inputTracksMS);
  m_outputTracksVT.initialize(m_cfg.outputTracksVT);
  m_outputTracksMS.initialize(m_cfg.outputTracksMS);
  m_outputTracksRefit.initialize(m_cfg.outputTracksRefit);
  m_outputMatchedTracks.initialize(m_cfg.outputMatchedTracks);
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

  Acts::GeometryIdentifier geoIdForPropagation(m_cfg.geoIdForPropagation);
  auto pSurfaceMatching = m_cfg.trackingGeometry->findSurface(geoIdForPropagation);

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

  ACTS_DEBUG("Start matching " << tracksMS.size() << " MS tracks with "
                              << tracksVT.size() << " VT tracks");
  for (const auto& trackMS : tracksMS) {
    ACTS_DEBUG("Matching MS track " << indexMS << " / "
                                   << tracksMS.size());
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
    /*
    const auto resMS =
        extrapolator.propagateToSurface(params1, *pSurfaceMatching, pOptions);

    if (!resMS.ok()) {
      continue;
    }
    const auto& endParamsMS = *resMS;
    */

    float chi2 = 1e12;
    Acts::BoundVector paramsMS = trackMS.parameters();
    const auto covMS = trackMS.covariance();
    
    int ntracksVT = 0;
    for (const auto& trackVT : tracksVT) {
      ntracksVT++;
      indexVT++;
      ACTS_DEBUG("Matching VS track " << indexVT << " / "
                                    << tracksVT.size());

      Acts::BoundTrackParameters params2(
        trackVT.referenceSurface().getSharedPtr(), trackVT.parameters(),
        trackVT.covariance(), trackVT.particleHypothesis());

      identifyContributingParticles(hitParticlesMapVT,
                                    tracksVT.getTrack(indexVT),
                                    particleHitCountsVT);

      if (particleHitCountsVT.empty()) {
        continue;
      }

      Acts::BoundVector paramsMatch = trackVT.parameters() - paramsMS;
      const auto covMatch = trackVT.covariance() + covMS;

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
      destProxyMS.parameters() = trackMS.parameters();
      destProxyMS.covariance() = trackMS.covariance();
    
      auto trackVTMatch = tracksVT.getTrack(trackPairIndex.second);
      auto destProxyVT = matchedTracksVT.makeTrack();
      auto srcProxyVT = trackPairTmp.second;
      destProxyVT.copyFrom(srcProxyVT, true);
      destProxyVT.tipIndex() = srcProxyVT.tipIndex();
      destProxyVT.parameters() = trackVTMatch.parameters();
      destProxyVT.covariance() = trackVTMatch.covariance(); 

    }
        
  }

  ACTS_DEBUG("Number of matched tracks: " << trackPairs.size());
  
  std::vector<Acts::SourceLink> trackSourceLinks;
  std::vector<const Acts::Surface*> surfSequence;


  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  ActsExamples::PassThroughCalibrator pcalibrator;
//ActsExamples::MeasurementCalibratorAdapter calibrator(pcalibrator,
//                                                        measurements);
  RefittingCalibrator calibrator;

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});


  for (auto& pair : trackPairs) {
    auto trackMS = pair.first;
    auto trackVT = pair.second;

    ACTS_DEBUG("Refitting matched track with tip index "
               << trackMS.tipIndex() << " (MS) and "
               << trackVT.tipIndex() << " (VT)");

    trackSourceLinks.clear();

    TrackFitterFunction::GeneralFitterOptions options{
        ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
        &trackVT.referenceSurface(),
        Acts::PropagatorPlainOptions(ctx.geoContext, ctx.magFieldContext)};
    options.doRefit = true;

    const Acts::BoundTrackParameters initialParams(
        trackVT.referenceSurface().getSharedPtr(), trackVT.parameters(),
        trackVT.covariance(), trackMS.particleHypothesis());

    for (auto state : trackMS.trackStatesReversed()) {
      surfSequence.push_back(&state.referenceSurface());

      if (!state.hasCalibrated()) {
        continue;
      }

      auto sl = RefittingCalibrator::RefittingSourceLink{state};
      trackSourceLinks.push_back(Acts::SourceLink{sl});
    }

    for (auto state : trackVT.trackStatesReversed()) {
      surfSequence.push_back(&state.referenceSurface());

      if (!state.hasCalibrated()) {
        continue;
      }

      auto sl = RefittingCalibrator::RefittingSourceLink{state};
      trackSourceLinks.push_back(Acts::SourceLink{sl});
    }

    if (trackSourceLinks.empty()) {
      ACTS_WARNING("Empty track found.");
      continue;
    }
      std::ranges::reverse(surfSequence);

//    std::reverse(trackSourceLinks.begin(), trackSourceLinks.end());
    //std::ranges::reverse(trackSourceLinks);

    ACTS_VERBOSE("Initial parameters: "
                 << initialParams.fourPosition(ctx.geoContext).transpose()
                 << " -> " << initialParams.direction().transpose());

    auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
                               calibrator, surfSequence, tracks);

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