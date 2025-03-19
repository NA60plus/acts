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
  if (m_cfg.inputTrackParametersMS.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackParametersMS) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackParametersMS.emplace_back(
        std::make_unique<ReadDataHandle<TrackParametersContainer>>(
            this, "inputTrackParametersMS#" +
                      std::to_string(m_inputTrackParametersMS.size())));
    handle->initialize(spName);
  }

  if (m_cfg.inputTrackContainerMS.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackContainerMS) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackContainerMS.emplace_back(
        std::make_unique<ReadDataHandle<ConstTrackContainer>>(
            this, "inputTrackContainerMS#" +
                      std::to_string(m_inputTrackContainerMS.size())));
    handle->initialize(spName);
  }

  if (m_cfg.inputTrackParametersVT.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackParametersVT) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackParametersVT.emplace_back(
        std::make_unique<ReadDataHandle<TrackParametersContainer>>(
            this, "inputTrackParametersVT#" +
                      std::to_string(m_inputTrackParametersVT.size())));
    handle->initialize(spName);
  }

  if (m_cfg.inputTrackContainerVT.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackContainerVT) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackContainerVT.emplace_back(
        std::make_unique<ReadDataHandle<ConstTrackContainer>>(
            this, "inputTrackContainerVT#" +
                      std::to_string(m_inputTrackContainerVT.size())));
    handle->initialize(spName);
  }

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
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

  const auto& measurements = m_inputMeasurements(ctx);

  Acts::Vector3 propPoint = Acts::Vector3{m_cfg.px, m_cfg.py, m_cfg.pz};

  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(propPoint);

  //const Acts::Surface* pSurface = m_cfg.trackingGeometry->findSurface(m_cfg.geoIdMatching);

  //auto pSurface = Acts::Surface::makeShared<Acts::Surface>(m_cfg.trackingGeometry->findSurface(m_cfg.geoIdMatching));

  Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> extrapolator(
      Acts::EigenStepper<>(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  // VOID PROPAGATOR
  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.magneticField);

  Acts::PropagatorOptions<Acts::ActionList<Acts::MaterialInteractor>,
                          Acts::AbortList<Acts::EndOfWorldReached>>
      extrapolationOptions(ctx.geoContext, ctx.magFieldContext);
  Acts::PropagatorOptions<> pOptions(ctx.geoContext, ctx.magFieldContext);

  Acts::FullBilloirVertexFitter::Config vertexFitterCfg;
  vertexFitterCfg.extractParameters
      .connect<&Acts::InputTrack::extractParameters>();

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

  auto it1MS = m_inputTrackParametersMS.begin();
  auto it2MS = m_inputTrackContainerMS.begin();
  std::vector<int> ntracksVTs;
  for (; it1MS != m_inputTrackParametersMS.end() &&
         it2MS != m_inputTrackContainerMS.end();
       ++it1MS, ++it2MS) {
    const auto& itrkparMS = *it1MS;
    const auto& itrkconMS = *it2MS;
    // for (const auto& itrkparMS : m_inputTrackParametersMS) {
    const auto& paramsMSinput = (*itrkparMS)(ctx);
    const auto& trackContainterMS = (*itrkconMS)(ctx);

    auto inputTracksMS = makeInputTracks(paramsMSinput);
    indexMS = 0;
    for (auto& trackMS : inputTracksMS) {
      indexVT = -1;
      std::pair<TrackProxyType, TrackProxyType> trackPairTmp(
          trackContainterMS.getTrack(indexMS),
          trackContainterMS.getTrack(indexMS));
      std::pair<int, int> trackPairIndex(indexMS, indexMS);

      identifyContributingParticles(hitParticlesMapMS,
                                    trackContainterMS.getTrack(indexMS),
                                    particleHitCountsMS);


      if (particleHitCountsMS.empty()) {
        ACTS_DEBUG(
            "No truth particle associated with this trajectory with tip index "
            "= "
            << trackContainterMS.getTrack(indexMS).tipIndex());
        continue;
      }
      ActsFatras::Barcode majorityParticleIdMS =
          particleHitCountsMS.front().particleId;
      ActsFatras::Barcode majorityParticleIdVT;

      indexMS += 1;

      Acts::BoundTrackParameters params1 =
          vertexFitterCfg.extractParameters(trackMS);
      const auto resMS =
          extrapolator.propagateToSurface(params1, *pSurface, pOptions);

      if (!resMS.ok()) {
        continue;
      }
      float chi2 = 1e12;
      const auto& endParamsMS = *resMS;

      Acts::BoundVector paramsMS = endParamsMS.parameters();
      const auto covMS = endParamsMS.covariance();

      Acts::ActsScalar covLoc0MS = (*covMS)(0, 0);
      Acts::ActsScalar covLoc1MS = (*covMS)(1, 1);
      Acts::ActsScalar covPhiMS = (*covMS)(2, 2);
      Acts::ActsScalar covThetaMS = (*covMS)(3, 3);
      Acts::ActsScalar covQOvPMS = (*covMS)(4, 4);

      Acts::ActsScalar loc0MS = paramsMS(Acts::BoundIndices::eBoundLoc0);
      Acts::ActsScalar loc1MS = paramsMS(Acts::BoundIndices::eBoundLoc1);
      Acts::ActsScalar phiMS = paramsMS(Acts::BoundIndices::eBoundPhi);
      Acts::ActsScalar thetaMS = paramsMS(Acts::BoundIndices::eBoundTheta);
      Acts::ActsScalar qOvPMS = paramsMS(Acts::BoundIndices::eBoundQOverP);
      
      int indexContVtBest = 0;
      int indexContVt = -1;
      auto it1VT = m_inputTrackParametersVT.begin();
      auto it2VT = m_inputTrackContainerVT.begin();

      for (; it1VT != m_inputTrackParametersVT.end() &&
             it2VT != m_inputTrackContainerVT.end();
           ++it1VT, ++it2VT) {
        const auto& itrkparVT = *it1VT;
        const auto& itrkconVT = *it2VT;

        const auto& paramsVTinput = (*itrkparVT)(ctx);
        const auto& trackContainterVT = (*itrkconVT)(ctx);
        indexContVt++;
        indexVT = -1;
        auto inputTracksVT = makeInputTracks(paramsVTinput);
        int ntracksVT = 0;
        for (auto& trackVT : inputTracksVT) {
          ntracksVT++;
          indexVT++;
          Acts::BoundTrackParameters params2 =
              vertexFitterCfg.extractParameters(trackVT);
          const auto resVT =
              extrapolator.propagateToSurface(params2, *pSurface, pOptions);

          identifyContributingParticles(hitParticlesMapVT,
                                        trackContainterVT.getTrack(indexVT),
                                        particleHitCountsVT);

          if (particleHitCountsVT.empty()) {
            ACTS_DEBUG(
                "No truth particle associated with this trajectory with tip "
                "index = "
                << trackContainterVT.getTrack(indexVT).tipIndex());
            continue;
          }

          if (!resVT.ok()) {
            continue;
          }

          const auto& endParamsVT = *resVT;
          Acts::BoundVector paramsVT = endParamsVT.parameters();
          const auto covVT = endParamsVT.covariance();

          Acts::ActsScalar covLoc0VT = (*covVT)(0, 0);
          Acts::ActsScalar covLoc1VT = (*covVT)(1, 1);
          Acts::ActsScalar covPhiVT = (*covVT)(2, 2);
          Acts::ActsScalar covThetaVT = (*covVT)(3, 3);
          Acts::ActsScalar covQOvPVT = (*covVT)(4, 4);

          Acts::ActsScalar loc0VT = paramsVT(Acts::BoundIndices::eBoundLoc0);
          Acts::ActsScalar loc1VT = paramsVT(Acts::BoundIndices::eBoundLoc1);
          Acts::ActsScalar phiVT = paramsVT(Acts::BoundIndices::eBoundPhi);
          Acts::ActsScalar thetaVT = paramsVT(Acts::BoundIndices::eBoundTheta);
          Acts::ActsScalar qOvPVT = paramsVT(Acts::BoundIndices::eBoundQOverP);

          float chi2tmp  = (loc0VT - loc0MS) * (loc0VT - loc0MS) / (covLoc0VT + covLoc0MS);
                chi2tmp += (loc1VT - loc1MS) * (loc1VT - loc1MS) / (covLoc1VT + covLoc1MS);
                chi2tmp += (phiVT - phiMS) * (phiVT - phiMS) / (covPhiVT + covPhiMS);
                chi2tmp += (thetaVT - thetaMS) * (thetaVT - thetaMS) / (covThetaVT + covThetaMS);
                chi2tmp += (qOvPVT - qOvPMS) * (qOvPVT - qOvPMS) / (covQOvPVT + covQOvPMS);

          if (chi2tmp < chi2) {
            indexContVtBest = indexContVt;
            chi2 = chi2tmp;
            trackPairTmp.second = trackContainterVT.getTrack(indexVT);
            trackPairIndex.second = indexVT;
            majorityParticleIdVT = particleHitCountsVT.front().particleId;
          }

        }
        ntracksVTs.push_back(ntracksVT);
      }

      if (chi2 != m_cfg.chi2max && trackPairTmp.first != trackPairTmp.second) {
        //if (majorityParticleIdVT == majorityParticleIdMS) {
        trackPairs.push_back(trackPairTmp);
          const auto& itrkparVT = m_inputTrackParametersVT[indexContVtBest];
          const auto& paramsVTinput = (*itrkparVT)(ctx);
          
          auto itVTTmp = m_inputTrackContainerVT.begin();
          const auto& itrkconVTTmp = *itVTTmp;
          const auto& trackContainterVT = (*itrkconVTTmp)(ctx);
          if (first) {
            matchedTracksVT.ensureDynamicColumns(trackContainterVT);
            matchedTracksMS.ensureDynamicColumns(trackContainterMS);
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
        
          auto inputTracksVT = makeInputTracks(paramsVTinput);  

          Acts::BoundTrackParameters params2Tmp =
              vertexFitterCfg.extractParameters(inputTracksVT[trackPairIndex.second]);
          const auto resVTki =
              extrapolator.propagateToSurface(params2Tmp, *pSurface, pOptions);
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
  }

  std::vector<Acts::SourceLink> trackSourceLinks;

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  ActsExamples::PassThroughCalibrator pcalibrator;
  ActsExamples::MeasurementCalibratorAdapter calibrator(pcalibrator,
                                                        measurements);

  TrackFitterFunction::GeneralFitterOptions options{
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, pSurface.get(),
      Acts::PropagatorPlainOptions()};

  for (auto& pair : trackPairs) {
    auto trackMS = pair.first;
    auto trackVT = pair.second;

    trackSourceLinks.clear();

    auto paramsFit = trackVT.covariance();
    Eigen::Matrix<double, 6, 6> tempParamsFit = paramsFit;  // Copy the data
    for(int i=0;i<6;i++)
      for(int j=0;j<6;j++)
        tempParamsFit(i, j) *= 100.0;  

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

  int counter = 0;
  for(auto ntracksVT : ntracksVTs){
    if(ntracksVT != 1){
      break;
    }
    counter++;
  }

  auto itVT = m_inputTrackContainerVT.begin();
  //itVT+=counter;
  const auto& itrkconVT = *itVT;
  const auto& trackContainterVT = (*itrkconVT)(ctx);
  ActsExamples::ConstTrackContainer outputTracksVT{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(matchedTracksVT.container())),
      trackContainterVT.trackStateContainerHolder()};

  m_outputTracksVT(ctx, std::move(outputTracksVT));


  auto itMS = m_inputTrackContainerMS.begin();
  const auto& itrkconMS = *itMS;
  const auto& trackContainterMS = (*itrkconMS)(ctx);
  ActsExamples::ConstTrackContainer outputTracksMS{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(matchedTracksMS.container())),
      trackContainterMS.trackStateContainerHolder()};

  m_outputTracksMS(ctx, std::move(outputTracksMS));


  
  return ActsExamples::ProcessCode::SUCCESS;
}