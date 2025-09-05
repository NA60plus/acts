// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/FilterMeasurements/FilterMeasurementsAlgorithm.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"  // For Acts::SourceLink
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"  // Contains MeasurementContainer definition
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <optional>  // For std::optional (if SourceLink is optional)
#include <set>       // For std::set
#include <stdexcept>
#include <unordered_map>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

ActsExamples::FilterMeasurementsAlgorithm::FilterMeasurementsAlgorithm(
    ActsExamples::FilterMeasurementsAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("FilterMeasurementsAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing tracks input collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if (m_cfg.outputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "Missing particle-to-measurements map output collection");
  }
  if (m_cfg.outputSimHitMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "Missing particle-to-simulated-hits map output collection");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementParticlesMap.initialize(
      m_cfg.outputMeasurementParticlesMap);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputParticleMeasurementsMap.initialize(
      m_cfg.outputParticleMeasurementsMap);
  m_outputSimHitMeasurementsMap.initialize(m_cfg.outputSimHitMeasurementsMap);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputSimHitMeasurementsMap);
  m_inputHits.initialize(m_cfg.inputSimHits);
}

ActsExamples::ProcessCode ActsExamples::FilterMeasurementsAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Get the input tracks and measurements
  const auto& tracks = m_inputTracks(ctx);
  const auto& inputMeasurements = m_inputMeasurements(ctx);
  const auto& inputMeasurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);
  const auto& simHits = m_inputHits(ctx);


  IndexMultimap<SimBarcode> measurementParticlesMap;
  IndexMultimap<Index> measurementSimHitsMap;
  measurementParticlesMap.reserve(simHits.size());
  measurementSimHitsMap.reserve(simHits.size());

  ACTS_VERBOSE("Number of input tracks: " << tracks.size());

  // --- Step 1: Identify measurements to be removed ---
  // Use a set to store the identifiers of measurements to be removed.
  // A pair of (GeometryIdentifier, Index) uniquely identifies a measurement.
  std::set<std::pair<Acts::GeometryIdentifier,
                     ActsExamples::MeasurementContainer::Index>>
      measurementsToRemove;

  for (const auto& track : tracks) {
    for (auto ts : track.trackStatesReversed()) {
      if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        // Access the source link from the track state proxy.
        // In newer Acts versions, `sourceLink()` usually returns an
        // `std::optional<const Acts::SourceLink&>`. You need to check if it has
        // a value.
        //
        // NOTE: The exact way to get `geometryId` and `index` from
        // `Acts::SourceLink` might depend on the specific type of SourceLink
        // (e.g., `IndexSourceLink`). `Acts::SourceLink` itself is a base class.
        // A common pattern is to use `dynamic_cast` or a visitor pattern if you
        // expect different SourceLink types, but for basic measurement
        // filtering, it's often an `IndexSourceLink`.
        //
        // Assuming SourceLink directly provides identifier() and index() or
        // you can cast it to IndexSourceLink or a similar type for access.
        // Let's use a robust way for Acts::SourceLink.
        Acts::SourceLink sourceLink = ts.getUncalibratedSourceLink();
        Acts::GeometryIdentifier geoId = ts.referenceSurface().geometryId();
        const auto& indexSourceLink = sourceLink.get<IndexSourceLink>();
        // Check if the SourceLink is of a type that contains an Index.
        // This usually means it's an Acts::MeasurementSourceLink or
        // Acts::IndexSourceLink. Acts::SourceLink itself doesn't directly
        // expose .index() or .geometryId().
        //
        // You might need a more specific way to extract this, e.g.:
        // if (const ActsExamples::IndexSourceLink* idxSourceLink =
        // dynamic_cast<const ActsExamples::IndexSourceLink*>(&sourceLink)) {
        //     measurementsToRemove.insert({idxSourceLink->geometryId(),
        //     idxSourceLink->index()});
        // } else {
        //     ACTS_WARNING("Track state has a MeasurementFlag but SourceLink is
        //     not an IndexSourceLink.");
        // }
        //
        // However, a common simplification if you *know* it's an
        // IndexSourceLink: Or if your Acts version has a helper that
        // casts/converts it:
        measurementsToRemove.insert({geoId, indexSourceLink.index()});
        // The above line `sourceLink.geometryId(), sourceLink.index()` is the
        // most common pattern for Acts::SourceLink in many use cases. If it
        // still complains about `no member named geometryId`, you might need
        // the `dynamic_cast` approach or check how your `SourceLink` type is
        // constructed/used. For Acts::SourceLink, `identifier()` usually
        // returns `GeometryIdentifier` and `index()` returns `Index`. Let's
        // stick with `identifier()` and `index()` for `Acts::SourceLink` if
        // these are indeed provided by the base `Acts::SourceLink`. If not, you
        // may need to use `static_cast` if you are sure it's an
        // `IndexSourceLink` from `ActsExamples`.
      }
    }
  }

  ACTS_VERBOSE("Identified " << measurementsToRemove.size()
                             << " unique measurements to remove.");

  // --- Step 2: Create a new MeasurementContainer with only the desired
  // measurements ---
  ActsExamples::MeasurementContainer filteredMeasurements;
  // Pre-allocate space if it makes sense for performance
  if (inputMeasurements.size() >= measurementsToRemove.size()) {
    filteredMeasurements.reserve(inputMeasurements.size() -
                                 measurementsToRemove.size());
  }

  for (const ActsExamples::MeasurementContainer::ConstVariableProxy
           currentMeasurement : inputMeasurements) {
    // Check if the current measurement's identifier (geometryId, index) is in
    // our 'toRemove' set.
    auto currentMeasurementId = std::make_pair(currentMeasurement.geometryId(),
                                               currentMeasurement.index());

    if (measurementsToRemove.find(currentMeasurementId) ==
        measurementsToRemove.end()) {
      // This measurement is NOT in the set of measurements to remove, so keep
      // it.
      auto measurement = filteredMeasurements.copyMeasurement(currentMeasurement);

      auto range = inputMeasurementSimHitsMap.equal_range(currentMeasurement.index());
      auto simHitIdx = range.first->second;  

      measurementParticlesMap.emplace_hint(
          measurementParticlesMap.end(), measurement.index(),
          simHits.nth(simHitIdx)->particleId());
      measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                         measurement.index(), simHitIdx);
    }
  }

  ACTS_VERBOSE("Original measurements count: " << inputMeasurements.size());
  ACTS_VERBOSE("Filtered measurements count: " << filteredMeasurements.size());

  // --- Step 3: Write the filtered measurements to the output collection ---
  m_outputMeasurements(ctx, std::move(filteredMeasurements));
  // invert them before they are moved
  m_outputParticleMeasurementsMap(
      ctx, invertIndexMultimap(measurementParticlesMap));
  m_outputSimHitMeasurementsMap(ctx,
                                invertIndexMultimap(measurementSimHitsMap));

  m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));

  return ActsExamples::ProcessCode::SUCCESS;
}