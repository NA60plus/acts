// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/SecondaryVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/SingleSeedVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/VertexFitterAlgorithm.hpp"
#include "ActsExamples/Vertexing/TrackletVertexingAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addVertexing(Context& ctx) {
  using Seeder = ActsExamples::AdaptiveMultiVertexFinderAlgorithm::SeedFinder;
  auto mex = ctx.get("examples");
  auto& m = ctx.get("main");

  py::enum_<Seeder>(m, "VertexSeedFinder")
      .value("TruthSeeder", Seeder::TruthSeeder)
      .value("GaussianSeeder", Seeder::GaussianSeeder)
      .value("AdaptiveGridSeeder", Seeder::AdaptiveGridSeeder);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::AdaptiveMultiVertexFinderAlgorithm, mex,
      "AdaptiveMultiVertexFinderAlgorithm", inputTrackParameters,
<<<<<<< HEAD
      inputTruthParticles, inputTruthVertices, outputProtoVertices,
      outputVertices, seedFinder, bField, minWeight, doSmoothing, maxIterations,
      useTime, tracksMaxZinterval, initialVariances, doFullSplitting,
      tracksMaxSignificance, maxMergeVertexSignificance, spatialBinExtent,
      temporalBinExtent);
=======
      outputProtoVertices, outputVertices, seedFinder, bField, useTime,
      spatialBinExtent, temporalBinExtent);
>>>>>>> main

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::IterativeVertexFinderAlgorithm,
                                mex, "IterativeVertexFinderAlgorithm",
                                inputTrackParameters, outputProtoVertices,
                                outputVertices, bField,
                                significanceCutSeeding,
                                maximumChi2cutForSeeding,
                                maxVertices,
                                createSplitVertices,
                                splitVerticesTrkInvFraction,
                                reassignTracksAfterFirstFit,
                                doMaxTracksCut,
                                maxTracks,
                                cutOffTrackWeight,
                                cutOffTrackWeightReassign,
                                rejectedFraction);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SecondaryVertexFinderAlgorithm,
                                mex, "SecondaryVertexFinderAlgorithm",
                                inputTrackParameters, outputProtoVertices,
                                outputVertices, bField, outputMasses);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::VertexFitterAlgorithm, mex,
                                "VertexFitterAlgorithm", inputTrackParameters,
                                inputProtoVertices, outputVertices, bField,
                                doConstrainedFit, constraintPos, constraintCov);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SingleSeedVertexFinderAlgorithm,
                                mex, "SingleSeedVertexFinderAlgorithm",
                                inputSpacepoints, outputVertices);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackletVertexingAlgorithm,
                                mex, "TrackletVertexingAlgorithm",
                                inputSpacePoints, inputSpacePointsMC, inputParticles, inputMeasurementParticlesMap,                            
                                zMaxTop, zMinTop, zMaxBot, zMinBot,deltaPhi, deltaThetaMin, deltaThetaMax, verbose, doMCtruth, outputRecPrimaryVertex, outputFitPrimaryVertex, outputGenPrimaryVertex, noGuessing,
                                useFit, nbins, zPerigee, outputFitFunction, outputZTracklets, outputZTrackletsPeak);
}

}  // namespace Acts::Python
