// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsExamples/FilterMeasurements/FilterMeasurementsAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {

void addFilterMeasurements(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::FilterMeasurementsAlgorithm, mex,
      "FilterMeasurementsAlgorithm", inputTracks, inputMeasurements,
      outputMeasurements, inputSimHits, outputMeasurementParticlesMap,
      outputMeasurementSimHitsMap, outputParticleMeasurementsMap,
      outputSimHitMeasurementsMap, inputSimHitMeasurementsMap);
}
}  // namespace Acts::Python
