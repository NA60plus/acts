// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

namespace Acts {

/// This is a fully guided navigator that progresses through a pre-given
/// sequence of surfaces.
///
/// This can either be used as a validation tool, for truth tracking, or track
/// refitting.
class DirectNavigator {
 public:
  /// The sequentially crossed surfaces
  using SurfaceSequence = std::vector<const Surface*>;
  using SurfaceIter = std::vector<const Surface*>::iterator;

  DirectNavigator(std::unique_ptr<const Logger> _logger =
                      getDefaultLogger("DirectNavigator", Logging::INFO))
      : m_logger{std::move(_logger)} {}

  /// @brief Nested Actor struct, called Initializer
  ///
  /// This is needed for the initialization of the surface sequence.
  struct Initializer {
    /// The Surface sequence
    SurfaceSequence navSurfaces = {};

    /// Actor result / state
    struct this_result {
      bool initialized = false;
    };

    using result_type = this_result;

    /// Defaulting the constructor
    Initializer() = default;

    /// Actor operator call
    /// @tparam statet Type of the full propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state the entire propagator state
    /// @param r the result of this Actor
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void operator()(propagator_state_t& state, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, result_type& r,
                    const Logger& /*logger*/) const {
      // Only act once
      if (!r.initialized) {
        // Initialize the surface sequence
        state.navigation.navSurfaces = navSurfaces;
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();

        // In case the start surface is in the list of nav surfaces
        // we need to correct the iterator to point to the next surface
        // in the vector
        if (state.navigation.startSurface) {
          auto surfaceIter = std::find(state.navigation.navSurfaces.begin(),
                                       state.navigation.navSurfaces.end(),
                                       state.navigation.startSurface);
          // if current surface in the list, point to the next surface
          if (surfaceIter != state.navigation.navSurfaces.end()) {
            state.navigation.navSurfaceIter = ++surfaceIter;
          }
        }

        r.initialized = true;
      }
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every
  /// propagation/extrapolation step and keep thread-local navigation
  /// information
  struct State {
    /// Externally provided surfaces - expected to be ordered along the path
    SurfaceSequence navSurfaces = {};

<<<<<<< HEAD
    /// Index of the next surface to try
    int surfaceIndex = 0;
=======
    /// Iterator the next surface
    SurfaceIter navSurfaceIter = navSurfaces.begin();
>>>>>>> main

    /// Navigation state - external interface: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state - external interface: the target surface
    const Surface* targetSurface = nullptr;
    /// Navigation state - starting layer
    const Layer* startLayer = nullptr;
    /// Navigation state - target layer
    const Layer* targetLayer = nullptr;
    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the target volume
    const TrackingVolume* targetVolume = nullptr;

    /// Navigation state - external interface: target is reached
    bool targetReached = false;
    /// Navigation state - external interface: a break has been detected
    bool navigationBreak = false;

    const Surface* navSurface() const {
      return options.surfaces.at(surfaceIndex);
    }

    void nextSurface(Direction direction) {
      if (direction == Direction::Forward) {
        ++surfaceIndex;
      } else {
        --surfaceIndex;
      }
    }

    bool endOfSurfaces() const {
      return surfaceIndex < 0 ||
             surfaceIndex >= static_cast<int>(options.surfaces.size());
    }

    int remainingSurfaces(Direction direction) const {
      if (direction == Direction::Forward) {
        return options.surfaces.size() - surfaceIndex;
      }
      return surfaceIndex + 1;
    }

    void resetSurfaceIndex(Direction direction) {
      surfaceIndex =
          direction == Direction::Forward ? 0 : options.surfaces.size() - 1;
    }
  };

  State makeState(const Surface* startSurface,
                  const Surface* targetSurface) const {
    State result;
    result.startSurface = startSurface;
    result.targetSurface = targetSurface;
    return result;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  bool targetReached(const State& state) const { return state.targetReached; }

  bool endOfWorldReached(State& state) const {
    return state.currentVolume == nullptr;
  }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  void currentSurface(State& state, const Surface* surface) const {
    state.currentSurface = surface;
  }

  void targetReached(State& state, bool targetReached) const {
    state.targetReached = targetReached;
  }

  void navigationBreak(State& state, bool navigationBreak) const {
    state.navigationBreak = navigationBreak;
  }

  /// @brief Initialize call - start of propagation
  ///
  /// @tparam propagator_state_t The state type of the propagator
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state,
                  const stepper_t& /*stepper*/) const {
<<<<<<< HEAD
    ACTS_VERBOSE("Initialize. Surface sequence for navigation:");
    for (const Surface* surface : state.navigation.options.surfaces) {
      ACTS_VERBOSE(surface->geometryId()
                   << " - " << surface->center(state.geoContext).transpose());
    }

    // We set the current surface to the start surface
    state.navigation.currentSurface = state.navigation.options.startSurface;
    if (state.navigation.currentSurface != nullptr) {
      ACTS_VERBOSE("Current surface set to start surface "
=======
    ACTS_VERBOSE(volInfo(state) << "initialize");

    // We set the current surface to the start surface
    state.navigation.currentSurface = state.navigation.startSurface;
    if (state.navigation.currentSurface) {
      ACTS_VERBOSE(volInfo(state)
                   << "Current surface set to start surface "
>>>>>>> main
                   << state.navigation.currentSurface->geometryId());
    } else {
      ACTS_VERBOSE("Current surface set to nullptr");
    }

    // Reset the surface index
    state.navigation.resetSurfaceIndex(state.options.direction);
    for (const Surface* surface : state.navigation.options.surfaces) {
      // make sure we skip over the start surface
      state.navigation.nextSurface(state.options.direction);
      if (surface == state.navigation.currentSurface) {
        break;
      }
    }
    ACTS_VERBOSE("Start surface index set to "
                 << state.navigation.surfaceIndex);
    if (state.navigation.endOfSurfaces()) {
      ACTS_DEBUG(
          "Did not find the start surface in the sequence. Assuming it is not "
          "part of the sequence. Trusting the correctness of the input "
          "sequence. Resetting the surface index.");
      state.navigation.resetSurfaceIndex(state.options.direction);
    }

    state.navigation.navigationBreak = false;
    state.navigation.targetReached = false;
  }

  /// @brief Navigator pre step call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& state, const stepper_t& stepper) const {
<<<<<<< HEAD
    if (state.navigation.navigationBreak) {
      return;
    }

    ACTS_VERBOSE("pre step");
=======
    ACTS_VERBOSE(volInfo(state) << "pre step");
>>>>>>> main

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
<<<<<<< HEAD
    ACTS_VERBOSE(state.navigation.remainingSurfaces(state.options.direction)
                 << " out of " << state.navigation.options.surfaces.size()
                 << " surfaces remain to try.");

    if (state.navigation.endOfSurfaces()) {
=======
    ACTS_VERBOSE(std::distance(state.navigation.navSurfaceIter,
                               state.navigation.navSurfaces.end())
                 << " out of " << state.navigation.navSurfaces.size()
                 << " surfaces remain to try.");

    if (state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      // Establish & update the surface status
      // TODO we do not know the intersection index - passing the closer one
      const auto& surface = **state.navigation.navSurfaceIter;
      const double farLimit = std::numeric_limits<double>::max();
      const auto index =
          chooseIntersection(
              state.geoContext, surface, stepper.position(state.stepping),
              state.options.direction * stepper.direction(state.stepping),
              BoundaryCheck(false), m_nearLimit, farLimit,
              state.options.surfaceTolerance)
              .index();
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, index, state.options.direction,
          BoundaryCheck(false), state.options.surfaceTolerance, *m_logger);
      if (surfaceStatus == Intersection3D::Status::unreachable) {
        ACTS_VERBOSE(
            "Surface not reachable anymore, switching to next one in "
            "sequence");
        // Move the sequence to the next surface
        ++state.navigation.navSurfaceIter;
      } else {
        ACTS_VERBOSE("Navigation stepSize set to "
                     << stepper.outputStepSize(state.stepping));
      }
    } else {
>>>>>>> main
      // Set the navigation break
      ACTS_VERBOSE("End of surfaces reached, navigation break.");
      state.navigation.navigationBreak = true;
      stepper.releaseStepSize(state.stepping, ConstrainedStep::actor);
      // If no externally provided target is given, the target is reached
      if (state.navigation.targetSurface == nullptr) {
        state.navigation.targetReached = true;
        // Announce it then
        ACTS_VERBOSE("No target Surface, job done.");
      }
<<<<<<< HEAD
      return;
    }

    // Establish & update the surface status
    // TODO we do not know the intersection index - passing the closer one
    const auto& surface = *state.navigation.navSurface();
    const double farLimit = std::numeric_limits<double>::max();
    const auto index =
        chooseIntersection(
            state.geoContext, surface, stepper.position(state.stepping),
            state.options.direction * stepper.direction(state.stepping),
            BoundaryTolerance::Infinite(), state.navigation.options.nearLimit,
            farLimit, state.options.surfaceTolerance)
            .index();
    auto surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, surface, index, state.options.direction,
        BoundaryTolerance::Infinite(), state.options.surfaceTolerance,
        *m_logger);
    if (surfaceStatus == Intersection3D::Status::unreachable) {
      ACTS_VERBOSE(
          "Surface not reachable anymore, switching to next one in "
          "sequence");
      // Move the sequence to the next surface
      state.navigation.nextSurface(state.options.direction);
    } else {
      ACTS_VERBOSE("Navigation stepSize set to "
                   << stepper.outputStepSize(state.stepping));
=======
>>>>>>> main
    }
  }

  /// @brief Navigator post step call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void postStep(propagator_state_t& state, const stepper_t& stepper) const {
<<<<<<< HEAD
    if (state.navigation.navigationBreak) {
      return;
    }

    ACTS_VERBOSE("post step");
=======
    ACTS_VERBOSE(volInfo(state) << "post step");
>>>>>>> main

    // Navigator post step always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
<<<<<<< HEAD
    ACTS_VERBOSE(state.navigation.remainingSurfaces(state.options.direction)
                 << " out of " << state.navigation.options.surfaces.size()
                 << " surfaces remain to try.");

    if (state.navigation.endOfSurfaces()) {
      return;
    }

    // Establish the surface status
    // TODO we do not know the intersection index - passing the closer one
    const auto& surface = *state.navigation.navSurface();
    const double farLimit = std::numeric_limits<double>::max();
    const auto index =
        chooseIntersection(
            state.geoContext, surface, stepper.position(state.stepping),
            state.options.direction * stepper.direction(state.stepping),
            BoundaryTolerance::Infinite(), state.navigation.options.nearLimit,
            farLimit, state.options.surfaceTolerance)
            .index();
    auto surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, surface, index, state.options.direction,
        BoundaryTolerance::Infinite(), state.options.surfaceTolerance,
        *m_logger);
    if (surfaceStatus == Intersection3D::Status::onSurface) {
      // Set the current surface
      state.navigation.currentSurface = state.navigation.navSurface();
      ACTS_VERBOSE("Current surface set to  "
                   << state.navigation.currentSurface->geometryId());
      // Move the sequence to the next surface
      state.navigation.nextSurface(state.options.direction);
      if (!state.navigation.endOfSurfaces()) {
        ACTS_VERBOSE("Next surface candidate is  "
                     << state.navigation.options.surfaces
                            .at(state.navigation.surfaceIndex)
                            ->geometryId());
=======
    ACTS_VERBOSE(std::distance(state.navigation.navSurfaceIter,
                               state.navigation.navSurfaces.end())
                 << " out of " << state.navigation.navSurfaces.size()
                 << " surfaces remain to try.");

    // Check if we are on surface
    if (state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      // Establish the surface status
      // TODO we do not know the intersection index - passing the closer one
      const auto& surface = **state.navigation.navSurfaceIter;
      const double farLimit = std::numeric_limits<double>::max();
      const auto index =
          chooseIntersection(
              state.geoContext, surface, stepper.position(state.stepping),
              state.options.direction * stepper.direction(state.stepping),
              BoundaryCheck(false), m_nearLimit, farLimit,
              state.options.surfaceTolerance)
              .index();
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, index, state.options.direction,
          BoundaryCheck(false), state.options.surfaceTolerance, *m_logger);
      if (surfaceStatus == Intersection3D::Status::onSurface) {
        // Set the current surface
        state.navigation.currentSurface = *state.navigation.navSurfaceIter;
        ACTS_VERBOSE("Current surface set to  "
                     << state.navigation.currentSurface->geometryId());
        // Move the sequence to the next surface
        ++state.navigation.navSurfaceIter;
        if (state.navigation.navSurfaceIter !=
            state.navigation.navSurfaces.end()) {
          ACTS_VERBOSE("Next surface candidate is  "
                       << (*state.navigation.navSurfaceIter)->geometryId());
        }
      } else if (surfaceStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE("Next surface reachable at distance  "
                     << stepper.outputStepSize(state.stepping));
>>>>>>> main
      }
    }
  }

 private:
  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.navigation.currentVolume != nullptr
                ? state.navigation.currentVolume->volumeName()
                : "No Volume") +
           " | ";
  }

  ObjectIntersection<Surface> chooseIntersection(
      const GeometryContext& gctx, const Surface& surface,
      const Vector3& position, const Vector3& direction,
      const BoundaryCheck& bcheck, double nearLimit, double farLimit,
      double tolerance) const {
    auto intersections =
        surface.intersect(gctx, position, direction, bcheck, tolerance);

    for (auto& intersection : intersections.split()) {
      if (detail::checkIntersection(intersection, nearLimit, farLimit,
                                    logger())) {
        return intersection;
      }
    }

    return ObjectIntersection<Surface>::invalid();
  }

  const Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Logger> m_logger;

  // TODO https://github.com/acts-project/acts/issues/2738
  /// Distance limit to discard intersections "behind us"
  /// @note this is only necessary because some surfaces have more than one
  ///       intersection
  double m_nearLimit = -100 * UnitConstants::um;
};

}  // namespace Acts
