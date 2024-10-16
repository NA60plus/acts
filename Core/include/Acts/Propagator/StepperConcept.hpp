// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
namespace Concepts {

/// @brief Concept that is satisfied by both single- and multi-steppers.
template <typename Stepper, typename State = typename Stepper::State>
concept CommonStepper = requires {
  typename Stepper::State;
  typename Stepper::Jacobian;
  typename Stepper::Covariance;
  typename Stepper::BoundState;
  typename Stepper::CurvilinearState;

  requires requires(const Stepper& s, State& t) {
    { s.transportCovarianceToCurvilinear(t) } -> std::same_as<void>;

    requires requires(const BoundVector& bv, const BoundSquareMatrix& bm,
                      const Surface& sf, const double d) {
      { s.resetState(t, bv, bm, sf, d) } -> std::same_as<void>;
    };

    requires requires(const Surface& sf, bool b,
                      const FreeToBoundCorrection& corr) {
      {
        s.boundState(t, sf, b, corr)
      } -> std::same_as<Result<typename Stepper::BoundState>>;
      { s.transportCovarianceToBound(t, sf, corr) } -> std::same_as<void>;
    };

    requires requires(bool b) {
      {
        s.curvilinearState(t, b)
      } -> std::same_as<typename Stepper::CurvilinearState>;
    };

    requires requires(const Surface& sf, std::uint8_t ui, Direction d,
                      const BoundaryTolerance& bt, ActsScalar sc,
                      const Logger& l) {
      { s.updateSurfaceStatus(t, sf, ui, d, bt, sc, l) };
    };

// clang-format off
template <typename S>
constexpr bool MultiStepperStateConcept= require<
  has_member<S, cov_transport_t, bool>,
  has_member<S, path_accumulated_t, double>
>;
// clang-format on

// clang-format off
    template <typename S, typename state = typename S::State>
      struct CommonStepperConcept {
        constexpr static bool state_exists = exists<state_t, S>;
        static_assert(state_exists, "State type not found");
        constexpr static bool jacobian_exists = exists<jacobian_t, S>;
        static_assert(jacobian_exists, "Jacobian type not found");
        constexpr static bool covariance_exists = exists<covariance_t, S>;
        static_assert(covariance_exists, "Covariance type not found");
        constexpr static bool bound_state_exists = exists<bound_state_t, S>;
        static_assert(bound_state_exists, "BoundState type not found");
        constexpr static bool curvilinear_state_exists = exists<curvilinear_state_t, S>;
        static_assert(curvilinear_state_exists, "CurvilinearState type not found");
        constexpr static bool reset_state_exists = has_method<const S, void, reset_state_t, state&, const BoundVector&, const BoundSquareMatrix&, const Surface&, const double>;
        static_assert(reset_state_exists, "resetState method not found");
        constexpr static bool position_exists = has_method<const S, Vector3, position_t, const state&>;
        static_assert(position_exists, "position method not found");
        constexpr static bool direction_exists = has_method<const S, Vector3, direction_t, const state&>;
        static_assert(direction_exists, "direction method not found");
        constexpr static bool qop_exists = has_method<const S, double, qop_t, const state&>;
        static_assert(qop_exists, "qOverP method not found");
        constexpr static bool absolute_momentum_exists = has_method<const S, double, absolute_momentum_t, const state&>;
        static_assert(absolute_momentum_exists, "absoluteMomentum method not found");
        constexpr static bool momentum_exists = has_method<const S, Vector3, momentum_t, const state&>;
        static_assert(momentum_exists, "momentum method not found");
        constexpr static bool charge_exists = has_method<const S, double, charge_t, const state&>;
        static_assert(charge_exists, "charge method not found");
        constexpr static bool time_exists = has_method<const S, double, time_t, const state&>;
        static_assert(time_exists, "time method not found");
        constexpr static bool bound_state_method_exists= has_method<const S, Result<typename S::BoundState>, bound_state_method_t, state&, const Surface&, bool, const FreeToBoundCorrection&>;
        static_assert(bound_state_method_exists, "boundState method not found");
        constexpr static bool curvilinear_state_method_exists = has_method<const S, typename S::CurvilinearState, curvilinear_state_method_t, state&, bool>;
        static_assert(curvilinear_state_method_exists, "curvilinearState method not found");
        constexpr static bool covariance_transport_exists = require<has_method<const S, void, covariance_transport_curvilinear_t, state&>,
                                                                    has_method<const S, void, covariance_transport_bound_t, state&, const Surface&, const FreeToBoundCorrection&>>;
        static_assert(covariance_transport_exists, "covarianceTransport method not found");
        constexpr static bool update_surface_exists = has_method<const S, Intersection3D::Status, update_surface_status_t, state&, const Surface&, std::uint8_t, Direction, const BoundaryCheck&, ActsScalar, const Logger&>;
        static_assert(update_surface_exists, "updateSurfaceStatus method not found");
        constexpr static bool update_step_size_exists = has_method<const S, void, update_step_size_t, state&, double, ConstrainedStep::Type, bool>;
        static_assert(update_step_size_exists, "updateStepSize method not found");
        constexpr static bool get_step_size_exists = has_method<const S, double, get_step_size_t, const state &, ConstrainedStep::Type>;
        static_assert(get_step_size_exists, "getStepSize method not found");
        constexpr static bool release_step_size_exists = has_method<const S, void, release_step_size_t, state&, ConstrainedStep::Type>;
        static_assert(release_step_size_exists, "releaseStepSize method not found");
        constexpr static bool output_step_size_exists = has_method<const S, std::string, output_step_size_t, const state&>;
        static_assert(output_step_size_exists, "outputStepSize method not found");

        constexpr static bool value = require<state_exists,
                                              jacobian_exists,
                                              covariance_exists,
                                              bound_state_exists,
                                              curvilinear_state_exists,
                                              position_exists,
                                              direction_exists,
                                              qop_exists,
                                              absolute_momentum_exists,
                                              momentum_exists,
                                              charge_exists,
                                              time_exists,
                                              bound_state_method_exists,
                                              curvilinear_state_method_exists,
                                              covariance_transport_exists,
                                              update_surface_exists,
                                              update_step_size_exists,
                                              release_step_size_exists,
                                              output_step_size_exists>;

      requires requires(double d, bool b) {
        { s.updateStepSize(t, d, st, b) } -> std::same_as<void>;
      };
    };
  };

  requires requires(const Stepper& s, const State& t) {
    { s.position(t) } -> std::same_as<Vector3>;
    { s.direction(t) } -> std::same_as<Vector3>;
    { s.qOverP(t) } -> std::same_as<double>;
    { s.absoluteMomentum(t) } -> std::same_as<double>;
    { s.momentum(t) } -> std::same_as<Vector3>;
    { s.charge(t) } -> std::same_as<double>;
    { s.time(t) } -> std::same_as<double>;
    { s.outputStepSize(t) } -> std::same_as<std::string>;

    requires requires(const ConstrainedStep::Type st) {
      { s.getStepSize(t, st) } -> std::same_as<double>;
    };
  };
};

/// @brief Concept that is satisfied by single-steppers.
template <typename Stepper, typename State = typename Stepper::State>
concept SingleStepper =
    CommonStepper<Stepper, State> && requires(const Stepper& s, State& t) {
      requires requires(const FreeVector& fv, const BoundVector& bv,
                        const BoundSquareMatrix& bm, const Surface& sf) {
        { s.update(t, fv, bv, bm, sf) } -> std::same_as<void>;
      };

      requires requires(const Vector3& v1, const Vector3& v2, double d1,
                        double d2) {
        { s.update(t, v1, v2, d1, d2) } -> std::same_as<void>;
        { s.getField(t, v1) } -> std::same_as<Result<Vector3>>;
      };
    };

/// @brief Concept that is satisfied by multi-steppers.
template <typename Stepper, typename State = typename Stepper::State>
concept MultiStepper = CommonStepper<Stepper, State> && requires {
  // TODO for now we do not check if the ComponentProxy does fulfill a concept
  typename Stepper::ComponentProxy;

  // TODO for now we do not check if the ConstComponentProxy does fulfill a
  // concept
  typename Stepper::ConstComponentProxy;

  requires requires(const Stepper& s, State& t) {
    { s.numberComponents(t) } -> std::same_as<std::size_t>;
    { s.clearComponents(t) } -> std::same_as<void>;
    { s.removeMissedComponents(t) } -> std::same_as<void>;
  };
};
}  // namespace Concepts

/// @brief Concept that is satisfied by steppers.
template <typename _Stepper, typename State = typename _Stepper::State>
concept StepperConcept = Concepts::SingleStepper<_Stepper, State> ||
                         Concepts::MultiStepper<_Stepper, State>;

/// @brief Concept that is satisfied by stepper states.
template <typename State>
concept StepperStateConcept = requires(const State& t) {
  { t.covTransport } -> Concepts::decayed_same_as<const bool&>;
  { t.pathAccumulated } -> Concepts::decayed_same_as<const double&>;
};

}  // namespace Acts
