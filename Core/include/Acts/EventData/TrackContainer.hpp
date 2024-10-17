// This file is part of the ACTS project.
//
<<<<<<< HEAD
// Copyright (C) 2016 CERN for the benefit of the ACTS project
=======
// Copyright (C) 2023 CERN for the benefit of the Acts project
>>>>>>> origin/clone_of_main
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainerBackendConcept.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/Utils.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/Iterator.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <any>
#include <cstddef>
#include <iterator>

namespace Acts {

template <typename T>
struct IsReadOnlyTrackContainer;

/// Track container interface class. This type represents a collections of
/// tracks. It uses a backend to store both the actual tracks and the
/// associated track states.
/// @tparam track_container_t the track container backend
/// @tparam traj_t the track state container backend
/// @tparam holder_t ownership management class for the backend
template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t,
          template <typename> class holder_t = detail::RefHolder>
class TrackContainer {
 public:
  /// Indicates if this track container is read-only, or if it can be modified
  static constexpr bool ReadOnly =
      IsReadOnlyTrackContainer<track_container_t>::value;

  /// Indicates if the track state container is read-only, or if it can be
  /// modified
  static constexpr bool TrackStateReadOnly =
      IsReadOnlyMultiTrajectory<traj_t>::value;

  static_assert(ReadOnly == TrackStateReadOnly,
                "Either both track container and track state container need to "
                "be readonly or both have to be readwrite");

  /// The index type of the track container, taken from the container backend
  using IndexType = TrackIndexType;

  /// Sentinel value that indicates an invalid index
  static constexpr IndexType kInvalid = kTrackIndexInvalid;

  /// Alias for the mutable version of a track proxy, with the same backends as
  /// this container
  using TrackProxy =
      Acts::TrackProxy<track_container_t, traj_t, holder_t, false>;

  /// Alias for the const version of a track proxy, with the same backends as
  /// this container
  using ConstTrackProxy =
      Acts::TrackProxy<track_container_t, traj_t, holder_t, true>;

<<<<<<< HEAD
  using TrackContainerBackend = track_container_t;
  using TrackStateContainerBackend = traj_t;

  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;

  using size_type = IndexType;
  using iterator =
      Acts::ContainerIndexIterator<TrackContainer, TrackProxy, false>;
  using const_iterator =
      Acts::ContainerIndexIterator<TrackContainer, ConstTrackProxy, true>;

=======
>>>>>>> origin/clone_of_main
#ifndef DOXYGEN
  friend TrackProxy;
  friend ConstTrackProxy;
#endif

  /// @anchor track_container_construction
  /// @name TrackContainer construction
  ///
  /// Constructors for the track container by using a set of backends
  /// (track + track state). The container can either take ownership of the
  /// backends or just hold references to them. This is driven by the @c
  /// holder_t template parameter, where you can also supply a custom holder
  /// type. Template deduction is used to try to guess the correct holder type.
  ///
  /// @{

  /// Constructor from a track container backend and a track state container
  /// backend
  /// @param container the track container backend
  /// @param traj the track state container backend
  TrackContainer(holder_t<track_container_t> container, holder_t<traj_t> traj)
      : m_container{std::move(container)}, m_traj{std::move(traj)} {}

  /// Constructor from references to a track container backend and to a track
  /// state container backend
  /// @note The track container will not assume ownership over the backends in this case.
  ///       You need to ensure suitable lifetime
  /// @param container the track container backend
  /// @param traj the track state container backend
<<<<<<< HEAD
  TrackContainer(auto& container, auto& traj)
    requires(detail::is_same_template<holder_t, detail::RefHolder>::value)
=======
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<
                detail::is_same_template<H, detail::RefHolder>::value>>
  TrackContainer(track_container_t& container, traj_t& traj)
>>>>>>> origin/clone_of_main
      : m_container{&container}, m_traj{&traj} {}

  /// Constructor from const references to a track container backend and to a
  /// track state container backend
  /// @note The track container will not assume ownership over the backends in this case.
  ///       You need to ensure suitable lifetime
  /// @param container the track container backend
  /// @param traj the track state container backend
<<<<<<< HEAD
  TrackContainer(const auto& container, const auto& traj)
    requires(detail::is_same_template<holder_t,
                                      detail::ConstRefHolder>::value &&
             ReadOnly && TrackStateReadOnly)
=======
  template <
      template <typename> class H = holder_t,
      bool RO = (IsReadOnlyTrackContainer<track_container_t>::value &&
                 IsReadOnlyMultiTrajectory<traj_t>::value),
      typename = std::enable_if_t<
          detail::is_same_template<H, detail::ConstRefHolder>::value && RO>>
  TrackContainer(const track_container_t& container, const traj_t& traj)
>>>>>>> origin/clone_of_main
      : m_container{&container}, m_traj{&traj} {}

  /// @}

  /// @anchor track_container_track_access
  /// @name TrackContainer track (proxy) access and manipulation
  ///
  /// These methods allow accessing tracks, i.e. adding or retrieving a track
  /// proxy that points at a specific track in the container.
  ///
  /// @{

  /// Get a const track proxy for a track index
  /// @param itrack the track index in the container
  /// @return A const track proxy for the index
  ConstTrackProxy getTrack(IndexType itrack) const { return {*this, itrack}; }

  /// Get a mutable track proxy for a track index
  /// @note Only available if the track container is not read-only
  /// @param itrack the track index in the container
  /// @return A mutable track proxy for the index
<<<<<<< HEAD
  TrackProxy getTrack(IndexType itrack)
    requires(!ReadOnly)
  {
    return {*this, itrack};
  }

  /// Get a const track proxy for a track index
  /// @param itrack the track index in the container
  /// @return A const track proxy for the index
  ConstTrackProxy at(IndexType itrack) const { return getTrack(itrack); }

  /// Get a mutable track proxy for a track index
  /// @note Only available if the track container is not read-only
  /// @param itrack the track index in the container
  /// @return A mutable track proxy for the index
  TrackProxy at(IndexType itrack)
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy getTrack(IndexType itrack) {
>>>>>>> origin/clone_of_main
    return {*this, itrack};
  }

  /// Add a track to the container. Note this only creates the logical track and
  /// allocates memory. You can combine this with @c getTrack to obtain a track proxy
  /// @note Only available if the track container is not read-only
  /// @return the index to the newly added track
<<<<<<< HEAD
  IndexType addTrack()
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType addTrack() {
>>>>>>> origin/clone_of_main
    auto track = getTrack(m_container->addTrack_impl());
    track.tipIndex() = kInvalid;
    return track.index();
  }

  /// Add a track to the container and return a track proxy to it
  /// This effectively calls @c addTrack and @c getTrack
  /// @note Only available if the track container is not read-only
  /// @return a track proxy to the newly added track
<<<<<<< HEAD
  TrackProxy makeTrack()
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy makeTrack() {
>>>>>>> origin/clone_of_main
    return getTrack(addTrack());
  }

  /// Remove a track at index @p itrack from the container
  /// @note Only available if the track container is not read-only
  /// @note This invalidates track proxies that point to tracks with larger
  ///       indices than @p itrack!
  /// @param itrack The index of the track to remove
<<<<<<< HEAD
  void removeTrack(IndexType itrack)
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void removeTrack(IndexType itrack) {
>>>>>>> origin/clone_of_main
    m_container->removeTrack_impl(itrack);
  }

  /// Get a mutable iterator to the first track in the container
  /// @note Only available if the track container is not read-only
  /// @return a mutable iterator to the first track
<<<<<<< HEAD
  iterator begin()
    requires(!ReadOnly)
  {
    return iterator{*this, 0};
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto begin() {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         TrackProxy, false>{*this, 0};
>>>>>>> origin/clone_of_main
  }

  /// Get a past-the-end iterator for this container
  /// @note Only available if the track container is not read-only
  /// @return a past-the-end iterator
<<<<<<< HEAD
  iterator end()
    requires(!ReadOnly)
  {
    return iterator{*this, size()};
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto end() {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         TrackProxy, false>{*this, size()};
>>>>>>> origin/clone_of_main
  }

  /// Get an const iterator to the first track in the container
  /// @return a const iterator to the first track
  const_iterator begin() const { return const_iterator{*this, 0}; }

  /// Get a past-the-end iterator for this container
  /// @return a past-the-end iterator
  const_iterator end() const { return const_iterator{*this, size()}; }

  /// @}

  /// @anchor track_container_columns
  /// @name TrackContainer column management
  /// TrackContainer can manage a set of common static columns, and dynamic
  /// columns that can be added at runtime. This set of methods allows you to
  /// manage the dynamic columns.
  /// @{

  /// Add a dymanic column to the track container
  /// @note Only available if the track container is not read-only
  /// @param key the name of the column to be added
<<<<<<< HEAD
  template <typename T>
  constexpr void addColumn(std::string_view key)
    requires(!ReadOnly)
  {
=======
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(std::string_view key) {
>>>>>>> origin/clone_of_main
    m_container->template addColumn_impl<T>(key);
  }

  /// Check if this track container has a specific dynamic column
  /// @param key the key to check for
  constexpr bool hasColumn(const std::string& key) const {
    return m_container->hasColumn_impl(hashString(key));
  }

  /// Check if a this track container has a specific dynamic column
  /// @param key the key to check for
  constexpr bool hasColumn(HashedString key) const {
    return m_container->hasColumn_impl(key);
  }

  /// Helper function to make this track container match the dynamic columns of
  /// another one. This will only work if the track container supports this
  /// source, and depends on the implementation details of the dynamic columns
  /// of the container
  /// @note Only available if the track container is not read-only
  /// @tparam other_track_container_t Type of the other track container
  /// @param other The other track container
<<<<<<< HEAD
  template <typename other_track_container_t>
  void ensureDynamicColumns(const other_track_container_t& other)
    requires(!ReadOnly)
  {
=======
  template <typename other_track_container_t, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void ensureDynamicColumns(const other_track_container_t& other) {
>>>>>>> origin/clone_of_main
    container().ensureDynamicColumns_impl(other.container());
  }

  /// @}

  /// @anchor track_congtainer_backend_access
  /// @name TrackContainer backend access
  /// These methods allow accessing the backend of the track container. In most
  /// cases, this is not necessary for interacting with the container.
  /// @{

  /// Get a mutable reference to the track container backend
  /// @note Only available if the track container is not read-only
  /// @return a mutable reference to the backend
<<<<<<< HEAD
  auto& container()
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& container() {
>>>>>>> origin/clone_of_main
    return *m_container;
  }

  /// Get a const reference to the track container backend
  /// @return a const reference to the backend
  const auto& container() const { return *m_container; }

  /// Get a mutable reference to the track state container backend
  /// @note Only available if the track container is not read-only
  /// @return a mutable reference to the backend
<<<<<<< HEAD
  auto& trackStateContainer()
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& trackStateContainer() {
>>>>>>> origin/clone_of_main
    return *m_traj;
  }

  /// Retrieve the holder of the track state container
  /// @return The track state container including it's holder
  /// @note Only available if the track container is not read-only
<<<<<<< HEAD
  auto& trackStateContainerHolder()
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& trackStateContainerHolder() {
>>>>>>> origin/clone_of_main
    return m_traj;
  }

  /// Get a const reference to the track state container backend
  /// @return a const reference to the backend
  const auto& trackStateContainer() const { return *m_traj; }

  /// Retrieve the holder of the track state container
  /// @return The track state container including it's holder
  const auto& trackStateContainerHolder() const { return m_traj; }

  /// @}

  /// Get the size (number of tracks) of the track container
  /// @return the sixe
  constexpr IndexType size() const { return m_container->size_impl(); }

  /// Clear the content of the track container
  /// @note Only available if the track container is not read-only
<<<<<<< HEAD
  void clear()
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void clear() {
>>>>>>> origin/clone_of_main
    m_container->clear();
    m_traj->clear();
  }

 protected:
<<<<<<< HEAD
  template <typename T, HashedString key>
  constexpr T& component(IndexType itrack)
    requires(!ReadOnly)
  {
    return *std::any_cast<T*>(container().component_impl(key, itrack));
  }

  template <typename T>
  constexpr T& component(HashedString key, IndexType itrack)
    requires(!ReadOnly)
  {
=======
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component(IndexType itrack) {
    return *std::any_cast<T*>(container().component_impl(key, itrack));
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key, IndexType itrack) {
>>>>>>> origin/clone_of_main
    return *std::any_cast<T*>(container().component_impl(key, itrack));
  }

  template <typename T, HashedString key>
  constexpr const T& component(IndexType itrack) const {
    return *std::any_cast<const T*>(container().component_impl(key, itrack));
  }

  template <typename T>
  constexpr const T& component(HashedString key, IndexType itrack) const {
    return *std::any_cast<const T*>(container().component_impl(key, itrack));
  }

<<<<<<< HEAD
  constexpr typename TrackProxy::Parameters parameters(IndexType itrack)
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackProxy::Parameters parameters(IndexType itrack) {
>>>>>>> origin/clone_of_main
    return container().parameters(itrack);
  }

  constexpr typename ConstTrackProxy::ConstParameters parameters(
      IndexType itrack) const {
    return container().parameters(itrack);
  }

<<<<<<< HEAD
  constexpr typename TrackProxy::Covariance covariance(IndexType itrack)
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackProxy::Covariance covariance(IndexType itrack) {
>>>>>>> origin/clone_of_main
    return container().covariance(itrack);
  }

  constexpr typename ConstTrackProxy::ConstCovariance covariance(
      IndexType itrack) const {
    return container().covariance(itrack);
  }

<<<<<<< HEAD
  auto reverseTrackStateRange(IndexType itrack)
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto reverseTrackStateRange(IndexType itrack) {
>>>>>>> origin/clone_of_main
    auto tip = component<IndexType, hashString("tipIndex")>(itrack);
    return m_traj->reverseTrackStateRange(tip);
  }

  auto reverseTrackStateRange(IndexType itrack) const {
    auto tip = component<IndexType, hashString("tipIndex")>(itrack);
    return m_traj->reverseTrackStateRange(tip);
  }

<<<<<<< HEAD
  auto forwardTrackStateRange(IndexType itrack)
    requires(!ReadOnly)
  {
=======
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto forwardTrackStateRange(IndexType itrack) {
>>>>>>> origin/clone_of_main
    auto stem = component<IndexType, hashString("stemIndex")>(itrack);
    if (stem == kInvalid) {
      throw std::invalid_argument{"Track has no stem index"};
    }
    return m_traj->forwardTrackStateRange(stem);
  }

  auto forwardTrackStateRange(IndexType itrack) const {
    auto stem = component<IndexType, hashString("stemIndex")>(itrack);
    if (stem == kInvalid) {
      throw std::invalid_argument{"Track has no stem index"};
    }
    return m_traj->forwardTrackStateRange(stem);
  }

 private:
<<<<<<< HEAD
  template <typename T>
  void copyDynamicFrom(IndexType dstIdx, const T& src, IndexType srcIdx)
    requires(!ReadOnly)
  {
=======
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void copyDynamicFrom(IndexType dstIdx, const T& src, IndexType srcIdx) {
>>>>>>> origin/clone_of_main
    const auto& dynamicKeys = src.dynamicKeys_impl();
    for (const auto key : dynamicKeys) {
      std::any srcPtr = src.component_impl(key, srcIdx);
      container().copyDynamicFrom_impl(dstIdx, key, srcPtr);
    }
  }

  const_if_t<ReadOnly, holder_t<track_container_t>> m_container;
  const_if_t<ReadOnly, holder_t<traj_t>> m_traj;
};

template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t>
TrackContainer(track_container_t& container, traj_t& traj)
    -> TrackContainer<track_container_t, traj_t, detail::RefHolder>;

template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t>
TrackContainer(const track_container_t& container, const traj_t& traj)
    -> TrackContainer<track_container_t, traj_t, detail::ConstRefHolder>;

template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t>
TrackContainer(track_container_t&& container, traj_t&& traj)
    -> TrackContainer<track_container_t, traj_t, detail::ValueHolder>;

}  // namespace Acts
