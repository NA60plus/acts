add_library(
<<<<<<< HEAD
    ActsExamplesTrackFindingExaTrkX
    SHARED
    src/TrackFindingAlgorithmExaTrkX.cpp
    src/PrototracksToParameters.cpp
    src/TrackFindingFromPrototrackAlgorithm.cpp
    src/TruthGraphBuilder.cpp
=======
  ActsExamplesTrackFindingExaTrkX SHARED
  src/TrackFindingAlgorithmExaTrkX.cpp
  src/PrototracksToParameters.cpp
  src/TrackFindingFromPrototrackAlgorithm.cpp
>>>>>>> origin/clone_of_main
)

target_include_directories(
    ActsExamplesTrackFindingExaTrkX
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTrackFindingExaTrkX
    PUBLIC ActsPluginExaTrkX ActsExamplesFramework ActsExamplesTrackFinding
)

install(
    TARGETS ActsExamplesTrackFindingExaTrkX
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
