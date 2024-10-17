add_library(
<<<<<<< HEAD
    ActsExamplesVertexing
    SHARED
    src/AdaptiveMultiVertexFinderAlgorithm.cpp
    src/IterativeVertexFinderAlgorithm.cpp
    src/VertexFitterAlgorithm.cpp
    src/SingleSeedVertexFinderAlgorithm.cpp
    src/TruthVertexSeeder.cpp
=======
  ActsExamplesVertexing SHARED
  src/AdaptiveMultiVertexFinderAlgorithm.cpp
  src/IterativeVertexFinderAlgorithm.cpp
  src/SecondaryVertexFinderAlgorithm.cpp
  src/VertexFitterAlgorithm.cpp
  src/SingleSeedVertexFinderAlgorithm.cpp
  src/TrackletVertexingAlgorithm.cpp
  #src/VertexingHelpers.hpp
>>>>>>> origin/clone_of_main
)
target_include_directories(
    ActsExamplesVertexing
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesVertexing
    PUBLIC ActsCore ActsExamplesFramework
    PRIVATE ActsExamplesTruthTracking
)

install(
    TARGETS ActsExamplesVertexing
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
