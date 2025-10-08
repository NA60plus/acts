acts_add_library(
    ExamplesMatching
    SHARED
    src/MatchingAlgorithm.cpp

)

target_include_directories(
    ActsExamplesMatching
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesMatching
    PUBLIC Acts::ExamplesFramework Acts::ExamplesTrackFitting Acts::ExamplesVertexing Acts::PluginJson
)
