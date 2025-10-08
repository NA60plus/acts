acts_add_library(
    ExamplesFilterMeasurements
    SHARED
    src/FilterMeasurementsAlgorithm.cpp

)

target_include_directories(
    ActsExamplesFilterMeasurements
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesFilterMeasurements
    PUBLIC Acts::ExamplesFramework Acts::PluginJson
)
