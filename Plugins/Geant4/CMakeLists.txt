include(FetchContent)

add_library(
<<<<<<< HEAD
    ActsPluginGeant4
    SHARED
    src/Geant4Converters.cpp
    src/Geant4DetectorElement.cpp
    src/Geant4DetectorSurfaceFactory.cpp
    src/Geant4PhysicalVolumeSelectors.cpp
)
=======
  ActsPluginGeant4 SHARED
  src/Geant4Converters.cpp
  src/Geant4DetectorElement.cpp
  src/Geant4DetectorSurfaceFactory.cpp)
>>>>>>> main
target_include_directories(
    ActsPluginGeant4
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_include_directories(
<<<<<<< HEAD
=======
      ActsPluginGeant4
      SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})    
target_link_libraries(
>>>>>>> main
    ActsPluginGeant4
    SYSTEM
    PRIVATE ${Geant4_INCLUDE_DIRS}
)
target_link_libraries(ActsPluginGeant4 PUBLIC ActsCore ${Geant4_LIBRARIES})

install(
    TARGETS ActsPluginGeant4
    EXPORT ActsPluginGeant4Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
