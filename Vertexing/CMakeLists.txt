# Displaced vertexing with muon-system tracks only (LLP decays behind the solenoid).
# Author: Mahmoud Althakeel

gaudi_add_module(Vertexing
  SOURCES
    src/MuonArcVertexing.cpp
    components/MuonSystemArcVertexFinder.cpp
  LINK
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    ROOT::Core
    ROOT::MathCore
    ROOT::Matrix
    ROOT::Physics
)

target_include_directories(Vertexing PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install(TARGETS Vertexing
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
)
