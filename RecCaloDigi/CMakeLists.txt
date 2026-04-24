################################################################################
# Package: RecCaloDigi
################################################################################

gaudi_add_library(RecCaloDigi
  SOURCES
    src/CalorimeterHitType.cc
    src/RealisticCaloDigi.cc
    src/RealisticCaloReco.cc
  LINK
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    DD4hep::DDCore
    DD4hep::DDRec
    ROOT::Core
    ROOT::Hist
    CLHEP::Random
)

target_include_directories(RecCaloDigi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Concrete algorithm plugins — loaded at runtime by Gaudi
gaudi_add_module(k4RecCaloDigiPlugins
  SOURCES
    src/FilterDoubleLayerHits.cc
    src/RealisticCaloDigiScinPpd.cc
    src/RealisticCaloDigiSilicon.cc
    src/RealisticCaloRecoScinPpd.cc
    src/RealisticCaloRecoSilicon.cc
  LINK
    RecCaloDigi
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    DD4hep::DDCore
    DD4hep::DDRec
    ROOT::Core
    ROOT::Hist
    CLHEP::Random
)

target_include_directories(k4RecCaloDigiPlugins PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
)

install(TARGETS RecCaloDigi k4RecCaloDigiPlugins
  EXPORT k4RecCalorimeterTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
                                                COMPONENT dev
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
