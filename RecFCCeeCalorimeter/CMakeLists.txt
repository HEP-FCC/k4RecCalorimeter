################################################################################
# Package: RecFCCeeCalorimeter
################################################################################

file(GLOB _module_sources src/components/*.cpp)
gaudi_add_module(k4RecFCCeeCalorimeterPlugins
                 SOURCES ${_module_sources}
                 LINK k4FWCore::k4FWCorePlugins GaudiAlgLib  GaudiKernel DD4hep::DDCore EDM4HEP::edm4hep  k4FWCore::k4Interface FCCDetectors::DetSegmentation DD4hep::DDG4 ROOT::Core ROOT::Hist)

install(TARGETS k4RecFCCeeCalorimeterPlugins
  EXPORT k4RecCalorimeterTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)


install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/tests DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/RecFCCeeCalorimeter)

#gaudi_add_test(reconstructFullCaloSystemCellPositions
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#           DEPENDS simulateFullCaloSystemForCellPositions
#           FRAMEWORK tests/options/recoPositions_fullCaloSystem.py)
