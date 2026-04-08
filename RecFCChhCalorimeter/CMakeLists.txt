################################################################################
# Package: RecFCChhCalorimeter
################################################################################

file(GLOB _module_sources src/components/*.cpp)
gaudi_add_module(k4RecFCChhCalorimeterPlugins
                 SOURCES ${_module_sources}
                 LINK k4FWCore::k4FWCore
                      k4FWCore::k4Interface
                      Gaudi::GaudiKernel
                      DD4hep::DDCore
                      EDM4HEP::edm4hep
                      k4geo::detectorSegmentations
                      k4geo::detectorCommon
                      RecCaloCommon
                      DD4hep::DDG4
                      ROOT::Core
                      ROOT::Hist
                      )

install(TARGETS k4RecFCChhCalorimeterPlugins
  EXPORT k4RecCalorimeterTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)


install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/tests DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/RecFCChhCalorimeter)

add_test(NAME FCChhBaselineSimDigiReco
         COMMAND ${PROJECT_SOURCE_DIR}/RecFCChhCalorimeter/tests/options/FCChhBaselineSimDigiReco.sh
         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/Testing/Temporary
)
set_test_env(FCChhBaselineSimDigiReco)
set_tests_properties(FCChhBaselineSimDigiReco PROPERTIES FIXTURES_SETUP FCChhBaseline_digi_files)

add_test(NAME FCChhBaselineCellPositions
         COMMAND ${PROJECT_SOURCE_DIR}/RecFCChhCalorimeter/tests/options/FCChhBaselineCellPositions.sh
         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/Testing/Temporary
)
set_test_env(FCChhBaselineCellPositions)
set_tests_properties(FCChhBaselineCellPositions PROPERTIES FIXTURES_REQUIRED FCChhBaseline_digi_files)

#
#gaudi_add_test(buildingCellNoiseMap
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#	       FRAMEWORK tests/options/noiseLevelPerCell.py)
#
#gaudi_add_test(buildingCellNeighboursMap
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#	       FRAMEWORK tests/options/neighbours.py)
