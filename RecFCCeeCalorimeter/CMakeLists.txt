################################################################################
# Package: RecFCCeeCalorimeter
################################################################################

file(GLOB _module_sources src/components/*.cpp)
gaudi_add_module(k4RecFCCeeCalorimeterPlugins
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
                      onnxruntime::onnxruntime
                      nlohmann_json::nlohmann_json
                      RecCaloCommon
                      )
install(TARGETS k4RecFCCeeCalorimeterPlugins
  EXPORT k4RecCalorimeterTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)


install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/tests DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/RecFCCeeCalorimeter)

add_test(NAME ALLEGRO_o1_v03_sim_reco
         COMMAND ${PROJECT_SOURCE_DIR}/RecFCCeeCalorimeter/tests/options/ALLEGRO_o1_v03.sh
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/build/Testing/Temporary
)
set_test_env(ALLEGRO_o1_v03_sim_reco)
set_tests_properties(ALLEGRO_o1_v03_sim_reco PROPERTIES FIXTURES_REQUIRED ALLEGRO_sim_files)

add_test(NAME FCCeeLAr_benchmarkCalibration
         COMMAND k4run RecFCCeeCalorimeter/tests/options/fcc_ee_caloBenchmarkCalibration.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_benchmarkCalibration)

add_test(NAME FCCeeLAr_xtalkNeighbours
  COMMAND k4run RecFCCeeCalorimeter/tests/options/runCaloXTalkNeighbours.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_xtalkNeighbours)

add_test(NAME FCCeeLAr_benchmarkCorrection
         COMMAND ${PROJECT_SOURCE_DIR}/RecFCCeeCalorimeter/tests/options/test_benchmarkCorr.sh
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/build/Testing/Temporary
)
set_test_env(FCCeeLAr_benchmarkCorrection)
set_tests_properties(FCCeeLAr_benchmarkCorrection PROPERTIES FIXTURES_REQUIRED ALLEGRO_sim_files)

add_test(NAME IDEA_o1_v03_reco
         COMMAND k4run RecFCCeeCalorimeter/tests/options/IDEA_o1_v03_reco.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(IDEA_o1_v03_reco)
set_tests_properties(IDEA_o1_v03_reco PROPERTIES DEPENDS simulateSiPM)

add_test(NAME FCCeeCalo_DigitizationWhiteningFilter
         COMMAND k4run RecFCCeeCalorimeter/tests/options/runDigitizationAndMatchedFilter.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeCalo_DigitizationWhiteningFilter)
set_tests_properties(ALLEGRO_o1_v03_sim_reco PROPERTIES FIXTURES_REQUIRED ALLEGRO_sim_files)
