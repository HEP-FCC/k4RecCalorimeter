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
                      DD4hep::DDG4
                      ROOT::Core
                      ROOT::Hist
                      onnxruntime::onnxruntime
                      nlohmann_json::nlohmann_json
                      )
install(TARGETS k4RecFCCeeCalorimeterPlugins
  EXPORT k4RecCalorimeterTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)


install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/tests DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/RecFCCeeCalorimeter)

add_test(NAME ALLEGRO_o1_v03_sim_reco
         COMMAND ../../../RecFCCeeCalorimeter/tests/options/ALLEGRO_o1_v03.sh
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/build/Testing/Temporary
)

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
         COMMAND k4run RecFCCeeCalorimeter/tests/options/run_thetamodulemerged_SW_benchmarkCorr.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_benchmarkCorrection)

add_test(NAME FCCeeLAr_runxtalk
  COMMAND k4run RecCalorimeter/tests/options/runEcalBarrel_ReconstructionTopoClusters_crosstalk.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_runxtalk)
