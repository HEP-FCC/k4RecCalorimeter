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

add_test(NAME FCCeeLAr_simulateForReco
         COMMAND k4run RecFCCeeCalorimeter/tests/options/runCaloSim.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_simulateForReco)

add_test(NAME FCCeeLAr_simRecoAllegroV3
         COMMAND k4run RecFCCeeCalorimeter/tests/options/run_thetamodulemerged.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_simRecoAllegroV3)

add_test(NAME FCCeeLAr_slidingWindowClustering
         COMMAND k4run RecFCCeeCalorimeter/tests/options/runFullCaloSystem_ReconstructionSW_noiseFromFile.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_slidingWindowClustering)
set_tests_properties(FCCeeLAr_slidingWindowClustering PROPERTIES DEPENDS "FCCeeLAr_simulateForReco")

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

add_test(NAME FCCeeLAr_ecalNoise
  COMMAND k4run RecCalorimeter/tests/options/runEcalBarrel_addNoise_thetamodulemerged.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
set_test_env(FCCeeLAr_ecalNoise)
