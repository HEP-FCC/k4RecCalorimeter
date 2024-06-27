################################################################################
# Package: RecCalorimeter
################################################################################


file(GLOB _module_sources src/components/*.cpp )

gaudi_add_module(k4RecCalorimeterPlugins
                 SOURCES ${_module_sources} 
                 LINK k4FWCore::k4FWCore
                      k4FWCore::k4Interface
                      Gaudi::GaudiAlgLib
                      Gaudi::GaudiKernel
                      DD4hep::DDCore
                      DD4hep::DDG4
                      EDM4HEP::edm4hep
                      ROOT::Core
                      ROOT::Hist
                      k4geo::detectorSegmentations
                      k4geo::detectorCommon
                      ${FASTJET_LIBRARIES}
                      )

install(TARGETS k4RecCalorimeterPlugins
  EXPORT k4RecCalorimeterTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)



include_directories(${FASTJET_INCLUDE_DIRS} include/RecCalorimeter)
link_directories(${FASTJET_LIBRARIES})

include(CTest)

add_test(NAME FCC_createJets
  COMMAND k4run RecCalorimeter/tests/options/runJetReco.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

add_test(NAME FCC_createTruthJets
  COMMAND k4run RecCalorimeter/tests/options/runTruthJetReco.py
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)


#install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/tests/options DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/Reconstruction/RecCalorimeter)
#
#gaudi_add_test(genJetClustering
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/genJetClustering.py)
#
#gaudi_add_test(simulateECalSimple
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/geant_fullsim_ecalSimple_singleparticles.py)
#
#gaudi_add_test(EcalReconstructionNoiseFromFile
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               DEPENDS simulateECalSimple
#               FRAMEWORK tests/options/runEcalSimple_ReconstructionSW_noiseFromFile.py
#               PASSREGEX "Segmentation cells  \\(Nphi, Neta, minEta\\): \\[629, 331, 1\\]" )
#
## TODO: causes issues with CI, investigate
##gaudi_add_test(EcalReconstructionCheckNumClusters
##               ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
##               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
##               COMMAND python Reconstruction/RecCalorimeter/tests/scripts/ecalRecoCheckNumClusters.py
##               DEPENDS EcalReconstructionNoiseFromFile)
#
#gaudi_add_test(simulateHCal
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/geant_fullsim_hcal_singleparticles.py)
#
#gaudi_add_test(HcalDigitisationWithoutNoise
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               DEPENDS simulateHCal
#               FRAMEWORK tests/options/runHcalDigitisationWithoutNoise.py)
#
#gaudi_add_test(HcalDigitisationFlatNoise
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               DEPENDS simulateHCal
#               FRAMEWORK tests/options/runHcalDigitisationFlatNoise.py)
#
#gaudi_add_test(simulateFullCaloSystem
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/runFullCaloSystem_SimAndDigitisation.py)
#
#gaudi_add_test(reconstructionFullCaloSystem
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                DEPENDS simulateFullCaloSystem
#                FRAMEWORK tests/options/runFullCaloSystem_ReconstructionSW_noNoise.py)
#
#gaudi_add_test(reconstructionFullCaloSystemWithNoise
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                DEPENDS simulateFullCaloSystem
#                FRAMEWORK tests/options/runFullCaloSystem_ReconstructionSW_noiseFromFile.py
#                PASSREGEX "Segmentation cells  \\(Nphi, Neta, minEta\\): \\[704, 335, 1\\]")
#
#gaudi_add_test(correctionFullCaloSystem
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                DEPENDS reconstructionFullCaloSystem
#                FRAMEWORK tests/options/runFullCaloSystem_CorrectECalBarrelSliWinClusters.py)
#
#gaudi_add_test(estimatePileup
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                FRAMEWORK tests/options/runPreparePileup.py)
#
#gaudi_add_test(recTopoClustersBarrelCaloSystem
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                DEPENDS simulateFullCaloSystem
#                FRAMEWORK tests/options/runBarrelCaloSystem_ReconstructionTopoClusters_noNoise.py)
#
#gaudi_add_test(recTopoClustersBarrelCaloSystemElecNoise
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                DEPENDS simulateFullCaloSystem
#                FRAMEWORK tests/options/runBarrelCaloSystem_ReconstructionTopoClusters_electrNoise.py)
#
#gaudi_add_test(recSplitTopoClustersBarrelCaloSystemElecNoiseConeSelected
#                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                DEPENDS simulateFullCaloSystem
#		FRAMEWORK tests/options/runBarrelCaloSystem_ReconstructionTopoClusters_electrNoise_coneSelection.py)
