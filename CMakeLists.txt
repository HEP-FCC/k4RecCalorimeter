CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

project(k4RecCalorimeter)

include(GNUInstallDirs)
include(CTest)

# Set default for install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/InstallArea/ CACHE PATH
    "Install path prefix, prepended onto install directories." FORCE )
endif()

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
if(NOT CMAKE_CXX_STANDARD MATCHES "17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})


#---------------------------------------------------------------
find_package(ROOT COMPONENTS RIO Tree)
find_package(Gaudi)
find_package(k4FWCore) # implicit: Gaudi
#find_package(EDM4HEP 0.10.1) # implicit: Podio
find_package(DD4hep)
find_package(k4geo)
find_package(ONNXRuntime)

# need to use our own FindFastJet.cmake
#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
find_package(FastJet)
#---------------------------------------------------------------


#--- The genConf directory has been renamed to genConfDir in Gaudi 35r1
#--- See https://gitlab.cern.ch/gaudi/Gaudi/-/merge_requests/1158
set(GAUDI_GENCONF_DIR "genConfDir")
if (${Gaudi_VERSION} VERSION_LESS 35.1)
  set(GAUDI_GENCONF_DIR "genConf")
endif()

function(set_test_env _testname)
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$<TARGET_FILE_DIR:k4RecCalorimeterPlugins>:${CMAKE_BINARY_DIR}:$<TARGET_FILE_DIR:k4RecFCCeeCalorimeterPlugins>:${CMAKE_BINARY_DIR}:$<TARGET_FILE_DIR:k4RecFCChhCalorimeterPlugins>:$<TARGET_FILE_DIR:ROOT::Core>:$<TARGET_FILE_DIR:k4FWCore::k4FWCore>:$<TARGET_FILE_DIR:EDM4hep::edm4hep>:$<TARGET_FILE_DIR:podio::podio>:$ENV{LD_LIBRARY_PATH})

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               PYTHONPATH=${CMAKE_BINARY_DIR}/RecCalorimeter/${GAUDI_GENCONF_DIR}:${CMAKE_BINARY_DIR}/RecFCCeeCalorimeter/${GAUDI_GENCONF_DIR}:${CMAKE_BINARY_DIR}/RecFCChhCalorimeter/${GAUDI_GENCONF_DIR}:$ENV{PYTHONPATH}:$<TARGET_FILE_DIR:k4FWCore::k4FWCore>/../python:$ENV{PYTHONPATH})

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               PATH=$<TARGET_FILE_DIR:k4FWCore::k4FWCore>/../bin:$ENV{PATH}
  )

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               K4RECCALORIMETER=${CMAKE_CURRENT_LIST_DIR}/
  )
endfunction()


add_subdirectory(EDM4hep)
add_subdirectory(RecCalorimeter)
add_subdirectory(RecFCChhCalorimeter)
add_subdirectory(RecFCCeeCalorimeter)





install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${PROJECT_NAME}Targets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/"
  )

gaudi_install(CMAKE cmake/${PROJECT_NAME}Config.cmake)

