CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

project(k4RecCalorimeter)

include(GNUInstallDirs)
include(CTest)
include(cmake/Key4hepConfig.cmake)

# Set default for install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/InstallArea/ CACHE PATH
    "Install path prefix, prepended onto install directories." FORCE )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})


#---------------------------------------------------------------
find_package(ROOT COMPONENTS RIO Tree)
find_package(Gaudi)
find_package(k4FWCore) # implicit: Gaudi
find_package(EDM4HEP) # implicit: Podio
find_package(DD4hep)
find_package(k4geo)
find_package(FastJet)
find_package(sipm)
# New versions of ONNRuntime package provide onnxruntime-Config.cmake
# and use the name onnxruntime
find_package(onnxruntime)
if (NOT onnxruntime_FOUND)
  message(STATUS "Could not find onnxruntime (> 1.17.1). Looking for an older version")
  find_package(ONNXRuntime)
endif()

if(onnxruntime_FOUND OR ONNXRuntime_FOUND)
else()
  message(FATAL_ERROR "Failed to locate ONNXRuntime!")
endif()

#---------------------------------------------------------------

set(GAUDI_GENCONF_DIR "genConfDir")

function(set_test_env _testname)
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$<TARGET_FILE_DIR:k4RecCalorimeterPlugins>:${CMAKE_BINARY_DIR}:$<TARGET_FILE_DIR:k4RecFCCeeCalorimeterPlugins>:${CMAKE_BINARY_DIR}:$<TARGET_FILE_DIR:k4RecFCChhCalorimeterPlugins>:$<TARGET_FILE_DIR:ROOT::Core>:$<TARGET_FILE_DIR:k4FWCore::k4FWCore>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:podio::podio>:$ENV{LD_LIBRARY_PATH})

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               PYTHONPATH=${CMAKE_BINARY_DIR}/RecCalorimeter/${GAUDI_GENCONF_DIR}:${CMAKE_BINARY_DIR}/RecFCCeeCalorimeter/${GAUDI_GENCONF_DIR}:${CMAKE_BINARY_DIR}/RecFCChhCalorimeter/${GAUDI_GENCONF_DIR}:$ENV{PYTHONPATH}:$<TARGET_FILE_DIR:k4FWCore::k4FWCore>/../python:$ENV{PYTHONPATH})

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               PATH=$<TARGET_FILE_DIR:k4FWCore::k4FWCore>/../bin:$ENV{PATH}
  )

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
               K4RECCALORIMETER=${CMAKE_CURRENT_LIST_DIR}/
  )
endfunction()


add_subdirectory(RecCaloCommon)
add_subdirectory(RecCalorimeter)
add_subdirectory(RecFCChhCalorimeter)
add_subdirectory(RecFCCeeCalorimeter)


install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${PROJECT_NAME}Targets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/"
  )

gaudi_install(CMAKE cmake/${PROJECT_NAME}Config.cmake)
