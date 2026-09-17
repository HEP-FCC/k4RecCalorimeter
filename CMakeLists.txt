CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

project(k4RecCalorimeter LANGUAGES CXX)

cmake_policy(SET CMP0144 NEW)
cmake_policy(SET CMP0167 NEW)

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
find_package(ROOT COMPONENTS RIO Tree REQUIRED)
find_package(Gaudi REQUIRED)
find_package(k4FWCore 1.3.0 REQUIRED)
find_package(EDM4HEP REQUIRED) # implicit: Podio
find_package(DD4hep REQUIRED)
find_package(k4geo 0.23 REQUIRED)
find_package(FastJet REQUIRED)
find_package(sipm REQUIRED)
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



add_subdirectory(RecCaloCommon)
add_subdirectory(RecCalorimeter)
add_subdirectory(RecFCChhCalorimeter)
add_subdirectory(RecFCCeeCalorimeter)
add_subdirectory(RecCaloDigi)


install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${PROJECT_NAME}Targets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/"
  )

gaudi_install(CMAKE cmake/${PROJECT_NAME}Config.cmake)
