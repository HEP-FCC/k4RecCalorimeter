################################################################################
# Package: RecCaloCommon
################################################################################

find_package(boost_timer REQUIRED)
file(GLOB _sources src/*.cpp )

gaudi_add_library(RecCaloCommon
                  SOURCES ${_sources}
                  LINK Gaudi::GaudiKernel
                       EDM4HEP::edm4hep
                       k4FWCore::k4FWCore
                       k4FWCore::k4Interface
                       DD4hep::DDCore
                       ${FASTJET_LIBRARIES}
)

target_include_directories(RecCaloCommon PUBLIC ${FASTJET_INCLUDE_DIRS})
target_link_directories(RecCaloCommon PUBLIC ${FASTJET_LIBRARY_DIRS})

install(TARGETS RecCaloCommon
        EXPORT k4RecCalorimeterTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
                                                      COMPONENT dev
)


if(BUILD_TESTING)
  gaudi_add_executable(IDMap_test.exe
    SOURCES tests/IDMap_test.cpp
    LINK DD4hep::DDCore Boost::timer
    TEST)
  target_include_directories(IDMap_test.exe PUBLIC include)


  gaudi_add_executable(IDMapIndexer_test.exe
    SOURCES tests/IDMapIndexer_test.cpp
    LINK DD4hep::DDCore k4FWCore::k4Interface k4FWCore::k4FWCore
    TEST)
  target_include_directories(IDMapIndexer_test.exe AFTER PUBLIC include)


  gaudi_add_executable(MultiIndexer_test.exe
    SOURCES tests/MultiIndexer_test.cpp src/MultiIndexer.cpp
    LINK DD4hep::DDCore k4FWCore::k4Interface k4FWCore::k4FWCore
    TEST)
  target_include_directories(MultiIndexer_test.exe AFTER PUBLIC include)
endif()
