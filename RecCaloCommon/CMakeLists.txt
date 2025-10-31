################################################################################
# Package: RecCaloCommon
################################################################################

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
