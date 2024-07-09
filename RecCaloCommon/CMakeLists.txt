################################################################################
# Package: RecCaloCommon
################################################################################

file(GLOB _sources src/*.cpp )

include_directories(${FASTJET_INCLUDE_DIRS})
link_directories(${FASTJET_LIBRARIES})

gaudi_add_library(RecCaloCommon
                  SOURCES ${_sources}
                  LINK Gaudi::GaudiKernel
                       EDM4HEP::edm4hep
                       ${FASTJET_LIBRARIES}
)

install(TARGETS RecCaloCommon
        EXPORT k4RecCalorimeterTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
                                                      COMPONENT dev
)
