################################################################################
# Package: RecCaloCommon
################################################################################

file(GLOB _sources src/*.cpp )

# include_directories(${FASTJET_INCLUDE_DIRS})
# link_directories(${FASTJET_LIBRARIES})

gaudi_add_library(RecCaloCommon
                  SOURCES ${_sources}
                  LINK Gaudi::GaudiKernel
                       EDM4HEP::edm4hep
                       ${FASTJET_LIBRARIES}
)

message(STATUS "FastJet include directories:\n${FASTJET_INCLUDE_DIRS}")
message(STATUS "FastJet libraries:\n${FASTJET_LIBRARIES}")

target_include_directories(RecCaloCommon PUBLIC ${FASTJET_INCLUDE_DIRS})
target_link_directories(RecCaloCommon PUBLIC ${FASTJET_LIBRARIES})

install(TARGETS RecCaloCommon
        EXPORT k4RecCalorimeterTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
                                                      COMPONENT dev
)
