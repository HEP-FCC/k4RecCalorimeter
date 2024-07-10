################################################################################
# Package: RecCaloCommon
################################################################################

file(GLOB _sources src/*.cpp )

gaudi_add_library(RecCaloCommon
                  SOURCES ${_sources}
                  LINK Gaudi::GaudiKernel
                       EDM4HEP::edm4hep
                       ${FASTJET_LIBRARIES}
)

message(STATUS "FastJet include directories:\n${FASTJET_INCLUDE_DIRS}")
message(STATUS "FastJet libraries:\n${FASTJET_LIBRARIES}")
message(STATUS "FastJet library:\n${FASTJET_LIBRARY}")
message(STATUS "FastJet library directories:\n${FASTJET_LIBRARY_DIRS}")

target_include_directories(RecCaloCommon PUBLIC ${FASTJET_INCLUDE_DIRS})
target_link_directories(RecCaloCommon PUBLIC ${FASTJET_LIBRARY_DIRS})

install(TARGETS RecCaloCommon
        EXPORT k4RecCalorimeterTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
                                                      COMPONENT dev
)
