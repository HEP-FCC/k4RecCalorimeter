################################################################################
# Package: RecCaloInterface
################################################################################

gaudi_add_header_only_library(RecCaloInterface)

install(TARGETS RecCaloInterface
        EXPORT ${PROJECT_NAME}Targets
        COMPONENT dev)
