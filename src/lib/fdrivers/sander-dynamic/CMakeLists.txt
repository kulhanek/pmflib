# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(SANDER_DRV_SRC_D01
        pmf_sander_dat_d01.F90
        pmf_sander_control_d01.f90
        pmf_sander_d01.F90
        )

SET(SANDER_DRV_CLIENT_SRC_D01
        PMFLibSANDER_d01.F90
        )


# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(sander_drv_static_d01 STATIC ${SANDER_DRV_SRC_D01})
    SET_TARGET_PROPERTIES(sander_drv_static_d01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander_d01
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(sander_drv_shared_d01 SHARED ${SANDER_DRV_SRC_D01})
    ADD_DEPENDENCIES(sander_drv_shared_d01 fpmf_shared)
    TARGET_LINK_LIBRARIES(sander_drv_shared_d01 ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(sander_drv_shared_d01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander_d01
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                sander_drv_shared_d01
            DESTINATION
                lib/drivers
            )

    # client test
    ADD_LIBRARY(sander_drv_client_shared_d01 SHARED ${SANDER_DRV_CLIENT_SRC_D01})
    ADD_DEPENDENCIES(sander_drv_client_shared_d01 fpmf_shared)
    TARGET_LINK_LIBRARIES(sander_drv_client_shared_d01 ${FPMF_LIBRARIES})
    TARGET_COMPILE_DEFINITIONS(sander_drv_client_shared_d01 PUBLIC PMFLIB)
    SET_TARGET_PROPERTIES(sander_drv_client_shared_d01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander_client_d01
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
ENDIF(LIBS_SHARED)

