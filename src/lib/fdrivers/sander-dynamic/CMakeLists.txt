# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(SANDER_DRV_SRC_D02
        pmf_sander_dat_d02.F90
        pmf_sander_control_d02.f90
        pmf_sander_d02.F90
        )

SET(SANDER_DRV_CLIENT_SRC_D02
        PMFLibSANDER_d02.F90
        )


# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(sander_drv_static_d02 STATIC ${SANDER_DRV_SRC_D02})
    SET_TARGET_PROPERTIES(sander_drv_static_d02 PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander_d02
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(sander_drv_shared_d02 SHARED ${SANDER_DRV_SRC_D02})
    ADD_DEPENDENCIES(sander_drv_shared_d02 fpmf_shared)
    TARGET_LINK_LIBRARIES(sander_drv_shared_d02 ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(sander_drv_shared_d02 PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander_d02
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                sander_drv_shared_d02
            DESTINATION
                lib/drivers
            )

    # client test
    ADD_LIBRARY(sander_drv_client_shared_d02 SHARED ${SANDER_DRV_CLIENT_SRC_D02})
    ADD_DEPENDENCIES(sander_drv_client_shared_d02 fpmf_shared)
    TARGET_LINK_LIBRARIES(sander_drv_client_shared_d02 ${FPMF_LIBRARIES})
    TARGET_COMPILE_DEFINITIONS(sander_drv_client_shared_d02 PUBLIC PMFLIB)
    SET_TARGET_PROPERTIES(sander_drv_client_shared_d02 PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander_client_d02
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
ENDIF(LIBS_SHARED)

