# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(PMEMD_DRV_SRC_D02
        pmf_pmemd_dat_d02.F90
        pmf_pmemd_control_d02.f90
        pmf_pmemd_d02.F90
        )

SET(PMEMD_DRV_CLIENT_SRC_D02
        PMFLibPMEMD_d02.F90
        )


# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(pmemd_drv_static_d02 STATIC ${PMEMD_DRV_SRC_D02})
    SET_TARGET_PROPERTIES(pmemd_drv_static_d02 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_d02
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(pmemd_drv_shared_d02 SHARED ${PMEMD_DRV_SRC_D02})
    ADD_DEPENDENCIES(pmemd_drv_shared_d02 fpmf_shared)
    TARGET_LINK_LIBRARIES(pmemd_drv_shared_d02 ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(pmemd_drv_shared_d02 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_d02
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                pmemd_drv_shared_d02
            DESTINATION
                lib/drivers
            )

    # client test
    ADD_LIBRARY(pmemd_drv_client_shared_d02 SHARED ${PMEMD_DRV_CLIENT_SRC_D02})
    ADD_DEPENDENCIES(pmemd_drv_client_shared_d02 fpmf_shared)
    TARGET_LINK_LIBRARIES(pmemd_drv_client_shared_d02 ${FPMF_LIBRARIES})
    TARGET_COMPILE_DEFINITIONS(pmemd_drv_client_shared_d02 PUBLIC PMFLIB)
    SET_TARGET_PROPERTIES(pmemd_drv_client_shared_d02 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_client_d02
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
ENDIF(LIBS_SHARED)

