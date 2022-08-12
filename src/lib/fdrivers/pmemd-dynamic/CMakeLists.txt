# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(PMEMD_DRV_SRC_D01
        pmf_pmemd_dat_d01.F90
        pmf_pmemd_control_d01.f90
        pmf_pmemd_d01.F90
        )

SET(PMEMD_DRV_CLIENT_SRC_D01
        PMFLibPMEMD_d01.F90
        )


# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(pmemd_drv_static_d01 STATIC ${PMEMD_DRV_SRC_D01})
    SET_TARGET_PROPERTIES(pmemd_drv_static_d01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_d01
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(pmemd_drv_shared_d01 SHARED ${PMEMD_DRV_SRC_D01})
    ADD_DEPENDENCIES(pmemd_drv_shared_d01 fpmf_shared)
    TARGET_LINK_LIBRARIES(pmemd_drv_shared_d01 ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(pmemd_drv_shared_d01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_d01
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                pmemd_drv_shared_d01
            DESTINATION
                lib/drivers
            )

    # client test
    ADD_LIBRARY(pmemd_drv_client_shared_d01 SHARED ${PMEMD_DRV_CLIENT_SRC_D01})
    TARGET_LINK_LIBRARIES(pmemd_drv_client_shared_d01 ${FPMF_LIBRARIES})
    TARGET_COMPILE_DEFINITIONS(pmemd_drv_client_shared_d01 PUBLIC PMFLIB)
    SET_TARGET_PROPERTIES(pmemd_drv_client_shared_d01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_client_d01
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
ENDIF(LIBS_SHARED)

