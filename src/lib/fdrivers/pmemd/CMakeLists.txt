# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(PMEMD_DRV_SRC_V01
        pmf_pmemd_dat_v01.F90
        pmf_pmemd_control_v01.f90
        pmf_pmemd_v01.F90
        PMFLibPMEMD_v01.F90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(pmemd_drv_static_v01 STATIC ${PMEMD_DRV_SRC_V01})
    SET_TARGET_PROPERTIES(pmemd_drv_static_v01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_v01
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(pmemd_drv_shared_v01 SHARED ${PMEMD_DRV_SRC_V01})
    ADD_DEPENDENCIES(pmemd_drv_shared_v01 fpmf_shared)
    TARGET_LINK_LIBRARIES(pmemd_drv_shared_v01 ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(pmemd_drv_shared_v01 PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd_v01
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                pmemd_drv_shared_v01
            DESTINATION
                lib/drivers
            )
ENDIF(LIBS_SHARED)

