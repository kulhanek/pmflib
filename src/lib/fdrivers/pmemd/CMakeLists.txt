# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(PMEMD_DRV_SRC
        pmf_pmemd_dat.F90
        pmf_pmemd_control.f90
        pmf_pmemd.F90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(pmemd_drv_static STATIC ${PMEMD_DRV_SRC})
    SET_TARGET_PROPERTIES(pmemd_drv_static PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(pmemd_drv_shared SHARED ${PMEMD_DRV_SRC})
    ADD_DEPENDENCIES(pmemd_drv_shared fpmf_shared)
    TARGET_LINK_LIBRARIES(pmemd_drv_shared ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(pmemd_drv_shared PROPERTIES
                            OUTPUT_NAME fpmfdrv_pmemd
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                pmemd_drv_shared
            DESTINATION
                lib/drivers
            )
ENDIF(LIBS_SHARED)

