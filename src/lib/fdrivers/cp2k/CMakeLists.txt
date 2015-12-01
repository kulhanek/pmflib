# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(CP2K_DRV_SRC
        pmf_cp2k_dat.f90
        pmf_cp2k_control.f90
        pmf_cp2k.f90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(cp2k_drv_static STATIC ${CP2K_DRV_SRC})
    SET_TARGET_PROPERTIES(cp2k_drv_static PROPERTIES
                            OUTPUT_NAME fpmfdrv_cp2k
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(cp2k_drv_shared SHARED ${CP2K_DRV_SRC})
    TARGET_LINK_LIBRARIES(cp2k_drv_shared ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(cp2k_drv_shared PROPERTIES
                            OUTPUT_NAME fpmfdrv_cp2k
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                cp2k_drv_shared
            DESTINATION
                lib/drivers
            )
ENDIF(LIBS_SHARED)
