# ==============================================================================
# PRMFile CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(SANDER_DRV_SRC
        pmf_sander_dat.F90
        pmf_sander_control.f90
        pmf_sander.F90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(sander_drv_static STATIC ${SANDER_DRV_SRC})
    SET_TARGET_PROPERTIES(sander_drv_static PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(sander_drv_shared SHARED ${SANDER_DRV_SRC})
    SET_TARGET_PROPERTIES(sander_drv_shared PROPERTIES
                            OUTPUT_NAME fpmfdrv_sander
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                sander_drv_shared
            LIBRARY DESTINATION
                lib/drivers
            )
ENDIF(LIBS_SHARED)
