# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(XDYNBP_DRV_OBJ
        pmf_xdynbp.f90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(xdynbp_drv_static STATIC ${XDYNBP_DRV_OBJ})
    SET_TARGET_PROPERTIES(xdynbp_drv_static PROPERTIES
                            OUTPUT_NAME fpmfdrv_xdynbp
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(xdynbp_drv_shared SHARED ${XDYNBP_DRV_OBJ})
    ADD_DEPENDENCIES(xdynbp_drv_shared fpmf_shared)
    TARGET_LINK_LIBRARIES(xdynbp_drv_shared ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(xdynbp_drv_shared PROPERTIES
                            OUTPUT_NAME fpmfdrv_xdynbp
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    TARGET_LINK_LIBRARIES(xdynbp_drv_shared
            fpmf
            )
    INSTALL(TARGETS
                xdynbp_drv_shared
            DESTINATION
                lib/drivers
            )
ENDIF(LIBS_SHARED)
