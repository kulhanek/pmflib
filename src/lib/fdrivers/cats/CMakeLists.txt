# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(CATS_DRV_SRC
        PMFCATsDriver.cpp
        pmf_cats.f90
        pmf_cats_dat.f90
        pmf_cats_control.f90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(cats_drv_static STATIC ${CATS_DRV_SRC})
    SET_TARGET_PROPERTIES(cats_drv_static PROPERTIES
                            OUTPUT_NAME fpmfdrv_cats
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(cats_drv_shared SHARED ${CATS_DRV_SRC})
    TARGET_LINK_LIBRARIES(cats_drv_shared ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(cats_drv_shared PROPERTIES
                            OUTPUT_NAME fpmfdrv_cats
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                cats_drv_shared
            LIBRARY DESTINATION
                lib/drivers
            )
ENDIF(LIBS_SHARED)
