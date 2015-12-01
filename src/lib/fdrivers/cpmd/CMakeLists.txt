# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(CPMD_DRV_SRC
        pmf_cpmd_common.f90
        pmf_cpmd_dat.f90
        pmf_cpmd_control.f90
        pmf_cpmd.F90
        pmf_cpmd_qmmm.f90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(cpmd_drv_static STATIC ${CPMD_DRV_SRC})
    SET_TARGET_PROPERTIES(cpmd_drv_static PROPERTIES
                            OUTPUT_NAME fpmfdrv_cpmd
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(cpmd_drv_shared SHARED ${CPMD_DRV_SRC})
    TARGET_LINK_LIBRARIES(cpmd_drv_shared ${FPMF_LIBRARIES})
    SET_TARGET_PROPERTIES(cpmd_drv_shared PROPERTIES
                            OUTPUT_NAME fpmfdrv_cpmd
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                cpmd_drv_shared
            DESTINATION
                lib/drivers
            )
ENDIF(LIBS_SHARED)
