# ==============================================================================
# PMFLib CMake File
# ==============================================================================

SET(FPMF_DYN_SRC
        pmfdyn_dynamics.f90
        pmfdyn_system_control.f90
        pmfdyn_system_dat.f90
        pmfdyn_system.f90
        pmfdyn_restraints_types.f90
        pmfdyn_restraints_control.f90
        pmfdyn_restraints_dat.f90
        pmfdyn_restraints.f90
        pmfdyn_thermostat_control.f90
        pmfdyn_thermostat_dat.f90
        pmfdyn_thermostat.f90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(fpmfdyn_static STATIC ${FPMF_DYN_SRC})
    SET_TARGET_PROPERTIES(fpmfdyn_static PROPERTIES
                            OUTPUT_NAME fpmfdyn
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(fpmfdyn_shared SHARED ${FPMF_DYN_SRC})
    SET_TARGET_PROPERTIES(fpmfdyn_shared PROPERTIES
                            OUTPUT_NAME fpmfdyn
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION}
                            )
    TARGET_LINK_LIBRARIES(fpmfdyn_shared
            fpmf
            ${CPMF_LIBRARIES}
            ${XBPLIB_XCORE_LIB_NAME}
            ${XBPLIB_XFILES_LIB_NAME}
            ${XBPLIB_XENE_LIB_NAME}
            ${PRMFILE_FLIB_NAME}
            ${LAPACK_LIB_NAME}
            ${BLAS_LIB_NAME}
            ${SCIMAFIC_FLIB_NAME}
            ${RANLUX_LIB_NAME}
            ${CPMF_LIBRARIES}
            )
    INSTALL(TARGETS
                fpmfdyn_shared
            DESTINATION
                lib)
ENDIF(LIBS_SHARED)



