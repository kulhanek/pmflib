# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMF_EXPLORE_SRC
        system_dat.f90
        system_control.f90
        system.f90
        explore.f90
        pmf_path_explore.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-path-explore ${PMF_EXPLORE_SRC})

TARGET_LINK_LIBRARIES(pmf-path-explore
        fpmf
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${RANLUX_LIB_NAME}
        ${RANLUX_LIB_NAME}
        )

INSTALL(TARGETS
            pmf-path-explore
        DESTINATION
            bin
        )
