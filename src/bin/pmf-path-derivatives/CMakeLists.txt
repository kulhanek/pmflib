# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMF_PR_SRC
        pmf_path_derivatives.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-path-derivatives ${PMF_PR_SRC})

TARGET_LINK_LIBRARIES(pmf-path-derivatives
        fpmf
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            pmf-path-derivatives
        DESTINATION
            bin
        )
