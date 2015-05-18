# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMF_GP_SRC
        pmf_gprocess.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-test-gprocess ${PMF_GP_SRC})

TARGET_LINK_LIBRARIES(pmf-test-gprocess
        fpmf
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            pmf-test-gprocess
        DESTINATION
            bin
        )
