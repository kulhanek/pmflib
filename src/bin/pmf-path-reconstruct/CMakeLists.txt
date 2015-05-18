# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMF_PR_SRC
        pmf_path_reconstruct.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-path-reconstruct ${PMF_PR_SRC})

TARGET_LINK_LIBRARIES(pmf-path-reconstruct
        fpmf
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            pmf-path-reconstruct
        DESTINATION
            bin
        )
