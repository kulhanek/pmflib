# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(TEST_STAT_SRC
        test_stat.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-test-stat ${TEST_STAT_SRC})
ADD_DEPENDENCIES(pmf-test-stat fpmf_shared)

TARGET_LINK_LIBRARIES(pmf-test-stat
        fpmf
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${RANLUX_LIB_NAME}
        ${CPMF_LIBRARIES}
        )

INSTALL(TARGETS
            pmf-test-stat
        RUNTIME DESTINATION
            bin
        )
