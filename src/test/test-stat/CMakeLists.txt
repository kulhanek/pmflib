# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(TEST_STAT_SRC
        test_stat.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(test-stat ${TEST_STAT_SRC})
ADD_DEPENDENCIES(test-stat fpmf_shared)

TARGET_LINK_LIBRARIES(test-stat
        fpmf
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${RANLUX_LIB_NAME}
        ${CPMF_LIBRARIES}
        )

INSTALL(TARGETS
            test-stat
        RUNTIME DESTINATION
            bin/tests
        )
