# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(CVSPLINES_TEST_SRC
        test-cv-splines.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(test-cv-splines ${CVSPLINES_TEST_SRC})
ADD_DEPENDENCIES(test-cv-splines cpmf_shared)

TARGET_LINK_LIBRARIES(test-cv-splines
        cpmf
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${CPMF_LIBRARIES}
        )

INSTALL(TARGETS
            test-cv-splines
        RUNTIME DESTINATION
            bin/tests
        )
