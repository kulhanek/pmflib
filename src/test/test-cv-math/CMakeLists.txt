# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(CV_MATH_TEST_SRC
        cv_math_test.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(test-cv-math ${CV_MATH_TEST_SRC})
ADD_DEPENDENCIES(test-cv-math fpmf_shared)

TARGET_LINK_LIBRARIES(test-cv-math
        fpmf
        ${XBPLIB_XCORE_LIB_NAME}
        ${XBPLIB_XFILES_LIB_NAME}
        ${XBPLIB_XENE_LIB_NAME}
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${CPMF_LIBRARIES}
        )

INSTALL(TARGETS
            test-cv-math
        RUNTIME DESTINATION
            bin
        )
