# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(CV_MATH_TEST_SRC
        cv_math_test.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(cv-math-test ${CV_MATH_TEST_SRC})
ADD_DEPENDENCIES(cv-math-test fpmf_shared)

TARGET_LINK_LIBRARIES(cv-math-test
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
            cv-math-test
        RUNTIME DESTINATION
            bin
        )
