# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(TEST_COORDS_SRC
        test_coords.f90
        test_coords_dat.f90
        test_coords_control.f90
        test_coords_utils.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(test-cvs ${TEST_COORDS_SRC})
ADD_DEPENDENCIES(test-cvs fpmf_shared)

TARGET_LINK_LIBRARIES(test-cvs
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
            test-cvs
        RUNTIME DESTINATION
            bin/tests
        )
