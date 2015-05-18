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
ADD_EXECUTABLE(pmf-test-cvs ${TEST_COORDS_SRC})

TARGET_LINK_LIBRARIES(pmf-test-cvs
        fpmf
        ${XBPLIB_XCORE_LIB_NAME}
        ${XBPLIB_XFILES_LIB_NAME}
        ${XBPLIB_XENE_LIB_NAME}
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${RANLUX_LIB_NAME}
        ${CPMF_LIBRARIES}
        )

INSTALL(TARGETS
            pmf-test-cvs
        RUNTIME DESTINATION
            bin
        )
