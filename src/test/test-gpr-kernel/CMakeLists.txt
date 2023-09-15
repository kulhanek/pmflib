# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(KERNEL_TEST_SRC
        gpr-test-kernel.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(test-gpr-kernel ${KERNEL_TEST_SRC})
ADD_DEPENDENCIES(test-gpr-kernel fpmf_shared)

TARGET_LINK_LIBRARIES(test-gpr-kernel
        cpmf
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
            test-gpr-kernel
        RUNTIME DESTINATION
            bin
        )
