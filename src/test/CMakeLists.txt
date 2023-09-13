# ==============================================================================
# PMFLib CMake File
# ==============================================================================

IF(PMFLIB_BUILD_UTILS)
    ADD_SUBDIRECTORY(test-cv-math)
    ADD_SUBDIRECTORY(test-gpr-kernel)
ENDIF(PMFLIB_BUILD_UTILS)
