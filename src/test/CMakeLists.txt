# ==============================================================================
# PMFLib CMake File
# ==============================================================================

IF(PMFLIB_BUILD_UTILS)
    ADD_SUBDIRECTORY(test-cv-math)
    ADD_SUBDIRECTORY(test-gpr-kernel)
    ADD_SUBDIRECTORY(test-cv-splines)
ENDIF(PMFLIB_BUILD_UTILS)
