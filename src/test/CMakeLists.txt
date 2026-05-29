# ==============================================================================
# PMFLib CMake File
# ==============================================================================

IF(PMFLIB_BUILD_UTILS)
    SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/tests)
    ADD_SUBDIRECTORY(test-cv-math)
    ADD_SUBDIRECTORY(test-gpr-kernel)
    ADD_SUBDIRECTORY(test-cv-splines)
    ADD_SUBDIRECTORY(test-coords)
    ADD_SUBDIRECTORY(test-stat)
ENDIF(PMFLIB_BUILD_UTILS)
