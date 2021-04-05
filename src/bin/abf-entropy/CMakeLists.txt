# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_ENTROPY_INTEGRATE_SRC
        ABFEntropyIntegrate.cpp
        ABFEntropyIntOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _abf-entropy)
ELSE()
    SET(TARGET_NAME abf-entropy)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${ABF_ENTROPY_INTEGRATE_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

