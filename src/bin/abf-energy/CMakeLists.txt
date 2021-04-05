# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_INTEGRATE_SRC
        ABFIntegrate.cpp
        ABFIntOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _abf-integrate)
ELSE()
    SET(TARGET_NAME abf-integrate)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${ABF_INTEGRATE_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

