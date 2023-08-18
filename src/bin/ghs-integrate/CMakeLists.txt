# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(GHS_INTEGRATE_SRC
        GHSEnergyIntegrate.cpp
        GHSEnergyIntOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _ghs-integrate)
ELSE()
    SET(TARGET_NAME ghs-integrate)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${GHS_INTEGRATE_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

