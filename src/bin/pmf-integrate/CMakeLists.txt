# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMF_INTEGRATE_SRC
        PMFEnergyIntegrate.cpp
        PMFEnergyIntOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _pmf-integrate)
ELSE()
    SET(TARGET_NAME pmf-integrate)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${PMF_INTEGRATE_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

