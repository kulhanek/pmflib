# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(AUS_INTEGRATE_SRC
        AUSEnergyIntegrate.cpp
        AUSEnergyIntOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _aus-integrate)
ELSE()
    SET(TARGET_NAME aus-integrate)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${AUS_INTEGRATE_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

