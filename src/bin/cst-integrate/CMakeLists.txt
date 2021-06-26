# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(CST_INTEGRATE_SRC
        CSTEnergyIntegrate.cpp
        CSTEnergyIntOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _cst-integrate)
ELSE()
    SET(TARGET_NAME cst-integrate)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${CST_INTEGRATE_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

