# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMFENERGY_SRC
        PMFEnergy.cpp
        PMFEneOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _pmf-energy)
ELSE()
    SET(TARGET_NAME pmf-energy)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${PMFENERGY_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

