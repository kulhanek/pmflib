# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ENTHALPY_SRC
        Enthalpy.cpp
        EnthalpyOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _pmf-enthalpy)
ELSE()
    SET(TARGET_NAME pmf-enthalpy)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${ENTHALPY_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

