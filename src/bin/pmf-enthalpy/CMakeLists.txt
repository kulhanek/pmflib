# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ENTHALPY_SRC
        Enthalpy.cpp
        EnthalpyOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-enthalpy ${ENTHALPY_SRC})
ADD_DEPENDENCIES(pmf-enthalpy cpmf_shared)

TARGET_LINK_LIBRARIES(pmf-enthalpy  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            pmf-enthalpy
        RUNTIME DESTINATION
            bin
        )

