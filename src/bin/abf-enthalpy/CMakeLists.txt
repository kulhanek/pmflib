# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_ENTHALPY_SRC
        ABFEnthalpy.cpp
        ABFEnthalpyOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-enthalpy ${ABF_ENTHALPY_SRC})
ADD_DEPENDENCIES(abf-enthalpy cpmf_shared)

TARGET_LINK_LIBRARIES(abf-enthalpy  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-enthalpy
        RUNTIME DESTINATION
            bin
        )

