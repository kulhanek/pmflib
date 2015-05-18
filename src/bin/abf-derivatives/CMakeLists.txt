# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_DERIVATIVES_SRC
        ABFDerivatives.cpp
        ABFDeriOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-derivatives ${ABF_DERIVATIVES_SRC})

TARGET_LINK_LIBRARIES(abf-derivatives  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-derivatives
        RUNTIME DESTINATION
            bin
        )

