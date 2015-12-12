# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_INTEGRATE_SRC
        ABFIntegrate.cpp
        ABFIntOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-integrate ${ABF_INTEGRATE_SRC})
ADD_DEPENDENCIES(pmf-integrate cpmf_shared)

TARGET_LINK_LIBRARIES(abf-integrate ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-integrate
        RUNTIME DESTINATION
            bin
        )

