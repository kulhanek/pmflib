# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(BM_INTEGRATE_SRC
        BMIntegrate.cpp
        BMIntOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-integrate ${BM_INTEGRATE_SRC})
ADD_DEPENDENCIES(pmf-integrate cpmf_shared)

TARGET_LINK_LIBRARIES(pmf-integrate  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            pmf-integrate
        RUNTIME DESTINATION
            bin
        )

