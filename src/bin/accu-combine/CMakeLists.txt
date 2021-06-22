# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_COMBINE_SRC
        ABFCombine.cpp
        ABFCombOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-combine ${ABF_COMBINE_SRC})
ADD_DEPENDENCIES(abf-combine cpmf_shared)

TARGET_LINK_LIBRARIES(abf-combine  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-combine
        RUNTIME DESTINATION
            bin
        )

