# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ACCU_COMBINE_SRC
        ACCUCombine.cpp
        ACCUCombOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(accu-combine ${ACCU_COMBINE_SRC})
ADD_DEPENDENCIES(accu-combine cpmf_shared)

TARGET_LINK_LIBRARIES(accu-combine  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            accu-combine
        RUNTIME DESTINATION
            bin
        )

