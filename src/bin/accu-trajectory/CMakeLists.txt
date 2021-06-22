# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_TRAJECTORY_SRC
        ABFTrajectory.cpp
        ABFTrajOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-trajectory ${ABF_TRAJECTORY_SRC})
ADD_DEPENDENCIES(abf-trajectory cpmf_shared)

TARGET_LINK_LIBRARIES(abf-trajectory  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-trajectory
        RUNTIME DESTINATION
            bin
        )

