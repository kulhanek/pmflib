# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABP_TRAJECTORY_SRC
        ABPTrajectory.cpp
        ABPTrajOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abp-trajectory ${ABP_TRAJECTORY_SRC})
ADD_DEPENDENCIES(abp-trajectory cpmf_shared)

TARGET_LINK_LIBRARIES(abp-trajectory  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abp-trajectory
        RUNTIME DESTINATION
            bin
        )

