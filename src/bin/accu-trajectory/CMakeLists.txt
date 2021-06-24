# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ACCU_TRAJECTORY_SRC
        ACCUTrajectory.cpp
        ACCUTrajOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(accu-trajectory ${ACCU_TRAJECTORY_SRC})
ADD_DEPENDENCIES(accu-trajectory cpmf_shared)

TARGET_LINK_LIBRARIES(accu-trajectory  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            accu-trajectory
        RUNTIME DESTINATION
            bin
        )

