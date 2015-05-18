# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(STM_TRAJECTORY_SRC
        STMTrajectory.cpp
        STMTrajOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(stm-trajectory ${STM_TRAJECTORY_SRC})

TARGET_LINK_LIBRARIES(stm-trajectory  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            stm-trajectory
        RUNTIME DESTINATION
            bin
        )

