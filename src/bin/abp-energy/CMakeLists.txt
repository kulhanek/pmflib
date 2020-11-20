# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABP_ENERGY_SRC
        ABPEnergy.cpp
        ABPEneOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abp-energy ${ABP_ENERGY_SRC})
ADD_DEPENDENCIES(abp-energy cpmf_shared)

TARGET_LINK_LIBRARIES(abp-energy ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abp-energy
        RUNTIME DESTINATION
            bin
        )

