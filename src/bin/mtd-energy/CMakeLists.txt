# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(MTD_ENERGY_SRC
        MTDEnergy.cpp
        MTDEneOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(mtd-energy ${MTD_ENERGY_SRC})

TARGET_LINK_LIBRARIES(mtd-energy ${CPMF_LIBRARIES})

INSTALL(TARGETS
            mtd-energy
        RUNTIME DESTINATION
            bin
        )

