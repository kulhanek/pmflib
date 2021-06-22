# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMF_ACCUINFO_SRC
        PMFAccuInfo.cpp
        PMFAccuInfoOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(accu-info ${PMF_ACCUINFO_SRC})
ADD_DEPENDENCIES(accu-info cpmf_shared)

TARGET_LINK_LIBRARIES(accu-info  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            accu-info
        RUNTIME DESTINATION
            bin
        )

