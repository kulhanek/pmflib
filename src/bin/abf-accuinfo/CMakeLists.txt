# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_ACCUINFO_SRC
        ABFAccuInfo.cpp
        ABFAccuInfoOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-accuinfo ${ABF_ACCUINFO_SRC})
ADD_DEPENDENCIES(abf-accuinfo cpmf_shared)

TARGET_LINK_LIBRARIES(abf-accuinfo ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-accuinfo
        RUNTIME DESTINATION
            bin
        )

