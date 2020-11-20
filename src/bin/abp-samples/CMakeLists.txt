# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABP_SAMPLES_SRC
        ABPSamples.cpp
        ABPSampOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abp-samples ${ABP_SAMPLES_SRC})
ADD_DEPENDENCIES(abp-samples cpmf_shared)

TARGET_LINK_LIBRARIES(abp-samples  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abp-samples
        RUNTIME DESTINATION
            bin
        )

