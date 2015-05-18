# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_SAMPLES_SRC
        ABFSamples.cpp
        ABFSampOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-samples ${ABF_SAMPLES_SRC})

TARGET_LINK_LIBRARIES(abf-samples  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-samples
        RUNTIME DESTINATION
            bin
        )

