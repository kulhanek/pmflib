# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(RST_SAMPLES_SRC
        RSTSamples.cpp
        RSTSamplesOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(rst-samples ${RST_SAMPLES_SRC})

TARGET_LINK_LIBRARIES(rst-samples ${CPMF_LIBRARIES})

INSTALL(TARGETS
            rst-samples
        RUNTIME DESTINATION
            bin
        )

