# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_MASK_SRC
        ABFMask.cpp
        ABFMaskOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-mask ${ABF_MASK_SRC})
ADD_DEPENDENCIES(abf-mask cpmf_shared)

TARGET_LINK_LIBRARIES(abf-mask ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-mask
        RUNTIME DESTINATION
            bin
        )

