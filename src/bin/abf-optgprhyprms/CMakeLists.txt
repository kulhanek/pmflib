# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_OPTGPRHPRMS_SRC
        ABFOptGPRHyprms.cpp
        ABFOptGPRHyprmsOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-optgprhyprms ${ABF_OPTGPRHPRMS_SRC})
ADD_DEPENDENCIES(abf-optgprhyprms cpmf_shared)

TARGET_LINK_LIBRARIES(abf-optgprhyprms ${LBFGS_LIB_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-optgprhyprms
        RUNTIME DESTINATION
            bin
        )

