# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_OPTGPRHPRMS_SRC
        ABFEnergyOptGPRHyprms.cpp
        ABFEnergyOptGPRHyprmsOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _abf-energy-optgprhyprms)
ELSE()
    SET(TARGET_NAME abf-energy-optgprhyprms)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${ABF_OPTGPRHPRMS_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${LBFGS_LIB_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

