# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_OPTENTGPRHPRMS_SRC
        ABFEnthalpyOptGPRHyprms.cpp
        ABFEnthalpyOptGPRHyprmsOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _abf-enthalpy-optgprhyprms)
ELSE()
    SET(TARGET_NAME abf-enthalpy-optgprhyprms)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${ABF_OPTENTGPRHPRMS_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${LBFGS_LIB_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

