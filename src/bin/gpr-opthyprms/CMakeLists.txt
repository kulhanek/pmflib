# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(OPTGPRHPRMS_SRC
        OptGPRHyprms.cpp
        OptGPRHyprmsOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _gpr-opthyprms)
ELSE()
    SET(TARGET_NAME gpr-opthyprms)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${OPTGPRHPRMS_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${LBFGS_LIB_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

