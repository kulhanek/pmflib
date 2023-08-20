# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(OPTGPRHPRMS_MULTI_SRC
        OptGPRHyprmsMulti.cpp
        OptGPRHyprmsMultiOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _gpr-opthyprms-multi)
ELSE()
    SET(TARGET_NAME gpr-opthyprms-multi)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${OPTGPRHPRMS_MULTI_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${LBFGS_LIB_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

