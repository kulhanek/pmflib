# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ENTROPYDER_SRC
        EntropyDer.cpp
        EntropyDerOptions.cpp
        )

# final build ------------------------------------------------------------------

IF(OpenMP_FOUND)
    SET(TARGET_NAME _pmf-entropy-der)
ELSE()
    SET(TARGET_NAME pmf-entropy-der)
ENDIF()

ADD_EXECUTABLE(${TARGET_NAME} ${ENTROPYDER_SRC})
ADD_DEPENDENCIES(${TARGET_NAME} cpmf_shared)

TARGET_LINK_LIBRARIES(${TARGET_NAME} ${CPMF_LIBRARIES})

INSTALL(TARGETS
            ${TARGET_NAME}
        RUNTIME DESTINATION
            bin
        )

