# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(STRING_ADMIN_SRC
        StringAdmin.cpp
        StringAdmOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(stm-admin ${STRING_ADMIN_SRC})
ADD_DEPENDENCIES(stm-admin cpmf_shared)

TARGET_LINK_LIBRARIES(stm-admin  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            stm-admin
        DESTINATION
            bin
        )
