# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_ADMIN_SRC
        ABFAdmin.cpp
        ABFAdmOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-admin ${ABF_ADMIN_SRC})

TARGET_LINK_LIBRARIES(abf-admin  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-admin
        RUNTIME DESTINATION
            bin
        )
