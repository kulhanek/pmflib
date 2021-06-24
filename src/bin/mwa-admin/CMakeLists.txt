# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(MWA_ADMIN_SRC
        MWAAdmin.cpp
        MWAAdmOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(mwa-admin ${MWA_ADMIN_SRC})
ADD_DEPENDENCIES(mwa-admin  cpmf_shared)

TARGET_LINK_LIBRARIES(mwa-admin   ${CPMF_LIBRARIES})

INSTALL(TARGETS
            mwa-admin
        RUNTIME DESTINATION
            bin
        )
