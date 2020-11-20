# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABP_ADMIN_SRC
        ABPAdmin.cpp
        ABPAdmOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abp-admin ${ABP_ADMIN_SRC})
ADD_DEPENDENCIES(abp-admin cpmf_shared)

TARGET_LINK_LIBRARIES(abp-admin  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abp-admin
        RUNTIME DESTINATION
            bin
        )
