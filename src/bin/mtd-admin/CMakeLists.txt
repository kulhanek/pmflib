# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(MTD_ADMIN_SRC
        MTDAdmin.cpp
        MTDAdmOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(mtd-admin ${MTD_ADMIN_SRC})

TARGET_LINK_LIBRARIES(mtd-admin  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            mtd-admin
        RUNTIME DESTINATION
            bin
        )

