# ==============================================================================
# PMFLib CMake File
# ==============================================================================

ADD_SUBDIRECTORY(tests)

# OpenMP wrappers
IF(OpenMP_FOUND)
    INSTALL(FILES
                bin/abf-integrate
                bin/abf-optgprhyprms
            DESTINATION
                bin
            PERMISSIONS
                OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
            )
ENDIF()
