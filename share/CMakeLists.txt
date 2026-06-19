# ==============================================================================
# PMFLib CMake File
# ==============================================================================

ADD_SUBDIRECTORY(tests)

# OpenMP wrappers
IF(OpenMP_FOUND)
    INSTALL(FILES
                bin/pmf-integrate
                bin/aus-integrate
                bin/pmf-energy
                bin/gpr-opthyprms
            DESTINATION
                bin
            PERMISSIONS
                OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
            )
ENDIF()

# Python stuff
INSTALL(FILES
            python/analyse-act.py
            python/analyse-1D-surface.py
            python/analyse-2D-surface.py
            python/stm-path-2D-surface.py
        DESTINATION
            share/python
        PERMISSIONS
            OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
        )

INSTALL(FILES
            scripts/analyse-act
            scripts/analyse-1D-surface
            scripts/analyse-2D-surface
            scripts/stm-path-2D-surface
        DESTINATION
            share/scripts
        PERMISSIONS
            OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
        )