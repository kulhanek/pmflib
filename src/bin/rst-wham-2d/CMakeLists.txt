# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(RST_WHAM_2D_SRC
        bootstrap.c
        file_read.c
        histogram.c
        locate.c
        ran2.c
        wham-2d.c
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(rst-wham-2d ${RST_WHAM_2D_SRC})

TARGET_LINK_LIBRARIES(rst-wham-2d m)

INSTALL(TARGETS
            rst-wham-2d
        RUNTIME DESTINATION
            bin
        )

