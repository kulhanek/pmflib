# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(RST_WHAM_SRC
        bootstrap.c
        file_read.c
        histogram.c
        locate.c
        ran2.c
        wham.c
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(rst-wham ${RST_WHAM_SRC})

TARGET_LINK_LIBRARIES(rst-wham m)

INSTALL(TARGETS
            rst-wham
        RUNTIME DESTINATION
            bin
        )

