# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# test ABF accumulator reader --------------------------------------------------

CONFIGURE_FILE("_abf.rst.v3" "_abf.rst" COPYONLY)

ADD_TEST(ABFReader-ENE-v3 ${EXECUTABLE_OUTPUT_PATH}/abf-integrate _abf.rst ene)
ADD_TEST(ABFReader-RST-v3 ${EXECUTABLE_OUTPUT_PATH}/pmf-dyn restartabf1.in)


