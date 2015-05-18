# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# test ABF accumulator reader --------------------------------------------------

ADD_TEST(ABFReader-ENE  ${EXECUTABLE_OUTPUT_PATH}/abf-integrate _abf.rst ene)

ADD_TEST(ABFReader-RST1 ${EXECUTABLE_OUTPUT_PATH}/pmf-dyn restartabf1.in)

ADD_TEST(ABFReader-RST2 ${EXECUTABLE_OUTPUT_PATH}/pmf-dyn restartabf2.in)
SET_PROPERTY(TEST ABFReader-RST2 PROPERTY WILL_FAIL "true")

ADD_TEST(ABFReader-RST3 ${EXECUTABLE_OUTPUT_PATH}/pmf-dyn restartabf3.in)

