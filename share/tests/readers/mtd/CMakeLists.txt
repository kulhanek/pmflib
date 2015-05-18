# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# test MTD restart reader ------------------------------------------------------

ADD_TEST(MTDReader-ENE ${EXECUTABLE_OUTPUT_PATH}/mtd-energy _mtd.rst ene)
ADD_TEST(MTDReader-RST ${EXECUTABLE_OUTPUT_PATH}/pmf-dyn restartmtd.in)

