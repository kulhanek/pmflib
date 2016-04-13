# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# nucleic acid -----------------------------------------------------------------

ADD_TEST(CVDerRand01-ADD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in add.in)
ADD_TEST(CVDerRand02-ADD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in add.in)
ADD_TEST(CVDerRand03-ADD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in add.in)







