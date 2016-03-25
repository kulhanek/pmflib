# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# nucleic acid -----------------------------------------------------------------

ADD_TEST(CVDerRand01-NABO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in nabo.in)
ADD_TEST(CVDerRand02-NABO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in nabo.in)
ADD_TEST(CVDerRand03-NABO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in nabo.in)







