# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# paths -----------------------------------------------------------------

ADD_TEST(CVDerRand01-PATHS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in paths.in)
ADD_TEST(CVDerRand02-PATHS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in paths.in)
ADD_TEST(CVDerRand03-PATHS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in paths.in)

ADD_TEST(CVDerRand01-PATHZ ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in pathz.in)
ADD_TEST(CVDerRand02-PATHZ ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in pathz.in)
ADD_TEST(CVDerRand03-PATHZ ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in pathz.in)

ADD_TEST(CVDerRand01-PATHS2P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in paths2p.in)
ADD_TEST(CVDerRand02-PATHS2P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in paths2p.in)
ADD_TEST(CVDerRand03-PATHS2P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in paths2p.in)










