# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# WORM CVs ---------------------------------------------------------------------

ADD_TEST(CVDerRand01-WORMPOS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in wormpos.in)
ADD_TEST(CVDerRand02-WORMPOS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in wormpos.in)
ADD_TEST(CVDerRand03-WORMPOS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in wormpos.in)
ADD_TEST(CVDerWORMTRAJ-WORMPOS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        worm-I.in wormpos-I.in)

ADD_TEST(CVDerRand01-WORMANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in wormang.in)
ADD_TEST(CVDerRand02-WORMANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in wormang.in)
ADD_TEST(CVDerRand03-WORMANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in wormang.in)
ADD_TEST(CVDerWORMTRAJ-WORMANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        worm-I.in wormang-I.in)




