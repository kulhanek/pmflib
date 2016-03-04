# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# angle CVs --------------------------------------------------------------------

ADD_TEST(CVDerRand01-ANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in ang.in)
ADD_TEST(CVDerRand02-ANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in ang.in)
ADD_TEST(CVDerRand03-ANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in ang.in)

ADD_TEST(CVDerRand01-ANG2 ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in ang2.in)
ADD_TEST(CVDerRand02-ANG2 ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in ang2.in)
ADD_TEST(CVDerRand03-ANG2 ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in ang2.in)

ADD_TEST(CVDerRand01-DIH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in dih.in)
ADD_TEST(CVDerRand02-DIH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in dih.in)
ADD_TEST(CVDerRand03-DIH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in dih.in)

ADD_TEST(CVDerRand01-PANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in pang.in)
ADD_TEST(CVDerRand02-PANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in pang.in)
ADD_TEST(CVDerRand03-PANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in pang.in)

ADD_TEST(CVDerRand01-PVANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in pvang.in)
ADD_TEST(CVDerRand02-PVANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in pvang.in)
ADD_TEST(CVDerRand03-PVANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in pvang.in)

ADD_TEST(CVDerRand01-AXANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in axang.in)
ADD_TEST(CVDerRand02-AXANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in axang.in)
ADD_TEST(CVDerRand03-AXANG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in axang.in)

