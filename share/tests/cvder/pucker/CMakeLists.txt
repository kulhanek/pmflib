# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# ring puckering ---------------------------------------------------------------

ADD_TEST(CVDerRand01-PUCK6Q ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in puck6q.in)
ADD_TEST(CVDerRand02-PUCK6Q ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in puck6q.in)
ADD_TEST(CVDerRand03-PUCK6Q ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in puck6q.in)

ADD_TEST(CVDerRand01-PUCK6T ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in puck6t.in)
ADD_TEST(CVDerRand02-PUCK6T ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in puck6t.in)
ADD_TEST(CVDerRand03-PUCK6T ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in puck6t.in)

ADD_TEST(CVDerRand01-PUCK6P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in puck6p.in)
ADD_TEST(CVDerRand02-PUCK6P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in puck6p.in)
ADD_TEST(CVDerRand03-PUCK6P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in puck6p.in)

ADD_TEST(CVDerRand01-PUCK5Q ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in puck5q.in)
ADD_TEST(CVDerRand02-PUCK5Q ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in puck5q.in)
ADD_TEST(CVDerRand03-PUCK5Q ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in puck5q.in)

ADD_TEST(CVDerRand01-PUCK5P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in puck5p.in)
ADD_TEST(CVDerRand02-PUCK5P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in puck5p.in)
ADD_TEST(CVDerRand03-PUCK5P ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in puck5p.in)






