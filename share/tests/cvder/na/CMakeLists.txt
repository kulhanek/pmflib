# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# nucleic acid -----------------------------------------------------------------

ADD_TEST(CVDerRand01-NAPBO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in napbo.in)
ADD_TEST(CVDerRand02-NAPBO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in napbo.in)
ADD_TEST(CVDerRand03-NAPBO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in napbo.in)

ADD_TEST(CVDerRand01-NASBO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in nasbo.in)
ADD_TEST(CVDerRand02-NASBO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in nasbo.in)
ADD_TEST(CVDerRand03-NASBO ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in nasbo.in)

ADD_TEST(CVDerRand01-NABEND ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in nabend.in)
ADD_TEST(CVDerRand02-NABEND ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in nabend.in)
ADD_TEST(CVDerRand03-NABEND ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in nabend.in)

ADD_TEST(CVDerRand01-NASBPP ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in nasbpp.in)
ADD_TEST(CVDerRand02-NASBPP ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in nasbpp.in)
ADD_TEST(CVDerRand03-NASBPP ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in nasbpp.in)







