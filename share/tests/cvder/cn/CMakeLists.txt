# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# Coordination Number CVs ------------------------------------------------------

ADD_TEST(CVDerRand01-CNFF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in cnff.in)
ADD_TEST(CVDerRand02-CNFF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in cnff.in)
ADD_TEST(CVDerRand03-CNFF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in cnff.in)

ADD_TEST(CVDerRand01-CNGFF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in cngff.in)
ADD_TEST(CVDerRand02-CNGFF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in cngff.in)
ADD_TEST(CVDerRand03-CNGFF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in cngff.in)

ADD_TEST(CVDerRand01-CNGFA ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in cngfa.in)
ADD_TEST(CVDerRand02-CNGFA ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in cngfa.in)
ADD_TEST(CVDerRand03-CNGFA ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in cngfa.in)

ADD_TEST(CVDerRand01-CNSW ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in cnsw.in)
ADD_TEST(CVDerRand02-CNSW ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in cnsw.in)
ADD_TEST(CVDerRand03-CNSW ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in cnsw.in)

ADD_TEST(CVDerRand01-CNGSW ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in cngsw.in)
ADD_TEST(CVDerRand02-CNGSW ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in cngsw.in)
ADD_TEST(CVDerRand03-CNGSW ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in cngsw.in)

ADD_TEST(CVDerRand01-CNRF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in cnrf.in)
ADD_TEST(CVDerRand02-CNRF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in cnrf.in)
ADD_TEST(CVDerRand03-CNRF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in cnrf.in)

ADD_TEST(CVDerRand01-CNGRF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in cngrf.in)
ADD_TEST(CVDerRand02-CNGRF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in cngrf.in)
ADD_TEST(CVDerRand03-CNGRF ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in cngrf.in)





