# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# shape CVs --------------------------------------------------------------------

ADD_TEST(CVDerRand01-PLANE ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in plane.in)
ADD_TEST(CVDerRand02-PLANE ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in plane.in)
ADD_TEST(CVDerRand03-PLANE ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in plane.in)

ADD_TEST(CVDerRand01-RGYR ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in rgyr.in)
ADD_TEST(CVDerRand02-RGYR ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in rgyr.in)
ADD_TEST(CVDerRand03-RGYR ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in rgyr.in)

ADD_TEST(CVDerRand01-ACYL ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in acyl.in)
ADD_TEST(CVDerRand02-ACYL ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in acyl.in)
ADD_TEST(CVDerRand03-ACYL ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in acyl.in)

ADD_TEST(CVDerRand01-ASPH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in asph.in)
ADD_TEST(CVDerRand02-ASPH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in asph.in)
ADD_TEST(CVDerRand03-ASPH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in asph.in)

ADD_TEST(CVDerRand01-SANIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in sanis.in)
ADD_TEST(CVDerRand02-SANIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in sanis.in)
ADD_TEST(CVDerRand03-SANIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in sanis.in)

ADD_TEST(CVDerRand01-PMOGT ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in pmogt.in)
ADD_TEST(CVDerRand02-PMOGT ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in pmogt.in)
ADD_TEST(CVDerRand03-PMOGT ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in pmogt.in)

ADD_TEST(CVDerRand01-PMGTD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in pmgtd.in)
ADD_TEST(CVDerRand02-PMGTD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in pmgtd.in)
ADD_TEST(CVDerRand03-PMGTD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in pmgtd.in)

ADD_TEST(CVDerRand01-RMSDT ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in rmsdt.in)
ADD_TEST(CVDerRand02-RMSDT ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in rmsdt.in)
ADD_TEST(CVDerRand03-RMSDT ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in rmsdt.in)

ADD_TEST(CVDerRand01-EVEC ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in evec.in)
ADD_TEST(CVDerRand02-EVEC ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in evec.in)
ADD_TEST(CVDerRand03-EVEC ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in evec.in)

