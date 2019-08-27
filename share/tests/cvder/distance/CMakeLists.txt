# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# distance CVs -----------------------------------------------------------------

ADD_TEST(CVDerRand01-DIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in dis.in)
ADD_TEST(CVDerRand02-DIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in dis.in)
ADD_TEST(CVDerRand03-DIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in dis.in)

ADD_TEST(CVDerRand01-ODIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in odis.in)
ADD_TEST(CVDerRand02-ODIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in odis.in)
ADD_TEST(CVDerRand03-ODIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in odis.in)

ADD_TEST(CVDerRand01-DD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in dd.in)
ADD_TEST(CVDerRand02-DD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in dd.in)
ADD_TEST(CVDerRand03-DD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in dd.in)

ADD_TEST(CVDerRand01-POS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in pos.in)
ADD_TEST(CVDerRand02-POS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in pos.in)
ADD_TEST(CVDerRand03-POS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in pos.in)

ADD_TEST(CVDerRand01-RAD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in rad.in)
ADD_TEST(CVDerRand02-RAD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in rad.in)
ADD_TEST(CVDerRand03-RAD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in rad.in)

ADD_TEST(CVDerRand01-OPOS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in opos.in)
ADD_TEST(CVDerRand02-OPOS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in opos.in)
ADD_TEST(CVDerRand03-OPOS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in opos.in)

ADD_TEST(CVDerRand01-ORAD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in orad.in)
ADD_TEST(CVDerRand02-ORAD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in orad.in)
ADD_TEST(CVDerRand03-ORAD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in orad.in)

ADD_TEST(CVDerRand01-MDIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in mdis.in)
ADD_TEST(CVDerRand02-MDIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in mdis.in)
ADD_TEST(CVDerRand03-MDIS ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in mdis.in)

ADD_TEST(CVDerRand01-MDISG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random01.in mdisg.in)
ADD_TEST(CVDerRand02-MDISG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random02.in mdisg.in)
ADD_TEST(CVDerRand03-MDISG ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random03.in mdisg.in)





