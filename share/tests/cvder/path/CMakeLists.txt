# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# paths -----------------------------------------------------------------

ADD_TEST(CVDerRand01-PTPATHS ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random01.in ptpaths.in)
ADD_TEST(CVDerRand02-PTPATHS ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random02.in ptpaths.in)
ADD_TEST(CVDerRand03-PTPATHS ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random03.in ptpaths.in)

ADD_TEST(CVDerRand01-PTPATHZ ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random01.in ptpathz.in)
ADD_TEST(CVDerRand02-PTPATHZ ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random02.in ptpathz.in)
ADD_TEST(CVDerRand03-PTPATHZ ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random03.in ptpathz.in)

ADD_TEST(CVDerRand01-PTPATHS2P ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random01.in ptpaths2p.in)
ADD_TEST(CVDerRand02-PTPATHS2P ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random02.in ptpaths2p.in)
ADD_TEST(CVDerRand03-PTPATHS2P ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random03.in ptpaths2p.in)


ADD_TEST(CVDerRand01-PTPATHS2PN ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random01.in ptpaths2pn.in)
ADD_TEST(CVDerRand02-PTPATHS2PN ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random02.in ptpaths2pn.in)
ADD_TEST(CVDerRand03-PTPATHS2PN ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random03.in ptpaths2pn.in)



ADD_TEST(CVDerRand01-PATHS ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random01.in paths.in path.in)
ADD_TEST(CVDerRand02-PATHS ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random02.in paths.in path.in)
ADD_TEST(CVDerRand03-PATHS ${EXECUTABLE_OUTPUT_PATH}/tests/test-cvs
                        ../random03.in paths.in path.in)






