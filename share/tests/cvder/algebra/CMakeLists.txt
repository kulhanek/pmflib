# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# nucleic acid -----------------------------------------------------------------

ADD_TEST(CVDerRand04-ADD ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random04.in add.in)

ADD_TEST(CVDerRand04-SUB ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random04.in sub.in)

ADD_TEST(CVDerRand04-MUL ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random04.in mul.in)

ADD_TEST(CVDerRand04-DIV ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random04.in div.in)

ADD_TEST(CVDerRand04-FSWITCH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random04.in fswitch.in)

ADD_TEST(CVDerRand04-RSWITCH ${EXECUTABLE_OUTPUT_PATH}/pmf-test-cvs
                        ../random04.in rswitch.in)







