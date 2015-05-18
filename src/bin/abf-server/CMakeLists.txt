# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABF_SERVER_SRC
        ABFServer.cpp
        ABFSrvOptions.cpp
        ABFProcessor.cpp
        ABFFactory.cpp
        ABFServerAccu.cpp
        OpFlushServerData.cpp
        OpGetInitialData.cpp
        OpGetServerInfo.cpp
        OpRegisterClient.cpp
        OpUnregisterClient.cpp
        OpExchangeData.cpp
        OpGetData.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abf-server ${ABF_SERVER_SRC})

TARGET_LINK_LIBRARIES(abf-server  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abf-server
        RUNTIME DESTINATION
            bin
        )

