# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(MWA_SERVER_SRC
        MWAServer.cpp
        MWASrvOptions.cpp
        MWAProcessor.cpp
        MWAFactory.cpp
        MWAServerAccu.cpp
        OpFlushServerData.cpp
        OpGetInitialData.cpp
        OpGetServerInfo.cpp
        OpRegisterClient.cpp
        OpUnregisterClient.cpp
        OpExchangeData.cpp
        OpGetData.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(mwa-server ${MWA_SERVER_SRC})
ADD_DEPENDENCIES(mwa-server cpmf_shared)

TARGET_LINK_LIBRARIES(mwa-server  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            mwa-server
        RUNTIME DESTINATION
            bin
        )

