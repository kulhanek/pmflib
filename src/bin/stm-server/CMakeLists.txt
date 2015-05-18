# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(STRING_SERVER_SRC
        StringServer.cpp
        StringSrvOptions.cpp
        StringProcessor.cpp
        StringFactory.cpp
        OpFlushServerData.cpp
        OpGetServerInfo.cpp
        OpRegisterClient.cpp
        OpUnregisterClient.cpp
        OpExchangeData.cpp
        OpShutdownTerminateServer.cpp
        OpGetData.cpp
        Launcher.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(stm-server ${STRING_SERVER_SRC})

TARGET_LINK_LIBRARIES(stm-server ${CPMF_LIBRARIES})

INSTALL(TARGETS
            stm-server
        DESTINATION
            bin
        )

