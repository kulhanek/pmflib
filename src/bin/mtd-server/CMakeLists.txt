# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(MTD_SERVER_SRC
        MTDServer.cpp
        MTDServerHist.cpp
        MTDSrvOptions.cpp
        MTDProcessor.cpp
        MTDFactory.cpp
        MTDRegClient.cpp
        OpFlushServerData.cpp
        OpGetInitialData.cpp
        OpGetServerInfo.cpp
        OpRegisterClient.cpp
        OpUnregisterClient.cpp
        OpExchangeData.cpp
        OpGetData.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(mtd-server ${MTD_SERVER_SRC})
ADD_DEPENDENCIES(mtd-server cpmf_shared)

TARGET_LINK_LIBRARIES(mtd-server  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            mtd-server
        RUNTIME DESTINATION
            bin
        )

