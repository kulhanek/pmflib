# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(ABP_SERVER_SRC
        ABPServer.cpp
        ABPSrvOptions.cpp
        ABPProcessor.cpp
        ABPFactory.cpp
        ABPServerAccu.cpp
        OpFlushServerData.cpp
        OpGetInitialData.cpp
        OpGetServerInfo.cpp
        OpRegisterClient.cpp
        OpUnregisterClient.cpp
        OpExchangeData.cpp
        OpGetData.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(abp-server ${ABP_SERVER_SRC})
ADD_DEPENDENCIES(abp-server cpmf_shared)

TARGET_LINK_LIBRARIES(abp-server  ${CPMF_LIBRARIES})

INSTALL(TARGETS
            abp-server
        RUNTIME DESTINATION
            bin
        )

