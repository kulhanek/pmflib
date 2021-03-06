# ==============================================================================
# PMFLib CMake File
# ==============================================================================

SET(CPMF_METHODS_SRC
        PMFMainHeader.cpp
        methods/common/ColVariable.cpp
        methods/common/EnergySurface.cpp
        methods/common/ESPrinter.cpp
        methods/abp/ABPAccumulator.cpp
        methods/abf/ABFAccumulator.cpp
        methods/abf/ABFIntegratorRFD.cpp
        methods/abf/ABFIntegratorRFD2.cpp
        methods/abf/ABFIntegratorRBF.cpp
        methods/abf/ABFIntegratorGPR.cpp
        methods/abf/ABFEnthalpyGPR.cpp
        methods/mtd/MTDHistory.cpp
        methods/mtd/MTDBuffer.cpp
        methods/rst/RSTAccumulator.cpp
        methods/stm/spline.cpp
        methods/stm/Bead.cpp
        methods/stm/CVSpline.cpp
        methods/stm/BeadList.cpp
        )

SET(CPMF_MASK_SRC
        mask/PMFResidue.cpp
        mask/PMFAtom.cpp
        mask/PMFTopology.cpp
        mask/maskparser/MaskParser.cpp
        mask/maskparser/Mask.yy.c
        mask/maskparser/Mask.tab.c
        mask/PMFMaskSelection.cpp
        mask/PMFMask.cpp
        mask/EntriesPMFMask.cpp
        )

SET(CPMF_NETWORK_SRC
        network/PMFOperation.cpp
        network/common/Snapshot.cpp
        network/client/AdminClient.cpp
        network/client/COpFlushServerData.cpp
        network/client/COpGetServerInfo.cpp
    # ABF ----------------------------------------
        network/abf/ABFClient.cpp
        network/abf/EntriesABFClient.cpp
    # ABP ----------------------------------------
        network/abp/ABPClient.cpp
        network/abp/EntriesABPClient.cpp
    # MTD ----------------------------------------
        network/mtd/MTDClient.cpp
        network/mtd/EntriesMTDClient.cpp
    # REMD ---------------------------------------
        network/remd/REMDClient.cpp
        network/remd/EntriesREMDClient.cpp
    # STM ----------------------------------------
        network/stm/STMClient.cpp
        network/stm/EntriesSTMClient.cpp
        )

SET(CPMF_SRC
        ${CPMF_METHODS_SRC} ${CPMF_MASK_SRC} ${CPMF_NETWORK_SRC}
        )

IF(OpenMP_FOUND)
    ADD_DEFINITIONS(${OpenMP_CXX_FLAGS})
ENDIF(OpenMP_FOUND)

#-------------------------------------------------------------------------------

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(cpmf_static STATIC ${CPMF_SRC})
    SET_TARGET_PROPERTIES(cpmf_static PROPERTIES
                            OUTPUT_NAME cpmf
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_DEFINITIONS(-DPMFLIB_BUILDING_DLL)
    ADD_LIBRARY(cpmf_shared SHARED ${CPMF_SRC})
IF(OpenMP_FOUND)
    TARGET_LINK_LIBRARIES(cpmf_shared OpenMP::OpenMP_CXX ${CPMF_LIBRARIES_DEPS})
ELSE(OpenMP_FOUND)
    TARGET_LINK_LIBRARIES(cpmf_shared ${CPMF_LIBRARIES_DEPS})
ENDIF(OpenMP_FOUND)
    SET_TARGET_PROPERTIES(cpmf_shared PROPERTIES
                            OUTPUT_NAME cpmf
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION})
    INSTALL(TARGETS
                cpmf_shared
            DESTINATION
                lib
            )
ENDIF(LIBS_SHARED)

#-------------------------------------------------------------------------------
