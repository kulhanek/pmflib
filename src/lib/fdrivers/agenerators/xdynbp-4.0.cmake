# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# overwrite previous files
file(WRITE ${CMAKE_BINARY_DIR}/include/xdynbp-4.0.inc "")

# create inline Makefile description for xdynbp --------------------------------

IF(PMFLIB_XBPLIB_AVAILABLE AND PMFLIB_XDYNBP_DRV)

    IF(PMFLIB_XBPLIB_AVAILABLE)
        SET(FPMF_CORE fpmf_egap)
    ELSE(PMFLIB_XBPLIB_AVAILABLE)
        SET(FPMF_CORE fpmf)
    ENDIF(PMFLIB_XBPLIB_AVAILABLE)

    IF(PMFLIB_NETWORK)
        SET(PMFLIB_LIBS "fpmfdrv_xdynbp ${FPMF_CORE} cpmf
                        ${CSPARSE_LIB_NAME}
                        ${SCIMAFIC_CLIB_NAME}
                        ${HIPOLY_LIB_NAME}
                        ${SYSTEM_LIBS_CMAKE}" )
        SET(PMFLIB_LIBS_PATH "${CSPARSE_ROOT}/lib ${SCIMAFIC_ROOT}/lib ${HIPOLY_ROOT}/lib")
    ELSE(PMFLIB_NETWORK)
        SET(PMFLIB_LIBS "fpmfdrv_xdynbp ${FPMF_CORE}" )
        SET(PMFLIB_LIBS_PATH "")
    ENDIF(PMFLIB_NETWORK)

    file(WRITE ${CMAKE_BINARY_DIR}/include/xdynbp-4.0.inc "LINK_DIRECTORIES(${PMFLIB_LIBS_PATH})\n")
    file(APPEND ${CMAKE_BINARY_DIR}/include/xdynbp-4.0.inc "SET(PMFLIB_LIBS ${PMFLIB_LIBS})\n")

ENDIF(PMFLIB_XBPLIB_AVAILABLE AND PMFLIB_XDYNBP_DRV)