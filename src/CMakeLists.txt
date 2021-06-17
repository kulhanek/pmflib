# ==============================================================================
# PMFLib CMake File
# ==============================================================================

INCLUDE_DIRECTORIES(lib/cpmf SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/abf SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/abp SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/common SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/mtd SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/rst SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/stm SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/abf SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/abp SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/mtd SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/common SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/client SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/stm SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/remd SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/mask SYSTEM)

INCLUDE_DIRECTORIES(lib/fpmf SYSTEM)
INCLUDE_DIRECTORIES(lib/fpmfdyn SYSTEM)

# ------------------------------------------------
SET(CPMF_LIBRARIES_DEPS
        ${NETLIB_LIB_NAME}
        ${PRMFILE_CLIB_NAME}
        ${CSPARSE_LIB_NAME}
        ${SCIMAFIC_CLIB_NAME}
        ${HIPOLY_LIB_NAME}
        ${SYSTEM_LIBS_CMAKE}
    )
SET(CPMF_LIBRARIES
    cpmf
        ${CPMF_LIBRARIES_DEPS}
    )

SET(FPMF_LIBRARIES_DEPS
        ${XBPLIB_XCORE_LIB_NAME}
        ${XBPLIB_XFILES_LIB_NAME}
        ${XBPLIB_XENE_LIB_NAME}
        ${PRMFILE_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${RANLUX_LIB_NAME}
    )
SET(FPMF_LIBRARIES
    fpmf
        ${FPMF_LIBRARIES_DEPS}
    )

# all libraries ----------------------------------
ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(test)

# all binaries -----------------------------------
IF(PMFLIB_BUILD_UTILS)
    ADD_SUBDIRECTORY(bin)
ENDIF(PMFLIB_BUILD_UTILS)
