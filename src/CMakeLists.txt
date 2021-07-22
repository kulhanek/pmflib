# ==============================================================================
# PMFLib CMake File
# ==============================================================================

INCLUDE_DIRECTORIES(lib/cpmf)
INCLUDE_DIRECTORIES(lib/cpmf/common)
INCLUDE_DIRECTORIES(lib/cpmf/common/integrators)
INCLUDE_DIRECTORIES(lib/cpmf/common/smoothers)
INCLUDE_DIRECTORIES(lib/cpmf/methods/abf)
INCLUDE_DIRECTORIES(lib/cpmf/methods/abp)
INCLUDE_DIRECTORIES(lib/cpmf/methods/mtd)
INCLUDE_DIRECTORIES(lib/cpmf/methods/rst)
INCLUDE_DIRECTORIES(lib/cpmf/methods/stm)
INCLUDE_DIRECTORIES(lib/cpmf/methods/cst)
INCLUDE_DIRECTORIES(lib/cpmf/methods/pmf)
INCLUDE_DIRECTORIES(lib/cpmf/network)
INCLUDE_DIRECTORIES(lib/cpmf/network/abf)
INCLUDE_DIRECTORIES(lib/cpmf/network/abp)
INCLUDE_DIRECTORIES(lib/cpmf/network/mtd)
INCLUDE_DIRECTORIES(lib/cpmf/network/stm)
INCLUDE_DIRECTORIES(lib/cpmf/mask)

INCLUDE_DIRECTORIES(lib/fpmf)
INCLUDE_DIRECTORIES(lib/fpmfdyn)

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
