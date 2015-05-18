# ==============================================================================
# PMFLib CMake File
# ==============================================================================

INCLUDE_DIRECTORIES(lib/cpmf SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/abf SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/common SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/mtd SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/rst SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/methods/stm SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/abf SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/mtd SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/remd SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/common SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/client SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/network/stm SYSTEM)
INCLUDE_DIRECTORIES(lib/cpmf/mask SYSTEM)

INCLUDE_DIRECTORIES(lib/fpmf SYSTEM)
INCLUDE_DIRECTORIES(lib/fpmfdyn SYSTEM)
INCLUDE_DIRECTORIES(lib/fdrivers/common SYSTEM)

SET(CPMF_LIBRARIES
    cpmf ${NETLIB_LIB_NAME} ${PRMFILE_CLIB_NAME} ${CSPARSE_LIB_NAME}
         ${SCIMAFIC_CLIB_NAME} ${HIPOLY_LIB_NAME} ${SYSTEM_LIBS_CMAKE}
    )

# all libraries ----------------------------------
ADD_SUBDIRECTORY(lib)

# all binaries -----------------------------------
ADD_SUBDIRECTORY(bin)
