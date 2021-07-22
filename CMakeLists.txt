# ==============================================================================
# PMFLib CMake File
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
PROJECT(PMFLib)

INCLUDE(CTest)
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# ==============================================================================
# conditional setup ------------------------------------------------------------
# ==============================================================================

SET(LIBS_SHARED             ON  CACHE BOOL "Should the shared version of pmflib library be built?")
SET(LIBS_STATIC             OFF CACHE BOOL "Should the static version of pmflib library be built?")

SET(PMFLIB_HAVE_XBPLIB      OFF CACHE BOOL "Is XBPLib available?")
SET(PMFLIB_NETWORK          ON  CACHE BOOL "Is network available")
SET(PMFLIB_COMPILE_FPMF     ON  CACHE BOOL "Should the fortran part of pmflib be built?")
SET(PMFLIB_HAVE_RPATH       ON  CACHE BOOL "Is rpath option supported by linker?")
SET(PMFLIB_USE_INST_RPATH   OFF CACHE BOOL "Should CMAKE_INSTALL_RPATH be used as the single RPATH?")
SET(PMFLIB_BUILD_UTILS      ON  CACHE BOOL "Should utilities be build?")
SET(PMFLIB_OPENMP           OFF CACHE BOOL "Should OpenMP be enabled for some PMFLib utilities?")

# drivers options --------------------------------
SET(PMFLIB_SANDER_DRV   ON   CACHE BOOL "Should the driver for sander be built?")
SET(PMFLIB_PMEMD_DRV    ON   CACHE BOOL "Should the driver for pmemd be built?")
SET(PMFLIB_CATS_DRV     ON   CACHE BOOL "Should the driver for cats be built?")
SET(PMFLIB_CPMD_DRV     ON   CACHE BOOL "Should the driver for cpmd be built?")
SET(PMFLIB_CP2K_DRV     ON   CACHE BOOL "Should the driver for cp2k be built?")
SET(PMFLIB_XDYNBP_DRV   OFF  CACHE BOOL "Should the driver for xdynbp be built?")
SET(PMFLIB_LIBATOMS_DRV OFF  CACHE BOOL "Should the driver for libatoms be built?")

# ==============================================================================
# project setup ----------------------------------------------------------------
# ==============================================================================

# define project version -------------------------
SET(PMFLIB_VERSION 6.0.x)

# what languages we use --------------------------
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)

# ==============================================================================
# external dependencies --------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    SET(DEVELOPMENT_ROOT "$ENV{DEVELOPMENT_ROOT}")
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

IF(NOT DEVELOPMENT_ROOT)
    MESSAGE(FATAL_ERROR "DEVELOPMENT_ROOT variable is not set!")
ENDIF(NOT DEVELOPMENT_ROOT)

# ==============================================================================
# external libraries -----------------------------------------------------------
# ==============================================================================
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${DEVELOPMENT_ROOT}/projects/cmake-modules/1.0/share/modules")

# boost =======================
FIND_PACKAGE(MyBoost REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})

# HIPOLY ====================
SET(HIPOLY_ROOT ${DEVELOPMENT_ROOT}/projects/hipoly/1.0)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/compress)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/error)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/event)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/file)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/list)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/option)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/string)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/stream)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/thread)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/time)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/xml)
LINK_DIRECTORIES(${HIPOLY_ROOT}/lib)
SET(HIPOLY_LIB_NAME hipoly)

# SCIMAFIC ====================
SET(SCIMAFIC_ROOT ${DEVELOPMENT_ROOT}/projects/scimafic/1.0)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/files)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/maths)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/chemistry)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/physics)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/fscimafic)
LINK_DIRECTORIES(${SCIMAFIC_ROOT}/lib)
SET(SCIMAFIC_CLIB_NAME cscimafic)
SET(SCIMAFIC_FLIB_NAME fscimafic)

# PRMFile ===================
SET(PRMFILE_ROOT ${DEVELOPMENT_ROOT}/projects/prmfile/1.0)
INCLUDE_DIRECTORIES(${PRMFILE_ROOT}/src/lib/cprmfile)
INCLUDE_DIRECTORIES(${PRMFILE_ROOT}/src/lib/fprmfile)
LINK_DIRECTORIES(${PRMFILE_ROOT}/lib)
SET(PRMFILE_CLIB_NAME cprmfile)
SET(PRMFILE_FLIB_NAME fprmfile)

# NETLib ====================
SET(NETLIB_ROOT ${DEVELOPMENT_ROOT}/projects/netlib/1.0)
INCLUDE_DIRECTORIES(${NETLIB_ROOT}/src/lib/netlib)
INCLUDE_DIRECTORIES(${NETLIB_ROOT}/src/lib/netlib/common)
INCLUDE_DIRECTORIES(${NETLIB_ROOT}/src/lib/netlib/client)
INCLUDE_DIRECTORIES(${NETLIB_ROOT}/src/lib/netlib/server)
INCLUDE_DIRECTORIES(${NETLIB_ROOT}/src/lib/netlib/extra/server)
INCLUDE_DIRECTORIES(${NETLIB_ROOT}/src/lib/netlib/extra/client)
LINK_DIRECTORIES(${NETLIB_ROOT}/lib)
SET(NETLIB_LIB_NAME netlib)

# setup BLAS/LAPACK either from MKL or from fallback netlib blas and lapack
FIND_PACKAGE(BlasLapack)

# runlux =====================
SET(RANLUX_ROOT ${DEVELOPMENT_ROOT}/thirdparty/ranlux/1.0)
LINK_DIRECTORIES(${RANLUX_ROOT}/lib)
SET(RANLUX_LIB_NAME franlux)

# CSPARSE ======================
SET(CSPARSE_ROOT ${DEVELOPMENT_ROOT}/thirdparty/csparse/2.2.1)
INCLUDE_DIRECTORIES(${CSPARSE_ROOT}/src/lib/csparse)
LINK_DIRECTORIES(${CSPARSE_ROOT}/lib)
SET(CSPARSE_LIB_NAME csparse)

# L-BFGS =====================
SET(LBFGS_ROOT ${DEVELOPMENT_ROOT}/thirdparty/lbfgs/1.0)
LINK_DIRECTORIES(${LBFGS_ROOT}/lib)
SET(LBFGS_LIB_NAME flbfgs)

IF(PMFLIB_HAVE_XBPLIB)
    # XBPLib ===================
    SET(XBPLIB_ROOT ${DEVELOPMENT_ROOT}/projects/xbplib/4.0)
    INCLUDE_DIRECTORIES(${XBPLIB_ROOT}/src/lib/xcore)
    INCLUDE_DIRECTORIES(${XBPLIB_ROOT}/src/lib/xfiles)
    INCLUDE_DIRECTORIES(${XBPLIB_ROOT}/src/lib/xene)
    LINK_DIRECTORIES(${XBPLIB_ROOT}/lib)
    SET(XBPLIB_XCORE_LIB_NAME xcore)
    SET(XBPLIB_XFILES_LIB_NAME xfiles)
    SET(XBPLIB_XENE_LIB_NAME xene)
ENDIF(PMFLIB_HAVE_XBPLIB)

# setup IntelCDK if available
INCLUDE(IntelCDK)
# setup MPI if available
INCLUDE(OpenMPI)

# if parallel environment is available
IF( NOT MPI_ENABLED )  # do not mix MPI and OpenMP
    IF( PMFLIB_OPENMP OR (MKL_MODE STREQUAL "parallel") OR (MKL_MODE STREQUAL "parallel/ilp64"))
        # setup OpenMP for cpmf if available
        FIND_PACKAGE(OpenMP)
        MESSAGE("-- OpenMP MODE: enabled")
    ELSE()
        MESSAGE("-- OpenMP MODE: disabled")
    ENDIF()
ELSE()
    MESSAGE("-- OpenMP MODE: disabled")
ENDIF()

# ==============================================================================
# output directories -----------------------------------------------------------
# ==============================================================================

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib/drivers)

# ==============================================================================
# rpaths -----------------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    # the RPATH to be used when installing
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/drivers")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# compiler setups --------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # C++ options
    SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas -Werror -pedantic-errors -Wundef -Wno-long-long")

    # get the name of fortran compiler
    get_filename_component(FCOMP_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)
    IF(MPI_ENABLED)
        EXECUTE_PROCESS(COMMAND mpifort --showme:command OUTPUT_VARIABLE FCOMP_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)
    ENDIF(MPI_ENABLED)

    # default compiler setups - gfortran
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O0")
    SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Wall -fcheck=all")

    # ifort
    IF(FCOMP_NAME STREQUAL ifort)
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
        SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -stand f03 -traceback -check all -warn all")
    ENDIF(FCOMP_NAME STREQUAL ifort)

    #SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_BUILD_TYPE Debug)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# project subdirectories  ------------------------------------------------------
# ==============================================================================

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(share)

