# ==============================================================================
# PMFLib CMake File
# ==============================================================================

file(WRITE ${CMAKE_BINARY_DIR}/include/amber-9.0.inc "")

# create inline Makefile description for amber ---------------------------------
IF(PMFLIB_SANDER_DRV OR PMFLIB_PMEMD_DRV OR PMFLIB_SANDER_XDYNBP_DRV)

    # Network ========================================
    IF(PMFLIB_NETWORK)

        SET(PMFLIB_CORE_LIBS
"        -L${CMAKE_BINARY_DIR}/lib -L${CMAKE_BINARY_DIR}/lib/drivers \\
        -l${NETLIB_LIB_NAME} -L${NETLIB_ROOT}/lib \\
        -l${PRMFILE_CLIB_NAME} -l${PRMFILE_FLIB_NAME} -L${PRMFILE_ROOT}/lib \\
        -l${LAPACK_LIB_NAME} -L${LAPACK_ROOT}/lib \\
        -l${BLAS_LIB_NAME} -L${BLAS_ROOT}/lib \\
        -lcpmf -L${CMAKE_BINARY_DIR}/lib \\
        -l${CSPARSE_LIB_NAME} -L${CSPARSE_ROOT}/lib \\
        -l${SCIMAFIC_FLIB_NAME} -L${SCIMAFIC_ROOT}/lib \\
        -l${SCIMAFIC_CLIB_NAME} -L${SCIMAFIC_ROOT}/lib \\
        -l${HIPOLY_LIB_NAME} -L${HIPOLY_ROOT}/lib \\
        ${SYSTEM_LIBS}")

        IF(PMFLIB_HAVE_RPATH)
            IF(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_INSTALL_RPATH}")
            ELSE(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_BINARY_DIR}/lib:${CMAKE_BINARY_DIR}/lib/drivers:${NETLIB_ROOT}/lib:${PRMFILE_ROOT}/lib:${LAPACK_ROOT}/lib:${BLAS_ROOT}/lib:${CMAKE_BINARY_DIR}/lib:${CSPARSE_ROOT}/lib:${SCIMAFIC_ROOT}/lib:${HIPOLY_ROOT}/lib")
            ENDIF(PMFLIB_USE_INST_RPATH)
        ENDIF(PMFLIB_HAVE_RPATH)

    ELSE(PMFLIB_NETWORK)

        SET(PMFLIB_CORE_LIBS
"        -L${CMAKE_BINARY_DIR}/lib -L${CMAKE_BINARY_DIR}/lib/drivers \\
        -l${PRMFILE_FLIB_NAME} -L${PRMFILE_ROOT}/lib \\
        -l${LAPACK_LIB_NAME} -L${LAPACK_ROOT}/lib \\
        -l${BLAS_LIB_NAME} -L${BLAS_ROOT}/lib" )

        IF(PMFLIB_HAVE_RPATH)
            IF(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_INSTALL_RPATH}")
            ELSE(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_BINARY_DIR}/lib:${CMAKE_BINARY_DIR}/lib/drivers:${PRMFILE_ROOT}/lib:${LAPACK_ROOT}/lib:${BLAS_ROOT}/lib")
            ENDIF(PMFLIB_USE_INST_RPATH)
        ENDIF(PMFLIB_HAVE_RPATH)

    ENDIF(PMFLIB_NETWORK)
    # Network ========================================

    # XBPLib ========================================
    IF(PMFLIB_XBPLIB_AVAILABLE)
        SET(PMFLIB_SANDER_XDYNBP_LIBS
"        -lfpmfdrv_sander_xdynbp -lfpmf_egap \\
        -l${RANLUX_LIB_NAME} -L${RANLUX_ROOT}/lib  \\
        -l${XBPLIB_XENE_LIB_NAME} -l${XBPLIB_XFILES_LIB_NAME} \\
        -l${XBPLIB_XCORE_LIB_NAME} -L${XBPLIB_ROOT}/lib \\
${PMFLIB_CORE_LIBS}" )

        IF(PMFLIB_HAVE_RPATH)
            IF(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_XDYNBP_RPATHS "${PMFLIB_RPATHS}")
            ELSE(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_XDYNBP_RPATHS "${PMFLIB_RPATHS}:${XBPLIB_ROOT}/lib:${RANLUX_ROOT}/lib")
            ENDIF(PMFLIB_USE_INST_RPATH)
            SET(PMFLIB_SANDER_XDYNBP_LIBS
"${PMFLIB_SANDER_XDYNBP_LIBS} \\
        -Wl,-rpath=${PMFLIB_XDYNBP_RPATHS}" )
        ENDIF(PMFLIB_HAVE_RPATH)

    ENDIF(PMFLIB_XBPLIB_AVAILABLE)
    # XBPLib ========================================

    SET(PMFLIB_SANDER_LIBS
"        -lfpmfdrv_sander -lfpmf \\
${PMFLIB_CORE_LIBS}" )

    SET(PMFLIB_PMEMD_LIBS
"        -lfpmfdrv_pmemd -lfpmf \\
${PMFLIB_CORE_LIBS}" )

    IF(PMFLIB_HAVE_RPATH)
        SET(PMFLIB_SANDER_LIBS
"${PMFLIB_SANDER_LIBS} \\
        -Wl,-rpath=${PMFLIB_RPATHS}" )
        SET(PMFLIB_PMEMD_LIBS
"${PMFLIB_PMEMD_LIBS} \\
        -Wl,-rpath=${PMFLIB_RPATHS}" )
    ENDIF(PMFLIB_HAVE_RPATH)

    # MPI ========================================
    IF(PMFLIB_MPI)
        SET(PMFLIB_SANDER_MPI_LIBS
"        -lfpmfdrv_sander_mpi -lfpmf -lfpmf_mpi \\
${PMFLIB_CORE_LIBS}" )

        SET(PMFLIB_PMEMD_MPI_LIBS
"        -lfpmfdrv_pmemd_mpi -lfpmf -lfpmf_mpi \\
${PMFLIB_CORE_LIBS}" )

        IF(PMFLIB_HAVE_RPATH)
            SET(PMFLIB_SANDER_MPI_LIBS
"${PMFLIB_SANDER_MPI_LIBS} \\
        -Wl,-rpath=${PMFLIB_RPATHS}" )
            SET(PMFLIB_PMEMD_MPI_LIBS
"${PMFLIB_PMEMD_MPI_LIBS} \\
        -Wl,-rpath=${PMFLIB_RPATHS}" )

        ENDIF(PMFLIB_HAVE_RPATH)
    ENDIF(PMFLIB_MPI)
    # MPI ========================================

    file(WRITE ${CMAKE_BINARY_DIR}/include/amber-9.0.inc "PMFLIB_SANDER_EXT=\\\n${PMFLIB_SANDER_LIBS}\n\n")
    file(APPEND ${CMAKE_BINARY_DIR}/include/amber-9.0.inc "PMFLIB_SANDER_MPI_EXT=\\\n${PMFLIB_SANDER_MPI_LIBS}\n\n")
    file(APPEND ${CMAKE_BINARY_DIR}/include/amber-9.0.inc "PMFLIB_SANDER_XDYNBP_EXT=\\\n${PMFLIB_SANDER_XDYNBP_LIBS}\n\n")
    file(APPEND ${CMAKE_BINARY_DIR}/include/amber-9.0.inc "PMFLIB_PMEMD_EXT=\\\n${PMFLIB_PMEMD_LIBS}\n\n")
    file(APPEND ${CMAKE_BINARY_DIR}/include/amber-9.0.inc "PMFLIB_PMEMD_MPI_EXT=\\\n${PMFLIB_PMEMD_MPI_LIBS}\n\n")

ENDIF(PMFLIB_SANDER_DRV OR PMFLIB_PMEMD_DRV OR PMFLIB_SANDER_XDYNBP_DRV)
