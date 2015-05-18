# ==============================================================================
# PMFLib CMake File
# ==============================================================================

file(WRITE ${CMAKE_BINARY_DIR}/include/libatoms-dev.inc "")

# create inline Makefile description for amber ---------------------------------
IF(PMFLIB_LIBATOMS_DRV OR PMFLIB_LIBATOMS_XDYNBP_DRV)

    # Network ========================================
    IF(PMFLIB_NETWORK)

        SET(PMFLIB_CORE_LIBS
"        -L${CMAKE_BINARY_DIR}/lib \\
        -l${PRMFILE_FLIB_NAME} -L${PRMFILE_ROOT}/lib \\
        -l${LAPACK_LIB_NAME} -L${LAPACK_ROOT}/lib \\
        -l${BLAS_LIB_NAME} -L${BLAS_ROOT}/lib \\
        -lcpmf -L${CMAKE_BINARY_DIR}/lib \\
        -l${CSPARSE_LIB_NAME} -L${CSPARSE_ROOT}/lib \\
        -l${SCIMAFIC_CLIB_NAME} -L${SCIMAFIC_ROOT}/lib \\
        -l${HIPOLY_LIB_NAME} -L${HIPOLY_ROOT}/lib \\
        ${SYSTEM_LIBS}")

        IF(PMFLIB_HAVE_RPATH)
            IF(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_INSTALL_RPATH}")
            ELSE(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_BINARY_DIR}/lib::${PRMFILE_ROOT}/lib:${LAPACK_ROOT}/lib:${BLAS_ROOT}/lib:${CMAKE_BINARY_DIR}/lib:${CSPARSE_ROOT}/lib:${SCIMAFIC_ROOT}/lib:${HIPOLY_ROOT}/lib")
            ENDIF(PMFLIB_USE_INST_RPATH)
        ENDIF(PMFLIB_HAVE_RPATH)

    ELSE(PMFLIB_NETWORK)

        SET(PMFLIB_CORE_LIBS
"        -L${CMAKE_BINARY_DIR}/lib \\
        -l${PRMFILE_FLIB_NAME} -L${PRMFILE_ROOT}/lib \\
        -l${LAPACK_LIB_NAME} -L${LAPACK_ROOT}/lib \\
        -l${BLAS_LIB_NAME} -L${BLAS_ROOT}/lib" )

        IF(PMFLIB_HAVE_RPATH)
            IF(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_INSTALL_RPATH}")
            ELSE(PMFLIB_USE_INST_RPATH)
                SET(PMFLIB_RPATHS "${CMAKE_BINARY_DIR}/lib::${PRMFILE_ROOT}/lib:${LAPACK_ROOT}/lib:${BLAS_ROOT}/lib")
            ENDIF(PMFLIB_USE_INST_RPATH)
        ENDIF(PMFLIB_HAVE_RPATH)

    ENDIF(PMFLIB_NETWORK)
    # Network ========================================

    # XBPLib ========================================
    IF(PMFLIB_XBPLIB_AVAILABLE)
        SET(PMFLIB_LIBATOMS_XDYNBP_LIBS
"        -lfpmfdrv_libatoms_xdynbp -lfpmf_egap \\
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

            SET(PMFLIB_LIBATOMS_XDYNBP_LIBS
"${PMFLIB_LIBATOMS_XDYNBP_LIBS} \\
        -Wl,-rpath=${PMFLIB_XDYNBP_RPATHS}" )
        ENDIF(PMFLIB_HAVE_RPATH)

    ENDIF(PMFLIB_XBPLIB_AVAILABLE)
    # XBPLib ========================================

    SET(PMFLIB_LIBATOMS_LIBS
"        -lfpmfdrv_libatoms -lfpmf \\
${PMFLIB_CORE_LIBS}" )

    IF(PMFLIB_HAVE_RPATH)
        SET(PMFLIB_LIBATOMS_LIBS
"${PMFLIB_LIBATOMS_LIBS} \\
        -Wl,-rpath=${PMFLIB_RPATHS}" )
    ENDIF(PMFLIB_HAVE_RPATH)

    file(WRITE ${CMAKE_BINARY_DIR}/include/libatoms-dev.inc "PMFLIB_LIBATOMS_EXT=\\\n${PMFLIB_LIBATOMS_LIBS}\n\n")
    file(APPEND ${CMAKE_BINARY_DIR}/include/libatoms-dev.inc "PMFLIB_LIBATOMS_XDYNBP_EXT=\\\n${PMFLIB_LIBATOMS_XDYNBP_LIBS}\n\n")

ENDIF(PMFLIB_LIBATOMS_DRV OR PMFLIB_LIBATOMS_XDYNBP_DRV)
