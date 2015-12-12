# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PMF_DYN_SRC
        pmf-dyn.f90
        pmfdyn_mypotene.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(pmf-dyn ${PMF_DYN_SRC})
ADD_DEPENDENCIES(pmf-dyn fpmfdyn_shared fpmf_shared)

IF(PMFLIB_NETWORK)
    TARGET_LINK_LIBRARIES(pmf-dyn
                fpmfdyn fpmf
                ${PRMFILE_FLIB_NAME} ${LAPACK_LIB_NAME} ${BLAS_LIB_NAME}
                ${SCIMAFIC_FLIB_NAME} ${RANLUX_LIB_NAME} ${CPMF_LIBRARIES})
ELSE(PMFLIB_NETWORK)
    TARGET_LINK_LIBRARIES(pmf-dyn
                fpmfdyn fpmf
                ${PRMFILE_FLIB_NAME} ${LAPACK_LIB_NAME}
                ${BLAS_LIB_NAME} ${SCIMAFIC_FLIB_NAME} ${RANLUX_LIB_NAME})
ENDIF(PMFLIB_NETWORK)

INSTALL(TARGETS
            pmf-dyn
        DESTINATION
            bin
        )
