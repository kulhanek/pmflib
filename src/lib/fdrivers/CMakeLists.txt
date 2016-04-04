# ==============================================================================
# PMFLib CMake File
# ==============================================================================

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/drivers)

# include subdirectories -------------------------------------------------------

# sander drv =================
IF(PMFLIB_SANDER_DRV)
    ADD_SUBDIRECTORY(sander)
ENDIF(PMFLIB_SANDER_DRV)

# pmemd drv =================
IF(PMFLIB_PMEMD_DRV)
    ADD_SUBDIRECTORY(pmemd)
ENDIF(PMFLIB_PMEMD_DRV)

# cpmd drv ===================
IF(PMFLIB_CPMD_DRV)
    ADD_SUBDIRECTORY(cpmd)
ENDIF(PMFLIB_CPMD_DRV)

# cpmd drv ===================
IF(PMFLIB_CP2K_DRV)
   # ADD_SUBDIRECTORY(cp2k)
ENDIF(PMFLIB_CP2K_DRV)

# cpmd drv ===================
IF(PMFLIB_XDYNBP_DRV)
   ADD_SUBDIRECTORY(xdynbp)
ENDIF(PMFLIB_XDYNBP_DRV)

# cats drv ===================
IF(PMFLIB_CATS_DRV)
   ADD_SUBDIRECTORY(cats)
ENDIF(PMFLIB_CATS_DRV)

# # libatoms drv ===================
# IF(PMFLIB_LIBATOMS_DRV)
#     ADD_SUBDIRECTORY(libatoms)
# ENDIF(PMFLIB_LIBATOMS_DRV)

