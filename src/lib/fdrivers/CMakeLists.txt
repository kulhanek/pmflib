# ==============================================================================
# PMFLib CMake File
# ==============================================================================

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/drivers)

# include sub-directories ------------------------------------------------------

# sander drv =================
IF(PMFLIB_SANDER_DRV)
   # ADD_SUBDIRECTORY(sander)
ENDIF(PMFLIB_SANDER_DRV)

# sander drv =================
IF(PMFLIB_SANDER_DYNAMIC_DRV)
    ADD_SUBDIRECTORY(sander-dynamic)
ENDIF(PMFLIB_SANDER_DYNAMIC_DRV)

# pmemd drv =================
IF(PMFLIB_PMEMD_DYNAMIC_DRV)
    ADD_SUBDIRECTORY(pmemd-dynamic)
ENDIF(PMFLIB_PMEMD_DYNAMIC_DRV)

# cats drv ===================
IF(PMFLIB_CATS_DRV)
   ADD_SUBDIRECTORY(cats)
ENDIF(PMFLIB_CATS_DRV)


