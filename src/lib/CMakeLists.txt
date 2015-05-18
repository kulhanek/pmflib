# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# these is requird for cpmf library --------------------------------------------
INCLUDE_DIRECTORIES(cpmf SYSTEM)

# include subdirectories -------------------------------------------------------
# C++ part of PMFLib
ADD_SUBDIRECTORY(cpmf)

# fortran part of PMFLib
IF(PMFLIB_COMPILE_FPMF)
    ADD_SUBDIRECTORY(fpmf)
    ADD_SUBDIRECTORY(fdrivers)
    ADD_SUBDIRECTORY(fpmfdyn)
ENDIF(PMFLIB_COMPILE_FPMF)


