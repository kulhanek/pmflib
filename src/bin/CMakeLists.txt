# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# subprojects ------------------------------------------------------------------

# fortran utilities
IF(PMFLIB_COMPILE_FPMF)
    ADD_SUBDIRECTORY(pmf-dyn)
    ADD_SUBDIRECTORY(pmf-path-explore)
    ADD_SUBDIRECTORY(pmf-path-reconstruct)
    ADD_SUBDIRECTORY(pmf-path-derivatives)
    ADD_SUBDIRECTORY(pmf-test-coords)
ENDIF(PMFLIB_COMPILE_FPMF)

# classical utilities ----------------------------------------------------------

# PMFAccumulator
ADD_SUBDIRECTORY(accu-info)
ADD_SUBDIRECTORY(accu-combine)
ADD_SUBDIRECTORY(accu-trajectory)

# ABF
#ADD_SUBDIRECTORY(abf-energy)
#ADD_SUBDIRECTORY(abf-energy-optgprhyprms)
#ADD_SUBDIRECTORY(abf-enthalpy)
#ADD_SUBDIRECTORY(abf-enthalpy-optgprhyprms)

#
#ADD_SUBDIRECTORY(abf-samples)
#ADD_SUBDIRECTORY(abf-trajectory)
#ADD_SUBDIRECTORY(abf-mask)
#
## ABP utilities
#ADD_SUBDIRECTORY(abp-energy)
#
## general utilities
#ADD_SUBDIRECTORY(pmf-integrate)
#
## metadynamics
#ADD_SUBDIRECTORY(mtd-energy)

## string method
#ADD_SUBDIRECTORY(stm-trajectory)
#
## network utilities ------------------------------------------------------------
IF(PMFLIB_NETWORK)
    # PMFLib MWA utilities
    ADD_SUBDIRECTORY(mwa-server)
    ADD_SUBDIRECTORY(mwa-admin)

    # STM MWA utilities
    ADD_SUBDIRECTORY(stm-server)
    ADD_SUBDIRECTORY(stm-admin)
ENDIF(PMFLIB_NETWORK)
