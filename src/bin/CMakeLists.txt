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

# PMF Accumulator
ADD_SUBDIRECTORY(accu-info)
ADD_SUBDIRECTORY(accu-combine)
ADD_SUBDIRECTORY(accu-trajectory)

# ABF utilities
ADD_SUBDIRECTORY(abf-integrate)
ADD_SUBDIRECTORY(abf-enthalpy)

# MTD utilities
ADD_SUBDIRECTORY(mtd-energy)

# ABP utilities
ADD_SUBDIRECTORY(abp-energy)

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
