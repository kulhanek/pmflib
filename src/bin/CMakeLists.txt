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

ADD_SUBDIRECTORY(accu-info)

## ABF utilities
#ADD_SUBDIRECTORY(abf-accuinfo)
#ADD_SUBDIRECTORY(abf-combine)
#ADD_SUBDIRECTORY(abf-derivatives)
#
#ADD_SUBDIRECTORY(abf-energy)
#ADD_SUBDIRECTORY(abf-energy-optgprhyprms)
#ADD_SUBDIRECTORY(abf-enthalpy)
#ADD_SUBDIRECTORY(abf-enthalpy-optgprhyprms)
#ADD_SUBDIRECTORY(abf-entropy)
##ADD_SUBDIRECTORY(abf-entropy-optgprhyprms)
#
#ADD_SUBDIRECTORY(abf-samples)
#ADD_SUBDIRECTORY(abf-trajectory)
#ADD_SUBDIRECTORY(abf-mask)
#
## ABP utilities
#ADD_SUBDIRECTORY(abp-samples)
#ADD_SUBDIRECTORY(abp-trajectory)
#ADD_SUBDIRECTORY(abp-energy)
#
## general utilities
#ADD_SUBDIRECTORY(pmf-integrate)
#
## metadynamics
#ADD_SUBDIRECTORY(mtd-energy)
#
## restraint dynamics
#ADD_SUBDIRECTORY(rst-samples)
#
## string method
#ADD_SUBDIRECTORY(stm-trajectory)
#
## network utilities ------------------------------------------------------------
#IF(PMFLIB_NETWORK)
#    # ABF MWA utilities
#    ADD_SUBDIRECTORY(abf-server)
#    ADD_SUBDIRECTORY(abf-admin)
#
#    # ABP MWA utilities
#    ADD_SUBDIRECTORY(abp-server)
#    ADD_SUBDIRECTORY(abp-admin)
#
#    # MTD MWA utilities
#    #ADD_SUBDIRECTORY(mtd-server)
#    ADD_SUBDIRECTORY(mtd-admin)
#
#    # STM MWA utilities
#    ADD_SUBDIRECTORY(stm-server)
#    ADD_SUBDIRECTORY(stm-admin)
#ENDIF(PMFLIB_NETWORK)
