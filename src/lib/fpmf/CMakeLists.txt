# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# core objects -----------------------------------------------------------------

SET(PMF_COM_SRC
        common/pmf_sizes.f90
        common/pmf_constants.f90
        common/pmf_utils.F90
        common/pmf_dat.f90
        common/pmf_init.F90
        common/pmf_control.F90
        common/pmf_control_utils.f90
        common/pmf_core.f90
        common/pmf_core_lf.f90
        common/pmf_pbc.f90
        common/pmf_cvs.f90
        common/pmf_paths.f90
        common/pmf_mask.f90
        common/pmf_timers.f90
        common/pmf_finalize.f90
        common/pmf_unit.f90
        common/pmf_ver.f90
        common/pmf_spline.f90
        common/pmf_alloc_cv.F90
        common/pmf_exit.f90
        common/pmf_accu.f90
        common/pmf_kinene.f90
        )

# coordinates ------------------------------------------------------------------

SET(PMF_CVS_SRC
        common/cvs/cv_common.f90
        common/cvs/cv_math.f90

    # distances
        common/cvs/distance/cv_pos.f90
        common/cvs/distance/cv_rad.f90
        common/cvs/distance/cv_opos.f90
        common/cvs/distance/cv_opos2.f90
        common/cvs/distance/cv_orad.f90
        common/cvs/distance/cv_orad2.f90
        common/cvs/distance/cv_dis.f90
        common/cvs/distance/cv_odis.f90
        common/cvs/distance/cv_mdis.f90
        common/cvs/distance/cv_mdisg.f90
        common/cvs/distance/cv_ds.f90
        common/cvs/distance/cv_dd.f90
        common/cvs/distance/cv_ppdis.f90

    # angle based
        common/cvs/angle/cv_ang.f90
        common/cvs/angle/cv_ang2.f90
        common/cvs/angle/cv_pang.f90
        common/cvs/angle/cv_pvang.f90
        common/cvs/angle/cv_axang.f90
        common/cvs/angle/cv_axang2.f90

        common/cvs/angle/cv_cang.f90
        common/cvs/angle/cv_cang2.f90
        common/cvs/angle/cv_cpang.f90
        common/cvs/angle/cv_cpvang.f90
        common/cvs/angle/cv_caxang.f90
        common/cvs/angle/cv_caxang2.f90

        common/cvs/angle/cv_dih.f90
        common/cvs/angle/cv_dih2.f90

    # shape
        common/cvs/shape/cv_rmsdt.f90
        common/cvs/shape/cv_rmsds.f90
        common/cvs/shape/cv_rgyr.f90
        common/cvs/shape/cv_pmgtd.f90
        common/cvs/shape/cv_pmogt.f90
        common/cvs/shape/cv_asph.f90
        common/cvs/shape/cv_acyl.f90
        common/cvs/shape/cv_sanis.f90
        common/cvs/shape/cv_plane.f90
        common/cvs/shape/cv_evec.f90

    # energy based
        common/cvs/energy/cv_epot.f90

    # coordination numbers
        common/cvs/cn/cv_cnff.f90
        common/cvs/cn/cv_cngff.f90
        common/cvs/cn/cv_cngfa.f90
        common/cvs/cn/cv_cnrf.f90
        common/cvs/cn/cv_cngrf.f90
        common/cvs/cn/cv_cnsw.f90
        common/cvs/cn/cv_cngsw.f90

    # nucleic acids
        common/cvs/na/cv_napbo.f90
        common/cvs/na/cv_nasbo.f90
        common/cvs/na/cv_nalbpp.f90
        common/cvs/na/cv_nasbpp.f90
        common/cvs/na/cv_naybpp.f90
        common/cvs/na/cv_nasbppold.f90
        common/cvs/na/cv_naspsxy.f90
        common/cvs/na/cv_nalstp.f90
        common/cvs/na/cv_nasstp.f90
        common/cvs/na/cv_nasstpold.f90
        common/cvs/na/cv_nabend.f90

    # ring puckering
        common/cvs/pucker/cv_puck5q.f90
        common/cvs/pucker/cv_puck5p.f90
        common/cvs/pucker/cv_puck6q.f90
        common/cvs/pucker/cv_puck6t.f90
        common/cvs/pucker/cv_puck6p.f90

    # worm
        common/cvs/worm/cv_wormpos.f90
        common/cvs/worm/cv_wormang.f90

    # path
        common/cvs/path/cv_paths.f90
        common/cvs/path/cv_pathz.f90

    # algebra
        common/cvs/algebra/cv_add.f90
        common/cvs/algebra/cv_sub.f90
        common/cvs/algebra/cv_mul.f90
        common/cvs/algebra/cv_div.f90
        common/cvs/algebra/cv_fswitch.f90
        common/cvs/algebra/cv_rswitch.f90
       )

# ABF --------------------------------------------------------------------------

SET(PMF_ABF_CORE_SRC
        methods/abf/abf_dat.f90
        methods/abf/abf_init.f90
        methods/abf/abf_cvs.f90
        methods/abf/abf_control.F90
        methods/abf/abf_output.f90
        methods/abf/abf_restart.f90
        methods/abf/abf_accu.f90
        methods/abf/abf_core.f90
        methods/abf/abf_core_lf.f90
        methods/abf/abf_client.F90
        methods/abf/abf_trajectory.f90
        methods/abf/abf_finalize.f90
        )

# ABP -------------------------------------------------------------------------

SET(PMF_ABP_CORE_SRC
        methods/abp/abp_dat.f90
        methods/abp/abp_init.f90
        methods/abp/abp_cvs.f90
        methods/abp/abp_control.F90
        methods/abp/abp_output.f90
        methods/abp/abp_restart.f90
        methods/abp/abp_accu.f90
        methods/abp/abp_core.f90
        methods/abp/abp_client.F90
        methods/abp/abp_trajectory.f90
        methods/abp/abp_finalize.f90
        )

# constrained dynamics ----------------------------------------------------------

SET(PMF_CST_CORE_SRC
        methods/cst/cst_dat.f90
        methods/cst/cst_init.F90
        methods/cst/cst_constraints.f90
        methods/cst/cst_control.F90
        methods/cst/cst_output.f90
        methods/cst/cst_restart.f90
        methods/cst/cst_core.f90
        methods/cst/cst_lambdas.f90
        methods/cst/cst_velocities.f90
        methods/cst/cst_shake.f90
        methods/cst/cst_finalize.f90
        methods/cst/cst_trajectory.f90
        methods/cst/cst_accu.f90
        )

# metadynamics -----------------------------------------------------------------

SET(PMF_MTD_CORE_SRC
        methods/mtd/mtd_dat.f90
        methods/mtd/mtd_init.f90
        methods/mtd/mtd_cvs.f90
        methods/mtd/mtd_control.F90
        methods/mtd/mtd_client.F90
        methods/mtd/mtd_output.f90
        methods/mtd/mtd_restart.f90
        methods/mtd/mtd_core.f90
        methods/mtd/mtd_finalize.f90
        methods/mtd/mtd_accu.f90
        methods/mtd/mtd_trajectory.f90
        )

# restrained dynamics -----------------------------------------------------------

SET(PMF_RST_CORE_SRC
        methods/rst/rst_dat.f90
        methods/rst/rst_init.f90
        methods/rst/rst_restraints.f90
        methods/rst/rst_control.F90
        methods/rst/rst_output.f90
        methods/rst/rst_accu.f90
        methods/rst/rst_restart.f90
        methods/rst/rst_core.f90
        methods/rst/rst_finalize.f90
        )

# string method ----------------------------------------------------------------

SET(PMF_STM_CORE_SRC
        methods/stm/stm_dat.f90
        methods/stm/stm_init.f90
        methods/stm/stm_cvs.f90
        methods/stm/stm_control.F90
        methods/stm/stm_output.f90
        methods/stm/stm_core.f90
        methods/stm/stm_client.F90
        methods/stm/stm_finalize.f90
        )

# monitoring -------------------------------------------------------------------

SET(PMF_MON_CORE_SRC
        methods/mon/mon_dat.f90
        methods/mon/mon_init.f90
        methods/mon/mon_control.f90
        methods/mon/mon_output.f90
        methods/mon/mon_core.f90
        methods/mon/mon_finalize.f90
        )

# path driving -----------------------------------------------------------------

SET(PMF_PDRV_CORE_SRC
        methods/pdrv/pdrv_dat.f90
        methods/pdrv/pdrv_paths.f90
        methods/pdrv/pdrv_init.f90
        methods/pdrv/pdrv_control.f90
        methods/pdrv/pdrv_output.f90
        methods/pdrv/pdrv_core.f90
        methods/pdrv/pdrv_finalize.f90
        )

# full list of objects ---------------------------------------------------------

SET(FPMF_SRC
        ${PMF_COM_SRC}
        ${PMF_CVS_SRC}
        ${PMF_ABF_CORE_SRC}
        ${PMF_ABP_CORE_SRC}
        ${PMF_CST_CORE_SRC}
        ${PMF_MTD_CORE_SRC}
        ${PMF_RST_CORE_SRC}
        ${PMF_STM_CORE_SRC}
        ${PMF_MON_CORE_SRC}
        ${PMF_PDRV_CORE_SRC}
        )

IF(PMFLIB_NETWORK)
    ADD_DEFINITIONS(-DPMFLIB_NETWORK)
ENDIF(PMFLIB_NETWORK)

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(fpmf_static STATIC ${FPMF_SRC})
    SET_TARGET_PROPERTIES(fpmf_static PROPERTIES
                            OUTPUT_NAME fpmf
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(fpmf_shared SHARED ${FPMF_SRC})
    ADD_DEPENDENCIES(fpmf_shared cpmf_shared)

    SET_TARGET_PROPERTIES(fpmf_shared PROPERTIES
                            OUTPUT_NAME fpmf
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION}
                            )
    TARGET_LINK_LIBRARIES(fpmf_shared ${CPMF_LIBRARIES} ${FPMF_LIBRARIES_DEPS})
    INSTALL(TARGETS
                fpmf_shared
            DESTINATION
                lib)
ENDIF(LIBS_SHARED)



