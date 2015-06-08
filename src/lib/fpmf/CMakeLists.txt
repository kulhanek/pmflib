# ==============================================================================
# PMFLib CMake File
# ==============================================================================

# core objects -----------------------------------------------------------------

SET(PMF_COM_SRC
        common/pmf_sizes.f90
        common/pmf_constants.f90
        common/pmf_utils.f90
        common/pmf_dat.f90
        common/pmf_init.F90
        common/pmf_control.F90
        common/pmf_core.f90
        common/pmf_core_lf.f90
        common/pmf_core_vv.f90
        common/pmf_pbc.f90
        common/pmf_cvs.f90
        common/pmf_paths.f90
        common/pmf_mask.f90
        common/pmf_timers.f90
        common/pmf_finalize.f90
        common/pmf_unit.f90
        common/pmf_utils.f90
        common/pmf_ver.f90
        common/pmf_spline.f90
        common/pmf_alloc_cv.F90
        common/pmf_exit.f90
        )

# coordinates ------------------------------------------------------------------

SET(PMF_CVS_SRC
        common/cvs/cv_common.f90

    # distances
        common/cvs/distance/cv_pos.f90
        common/cvs/distance/cv_rad.f90
        common/cvs/distance/cv_opos.f90
        common/cvs/distance/cv_orad.f90
        common/cvs/distance/cv_dis.f90
        common/cvs/distance/cv_odis.f90
        common/cvs/distance/cv_mdis.f90
        common/cvs/distance/cv_mdisg.f90
        common/cvs/distance/cv_ds.f90
        common/cvs/distance/cv_dd.f90
 
    # angle based
        common/cvs/angle/cv_ang.f90
        common/cvs/angle/cv_pang.f90
        common/cvs/angle/cv_pvang.f90
        common/cvs/angle/cv_axang.f90
        common/cvs/angle/cv_dih.f90

    # shape
        common/cvs/shape/cv_rmsdt.f90
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

    # ring puckering
        common/cvs/pucker/cv_puck5q.f90
        common/cvs/pucker/cv_puck5p.f90
        common/cvs/pucker/cv_puck6q.f90
        common/cvs/pucker/cv_puck6t.f90
        common/cvs/pucker/cv_puck6p.f90

    # algebra
        common/cvs/algebra/cv_add.f90
        common/cvs/algebra/cv_sub.f90
        common/cvs/algebra/cv_mul.f90
        common/cvs/algebra/cv_div.f90
       )

# gaussian process -------------------------------------------------------------

SET(PMF_GP_SRC
        common/gp/positions.f90
        common/gp/clustering.f90
        common/gp/gp_dat.f90
        common/gp/gp_basic.f90
   )

# XBPLib binding ---------------------------------------------------------------

IF(PMFLIB_HAVE_XBPLIB)
    SET(PMF_CVS_SRC ${PMF_CVS_SRC}
        # EGAP based
            common/cvs/energy/cv_egap.f90
        )
    ADD_DEFINITIONS(-DPMFLIB_HAVE_XBPLIB)
ENDIF(PMFLIB_HAVE_XBPLIB)

# ABF --------------------------------------------------------------------------

SET(PMF_ABF_CORE_SRC
        methods/abf/abf_dat.f90
        methods/abf/abf_init.f90
        methods/abf/abf_cvs.f90
        methods/abf/abf_control.F90
        methods/abf/abf_output.f90
        methods/abf/abf_restart.f90
        methods/abf/abf_accumulator.f90
        methods/abf/abf_core.f90
        methods/abf/abf_client.F90
        methods/abf/abf_trajectory.f90
        methods/abf/abf_finalize.f90
        methods/abf/abf_gprocess.f90
        )

# constraint dynamics ----------------------------------------------------------

SET(PMF_CON_CORE_SRC
        methods/con/con_dat.f90
        methods/con/con_init.F90
        methods/con/con_constraints.f90
        methods/con/con_control.F90
        methods/con/con_output.f90
        methods/con/con_restart.f90
        methods/con/con_core.f90
        methods/con/con_lambdas.f90
        methods/con/con_velocities.f90
        methods/con/con_shake.f90
        methods/con/con_finalize.f90
        )

# metadynamics -----------------------------------------------------------------

SET(PMF_MTD_CORE_SRC
        methods/mtd/mtd_dat.f90
        methods/mtd/mtd_init.f90
        methods/mtd/mtd_cvs.f90
        methods/mtd/mtd_history.f90
        methods/mtd/mtd_control.F90
        methods/mtd/mtd_client.F90
        methods/mtd/mtd_output.f90
        methods/mtd/mtd_restart.f90
        methods/mtd/mtd_core.f90
        methods/mtd/mtd_finalize.f90
        methods/mtd/mtd_grid.f90
        )

# restraint dynamics -----------------------------------------------------------

SET(PMF_RST_CORE_SRC
        methods/rst/rst_dat.f90
        methods/rst/rst_init.f90
        methods/rst/rst_restraints.f90
        methods/rst/rst_control.F90
        methods/rst/rst_output.f90
        methods/rst/rst_accumulator.f90
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

# remd -------------------------------------------------------------------------

SET(PMF_REMD_CORE_SRC
        methods/remd/remd_dat.f90
        methods/remd/remd_init.f90
        methods/remd/remd_control.F90
        methods/remd/remd_client.F90
        methods/remd/remd_output.f90
        methods/remd/remd_core.f90
        methods/remd/remd_finalize.f90
        )

# gap --------------------------------------------------------------------------

SET(PMF_GAP_CORE_SRC
        methods/gap/gap_dat.f90
        methods/gap/gap_init.f90
        methods/gap/gap_cvs.f90
        methods/gap/gap_control.f90
        methods/gap/gap_output.f90
        methods/gap/gap_core.f90
        methods/gap/gap_finalize.f90
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
        ${PMF_GP_SRC}
        ${PMF_ABF_CORE_SRC}
        ${PMF_CON_CORE_SRC}
        ${PMF_MTD_CORE_SRC}
        ${PMF_RST_CORE_SRC}
        ${PMF_STM_CORE_SRC}
        ${PMF_MON_CORE_SRC}
        ${PMF_REMD_CORE_SRC}
        ${PMF_GAP_CORE_SRC}
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
    SET_TARGET_PROPERTIES(fpmf_shared PROPERTIES
                            OUTPUT_NAME fpmf
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PMFLIB_VERSION}
                            )
    TARGET_LINK_LIBRARIES(fpmf_shared ${CPMF_LIBRARIES} ${FPMF_LIBRARIES_DEPS})
    INSTALL(TARGETS
                fpmf_shared
            LIBRARY DESTINATION
                lib)
ENDIF(LIBS_SHARED)



