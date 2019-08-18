!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2016 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2012      Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011      Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module pmf_alloc_cv

implicit none
contains

!===============================================================================
! Subroutine:  pmf_alloc_cv_allocate
!===============================================================================

subroutine pmf_alloc_cv_allocate(cv_type,cv_item)

    use prmfile
    use pmf_utils

! distance based ---------------------------------
    use cv_pos
    use cv_rad
    use cv_opos
    use cv_orad
    use cv_dis
    use cv_odis
    use cv_mdis
    use cv_mdisg
    use cv_dd
    use cv_ds
    use cv_ppdis
    use cv_wormod

! angle based ------------------------------------
    use cv_ang
    use cv_ang2
    use cv_pang
    use cv_pvang
    use cv_axang
    use cv_axang2

    use cv_cang
    use cv_cang2
    use cv_cpang
    use cv_cpvang
    use cv_caxang
    use cv_caxang2

    use cv_dih
    use cv_dih2

! shape ------------------------------------------
    use cv_rgyr
    use cv_pmogt
    use cv_pmgtd
    use cv_acyl
    use cv_asph
    use cv_sanis
    use cv_rmsdt
    use cv_rmsds
    use cv_plane
    use cv_evec

! energy based CVs -------------------------------
    use cv_epot
#ifdef PMFLIB_HAVE_XBPLIB
    use cv_egap
#endif

! coordination number  ---------------------------
    use cv_cnff
    use cv_cngff
    use cv_cngfa
    use cv_cnrf
    use cv_cngrf
    use cv_cnsw
    use cv_cngsw

! nucleic acids ---------------------------------
    use cv_napbo
    use cv_nasbo
    use cv_nasbpp
    use cv_nabend

! ring puckering ---------------------------------
    use cv_puck5q
    use cv_puck5p
    use cv_puck6q
    use cv_puck6t
    use cv_puck6p

! algebra ----------------------------------------
    use cv_add
    use cv_sub
    use cv_mul
    use cv_div
    use cv_fswitch
    use cv_rswitch

    implicit none
    character(*)            :: cv_type
    class(CVType),pointer   :: cv_item
    ! --------------------------------------------------------------------------

    select case(cv_type)
    ! distance -----------------------------------
        case('POS')
            allocate(CVTypePOS::cv_item)
        case('RAD')
            allocate(CVTypeRAD::cv_item)
        case('OPOS')
            allocate(CVTypeOPOS::cv_item)
        case('ORAD')
            allocate(CVTypeORAD::cv_item)
        case('DIS')
            allocate(CVTypeDIS::cv_item)
        case('ODIS')
            allocate(CVTypeODIS::cv_item)
        case('MDIS')
            allocate(CVTypeMDIS::cv_item)
        case('MDISG')
            allocate(CVTypeMDISG::cv_item)
        case('DD')
            allocate(CVTypeDD::cv_item)
        case('PPDIS')
            allocate(CVTypePPDIS::cv_item)
        case('WORMOD')
            allocate(CVTypeWORMOD::cv_item)

    ! angle --------------------------------------
        case('ANG')
            allocate(CVTypeANG::cv_item)
        case('ANG2')
            allocate(CVTypeANG2::cv_item)
        case('PANG')
            allocate(CVTypePANG::cv_item)
        case('PVANG')
            allocate(CVTypePVANG::cv_item)
        case('AXANG')
            allocate(CVTypeAXANG::cv_item)
        case('AXANG2')
            allocate(CVTypeAXANG2::cv_item)

        case('CANG')
            allocate(CVTypeCANG::cv_item)
        case('CANG2')
            allocate(CVTypeCANG2::cv_item)
        case('CPANG')
            allocate(CVTypeCPANG::cv_item)
        case('CPVANG')
            allocate(CVTypeCPVANG::cv_item)
        case('CAXANG')
            allocate(CVTypeCAXANG::cv_item)
        case('CAXANG2')
            allocate(CVTypeCAXANG2::cv_item)

        case('DIH')
            allocate(CVTypeDIH::cv_item)
        case('DIH2')
            allocate(CVTypeDIH2::cv_item)

    ! shape --------------------------------------
        case('RGYR')
            allocate(CVTypeRGYR::cv_item)
        case('PMOGT')
            allocate(CVTypePMOGT::cv_item)
        case('PMGTD')
            allocate(CVTypePMGTD::cv_item)
        case('ACYL')
            allocate(CVTypeACYL::cv_item)
        case('ASPH')
            allocate(CVTypeASPH::cv_item)
        case('SANIS')
            allocate(CVTypeSANIS::cv_item)

        case('RMSDT')
            allocate(CVTypeRMSDT::cv_item)
        case('RMSDS')
            allocate(CVTypeRMSDS::cv_item)
        case('PLANE')
            allocate(CVTypePLANE::cv_item)
        case('EVEC')
            allocate(CVTypeEVEC::cv_item)

    ! energy -------------------------------------
        case('EPOT')
            allocate(CVTypeEPOT::cv_item)
#ifdef PMFLIB_HAVE_XBPLIB
        case('EGAP')
            allocate(CVTypeEGAP::cv_item)
#endif

    ! coordination numbers -----------------------
        case('CNFF')
            allocate(CVTypeCNFF::cv_item)
        case('CNGFF')
            allocate(CVTypeCNGFF::cv_item)
        case('CNGFA')
            allocate(CVTypeCNGFA::cv_item)
        case('CNRF')
            allocate(CVTypeCNRF::cv_item)
        case('CNGRF')
            allocate(CVTypeCNGRF::cv_item)
        case('CNSW')
            allocate(CVTypeCNSW::cv_item)
        case('CNGSW')
            allocate(CVTypeCNGSW::cv_item)

! nucleic acids ---------------------------------
        case('NAPBO')
            allocate(CVTypeNAPBO::cv_item)
        case('NASBO')
            allocate(CVTypeNASBO::cv_item)
        case('NASBPP')
            allocate(CVTypeNASBPP::cv_item)
        case('NABEND')
            allocate(CVTypeNABEND::cv_item)

    ! pucker -------------------------------------
        case('PUCK5Q')
            allocate(CVTypePUCK5Q::cv_item)
        case('PUCK5P')
            allocate(CVTypePUCK5P::cv_item)
        case('PUCK6Q')
            allocate(CVTypePUCK6Q::cv_item)
        case('PUCK6T')
            allocate(CVTypePUCK6T::cv_item)
        case('PUCK6P')
            allocate(CVTypePUCK6P::cv_item)

    ! algebra ------------------------------------
        case('ADD')
            allocate(CVTypeADD::cv_item)
        case('SUB')
            allocate(CVTypeSUB::cv_item)
        case('MUL')
            allocate(CVTypeMUL::cv_item)
        case('DIV')
            allocate(CVTypeDIV::cv_item)
        case('FSWITCH')
            allocate(CVTypeFSWITCH::cv_item)
        case('RSWITCH')
            allocate(CVTypeRSWITCH::cv_item)

        case default
            call pmf_utils_exit(PMF_OUT,1,&
                                'This coordinate type is not supported!')
    end select

end subroutine pmf_alloc_cv_allocate

!===============================================================================

end module pmf_alloc_cv
