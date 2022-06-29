!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cst_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_core_main_lf
!===============================================================================

subroutine cst_core_main_lf

    use cst_constraints
    use cst_lambdas
    use cst_output
    use cst_restart
    use cst_trajectory

    implicit none
    ! --------------------------------------------------------------------------

    call cst_constraints_increment
    call cst_lambdas_calculate
    call cst_core_analyze
    call cst_output_write
    call cst_restart_update
    call cst_trajectory_write_snapshot

end subroutine cst_core_main_lf

!===============================================================================
! Subroutine:  cst_core_main_vv_shake
!===============================================================================

subroutine cst_core_main_vv_shake

    use cst_lambdas
    use cst_output
    use cst_velocities

    implicit none
    ! --------------------------------------------------------------------------

    ! FIXME
    call cst_lambdas_calculate
    call cst_velocities_correct_a
    call cst_core_analyze
    call cst_output_write

end subroutine cst_core_main_vv_shake

!===============================================================================
! Subroutine:  cst_core_main_vv_rattle
!===============================================================================

subroutine cst_core_main_vv_rattle()

    use cst_constraints
    use cst_velocities

    implicit none
    ! --------------------------------------------------------------------------

    call cst_constraints_increment
    call cst_velocities_correct_b

end subroutine cst_core_main_vv_rattle

!===============================================================================
! Subroutine:  cst_core_analyze
!===============================================================================

subroutine cst_core_analyze

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,ci,j,cj,k,info
    real(PMFDP)            :: fzv,isrz,lam,mu,dval1,dval2,invn
    real(PMFDP)            :: epot, ekin, etot, erst
    real(PMFDP)            :: detot1, detot2
    real(PMFDP)            :: depot1, depot2
    real(PMFDP)            :: dekin1, dekin2
    real(PMFDP)            :: derst1, derst2
    ! --------------------------------------------------------------------------

    ! reset accumulators ---------------------------------------------------
    if ( faccurst .eq. 0 ) then
        faccumulation   = 0
        misrz           = 0
        m2isrz          = 0
        faccurst        = -1

        mlambda(:)      = 0.0d0
        m2lambda(:)     = 0.0d0

        if( has_lambdav ) then
            mlambdav(:) = 0.0d0
            m2lambdav(:) = 0.0d0
        end if

        lambda0(:)  = 0.0d0
        lambda1(:)  = 0.0d0

        if( fenthalpy .or. fentropy ) then
            c11hh(:)    = 0.0d0
            c11hp(:)    = 0.0d0
            c11hk(:)    = 0.0d0
            c11hr(:)    = 0.0d0
        end if

        metot       = 0.0d0
        m2etot      = 0.0d0
        mepot       = 0.0d0
        m2epot      = 0.0d0
        mekin       = 0.0d0
        m2ekin      = 0.0d0
        merst       = 0.0d0
        m2erst      = 0.0d0

        epothist0   = 0.0d0
        epothist1   = 0.0d0

        ersthist0   = 0.0d0
        ersthist1   = 0.0d0

        isrz0       = 0.0d0
        isrz1       = 0.0d0

        CONList(:)%sdevtot = 0.0d0

        write(CST_OUT,'(A)') '#-------------------------------------------------------------------------------'
        write(CST_OUT,'(A)') '# INFO: ALL ACCUMULATORS WERE RESETED                                           '
        write(CST_OUT,'(A)') '#       PRODUCTION STAGE OF ACCUMULATION IS STARTED                             '
        write(CST_OUT,'(A)') '#-------------------------------------------------------------------------------'
    end if

    ! accumulate results -----------------------------------------------------
    if( faccurst .gt. 0 ) then
        faccurst = faccurst - 1
    end if

! calculate final constraint deviations
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        CONList(i)%deviation = CONList(i)%cv%get_deviation(CVContextP%CVsValues(ci),CONList(i)%value)   ! t+dt
        CONList(i)%sdevtot = CONList(i)%sdevtot + CONList(i)%deviation**2                               ! t+dt
    end do

    ! calculate Z matrix at Crd (in t)
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfCONs
            cj = CONList(j)%cvindx
            fzv = 0.0
            do k=1,NumOfLAtoms
                fzv = fzv + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
            fz(i,j) = fzv
        end do
    end do

    ! calculate Z determinant ------------------------------------
    if( NumOfCONs .gt. 1 ) then
        ! LU decomposition
        call dgetrf(NumOfCONs,NumOfCONs,fz,NumOfCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_main!')
        end if
        fzdet = 1.0d0
        ! and finally determinant
        do i=1,NumOfCONs
            if( indx(i) .ne. i ) then
                fzdet = - fzdet * fz(i,i)
            else
                fzdet = fzdet * fz(i,i)
            end if
        end do
    else
        fzdet = fz(1,1)
    end if

    ! calculate metric tensor correction --------------------------------------------------------
    isrz    = 1.0d0/sqrt(fzdet)         ! t

! record history of isrz
    isrz0 = isrz1                       ! t-dt
    isrz1 = isrz                        ! t

! record history of lambda
    lambda0(:) = lambda1(:)             ! t-dt
    lambda1(:) = lambda(:)              ! t

! record history of Epot
    epothist0  = epothist1              ! t-dt
    epothist1  = PotEne + fepotoffset   ! t

    ersthist0 = ersthist1               ! t-dt
    ersthist1  = PMFEne                 ! t

    ! KinEne                            ! t-dt

     ! do we have enough samples?
    if( fstep .le. 2 ) return

    faccumulation = faccumulation + 1

    ! this will occur for velocity verlet algorithm
    if( faccumulation .le. 0 ) return

    invn = 1.0d0/real(faccumulation,PMFDP)

! values
    epot = epothist0                    ! t-dt
    ekin = KinEne%KinEneVV + fekinoffset         ! t-dt
    erst = ersthist0                    ! t-dt
    etot = epot + ekin + erst           ! t-dt
    isrz = isrz0                        ! t-dt
    ! lambda0(i)                        ! t-dt

! isrz
    dval1   = isrz - misrz
    misrz   = misrz  + dval1 * invn
    dval2   = isrz - misrz
    m2isrz  = m2isrz + dval1*dval2

! total energy
    detot1 = etot - metot
    metot  = metot  + detot1 * invn
    detot2 = etot - metot
    m2etot = m2etot + detot1 * detot2

! potential energy
    depot1 = epot - mepot
    mepot  = mepot  + depot1 * invn
    depot2 = epot - mepot
    m2epot = m2epot + depot1 * depot2

! potential energy
    dekin1 = ekin - mekin
    mekin  = mekin  + dekin1 * invn
    dekin2 = ekin - mekin
    m2ekin = m2ekin + dekin1 * dekin2

! restraint energy
    derst1 = erst - merst
    merst  = merst  + derst1 * invn
    derst2 = erst - merst
    m2erst = m2erst + derst1 * derst2

! lambda and entropy
    do i=1,NumOfCONs

        if( has_lambdav ) then
            ! FIXME - lambdav is not in correct time?
            mu              = lambdav(i)
            dval1           = mu - mlambdav(i)
            mlambdav(i)     = mlambdav(i)  + dval1 * invn
            dval2           = mu - mlambdav(i)
            m2lambdav(i)    = m2lambdav(i) + dval1 * dval2
        end if


        ! lambda ----------------------------------
        lam             = lambda0(i)        ! t-dt
        dval1           = lam - mlambda(i)
        mlambda(i)      = mlambda(i)  + dval1 * invn
        dval2           = lam - mlambda(i)
        m2lambda(i)     = m2lambda(i) + dval1 * dval2

        if( fentropy ) then
            c11hh(i)   = c11hh(i) + dval1 * detot2
            c11hp(i)   = c11hp(i) + dval1 * depot2
            c11hk(i)   = c11hk(i) + dval1 * dekin2
            c11hr(i)   = c11hr(i) + dval1 * derst2
        end if
    end do

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>> ', (lambda0(i), i=1,NumOfCONs), etot, epot, ekin, erst
    end if

end subroutine cst_core_analyze

!===============================================================================

end module cst_core

