!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

module cv_egap

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeEGAP

    integer                 :: nstates
    real(PMFDP),allocatable :: alphas(:)

    contains
        procedure :: load_cv        => load_egap
        procedure :: calculate_cv   => calculate_egap
end type CVTypeEGAP

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_egap
!===============================================================================

subroutine load_egap(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common
    use xbp_fepfile_dat
    use xbp_topfile_dat

    implicit none
    class(CVTypeEGAP)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: i, stat, alloc_failed
    character(len=80)                   :: value
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'EGAP'
    cv_item%unit          = EnergyUnit
    cv_item%gradforanycrd = .false.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1

    ! allocate grps array
    allocate(cv_item%grps(cv_item%ngrps),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory!')
    end if

    cv_item%natoms  = fnatoms
    cv_item%grps(1) = cv_item%natoms
    cv_item%nstates = nstates

    write(PMF_OUT,10) cv_item%natoms
    write(PMF_OUT,20) cv_item%nstates

    ! allocate arrays
    allocate(cv_item%rindexes(cv_item%natoms), &
             cv_item%lindexes(cv_item%natoms), &
             cv_item%alphas(cv_item%nstates), stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory!')
    end if

    do i=1,cv_item%natoms
        cv_item%rindexes(i) = i
    end do

     ! check if this is EVB or mapping calculation -----------------------------
     select case(SurfaceType)
        case(SURFACE_MAPPING,SURFACE_EVB,SURFACE_RST_EGAP)
            ! this is allowed
        case default
            call pmf_utils_exit(PMF_OUT,1,'EGAP can be used only for simulation on mapping/EVB/EGAP energy surface!')
     end select

     ! load mixing values ------------------------------------------------------
     if( .not. prmfile_init_field_by_key(prm_fin,'mixing') ) then
        call pmf_utils_exit(PMF_OUT,1,'mixing is not specified!')
     end if

     do i=1,nstates
        if( .not. prmfile_get_field_by_key(prm_fin,value) ) then
            write(PMF_OUT,'(/,A,I1,A)') 'Missing mixing record for state ',i,'!'
            call pmf_utils_exit(PMF_OUT,1,'Incorrect mixing record!')
        end if
        read(value,*,iostat=stat) cv_item%alphas(i)
        if( stat .ne. 0 ) then
            write(PMF_OUT,'(/,A,I1,A)') 'Mixing record for state ',i,' is not real number!'
            call pmf_utils_exit(PMF_OUT,1,'Incorrect mixing record!')
        end if
        write(PMF_OUT,30) i, cv_item%alphas(i)
     end do

 10 format('   ** number of atoms    : ',I9)
 20 format('   ** number of states   : ',I9)
 30 format('   ** mixing value (',I2,')   :',F10.3)

end subroutine load_egap

!===============================================================================
! Subroutine:  calculate_egap
!===============================================================================

subroutine calculate_egap(cv_item,x,ctx)

    use pmf_dat
    use xbp_fepfile_dat
    use xbp_nonbonded
    use xbp_qbonded
    use xbp_energy_utils
    use xbp_energy_dat

    implicit none
    class(CVTypeEGAP)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ii,istate,ji
    type(Q_ENERGIES)    :: EQ(nstates)
    ! --------------------------------------------------------------------------

    ! original coordinates
    qx(:,:) = md_x(:,:)
    qd(:,:,:) = 0.0

    ! make copy of x
    do i=1,cv_item%natoms
        qx(:,cv_item%rindexes(i)) = x(:,cv_item%lindexes(i))
    end do

    ! clear energies
    call clear_qenergy(EQ)

    ! calculate q-bonded part
    call potene_qbonded(qx,qd,EQ)

    ! calculate q-nonbonded part
    call potene_nb_qonly(qx,qd,EQ)

    ! sum all energies together
    do istate = 1, nstates
    EQ(istate)%qx%el  = EQ(istate)%qq%el +EQ(istate)%qp%el +EQ(istate)%qw%el +EQ(istate)%lr%el
    EQ(istate)%qx%lr  = EQ(istate)%qq%lr +EQ(istate)%qp%lr +EQ(istate)%qw%lr +EQ(istate)%lr%lr
    EQ(istate)%qx%vdw = EQ(istate)%qq%vdw+EQ(istate)%qp%vdw+EQ(istate)%qw%vdw+EQ(istate)%lr%vdw

    EQ(istate)%total =  EQ(istate)%q%bond       + EQ(istate)%q%angle    + &
                        EQ(istate)%q%torsion    + EQ(istate)%q%improper + &
                        EQ(istate)%qx%el        + EQ(istate)%qx%lr      + EQ(istate)%qx%vdw + &
                        EQ(istate)%restraint
    end do

    ! finalize energy and gradients
    ctx%CVsValues(cv_item%idx) = 0.0
    do istate = 1, nstates
        ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + cv_item%alphas(istate)*(EQ(istate)%total + Egalphas(istate))
    end do

    ! --------------------------------------------------------------------------

    ! calc derivatives
    do i=1,cv_item%natoms
        ii = cv_item%rindexes(i)
        ji = cv_item%lindexes(i)
        do istate = 1, nstates
            ctx%CVsDrvs(:,ji,cv_item%idx) = ctx%CVsDrvs(:,ji,cv_item%idx) + cv_item%alphas(istate)*qd(:,ii,istate)
        end do
    end do

    return

end subroutine calculate_egap

!===============================================================================

end module cv_egap


