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

module pmf_core

implicit none
contains

!===============================================================================
! Subroutine:  pmf_core_in_data_xvf
!===============================================================================

subroutine pmf_core_in_data_xvf(x,v,f)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: x(:,:)
    real(PMFDP),intent(in)     :: v(:,:)
    real(PMFDP),intent(in)     :: f(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        Crd(:,i) = x(:,ridx)*LengthConv
        Vel(:,i) = v(:,ridx)*VelocityConv
        Frc(:,i) = f(:,ridx)*ForceConv
    end do

end subroutine pmf_core_in_data_xvf

!===============================================================================
! Subroutine:  pmf_core_in_data_xf
!===============================================================================

subroutine pmf_core_in_data_xf(x,f)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: x(:,:)
    real(PMFDP),intent(in)     :: f(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        Crd(:,i) = x(:,ridx)*LengthConv
        Frc(:,i) = f(:,ridx)*ForceConv
    end do

end subroutine pmf_core_in_data_xf

!===============================================================================
! Subroutine:  pmf_core_in_data_xv
!===============================================================================

subroutine pmf_core_in_data_xv(x,v)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: x(:,:)
    real(PMFDP),intent(in)     :: v(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        Crd(:,i) = x(:,ridx)*LengthConv
        Vel(:,i) = v(:,ridx)*VelocityConv
    end do

end subroutine pmf_core_in_data_xv

!===============================================================================
! Subroutine:  pmf_core_in_data_xp
!===============================================================================

subroutine pmf_core_in_data_xp(xp)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: xp(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        CrdP(:,i) = xp(:,ridx)*LengthConv
    end do

end subroutine pmf_core_in_data_xp

!===============================================================================
! Subroutine:  pmf_core_in_data_vp
!===============================================================================

subroutine pmf_core_in_data_vp(vp)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: vp(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        VelP(:,i) = vp(:,ridx)*VelocityConv
    end do

end subroutine pmf_core_in_data_vp

!===============================================================================
! Subroutine:  pmf_core_in_data_xpvp
!===============================================================================

subroutine pmf_core_in_data_xpvp(xp,vp)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: xp(:,:)
    real(PMFDP),intent(in)     :: vp(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        CrdP(:,i) = xp(:,ridx)*LengthConv
        VelP(:,i) = vp(:,ridx)*VelocityConv
    end do

end subroutine pmf_core_in_data_xpvp

!===============================================================================
! Subroutine:  pmf_core_out_data_f
!===============================================================================

subroutine pmf_core_out_data_f(f)

    use pmf_dat

    implicit none
    real(PMFDP),intent(out)    :: f(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    real(PMFDP)                :: fconv
    ! --------------------------------------------------------------------------

    fconv = 1.0d0 / ForceConv

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        f(:,ridx) = Frc(:,i)*fconv
    end do

end subroutine pmf_core_out_data_f

!===============================================================================
! Subroutine:  pmf_core_out_data_xv
!===============================================================================

subroutine pmf_core_out_data_xv(x,v)

    use pmf_dat

    implicit none
    real(PMFDP),intent(out)    :: x(:,:)
    real(PMFDP),intent(out)    :: v(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    real(PMFDP)                :: xconv,vconv
    ! --------------------------------------------------------------------------

    xconv = 1.0d0 / LengthConv
    vconv = 1.0d0 / VelocityConv

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        x(:,ridx) = Crd(:,i)*xconv
        v(:,ridx) = Vel(:,i)*vconv
    end do

end subroutine pmf_core_out_data_xv

!===============================================================================
! Subroutine:  pmf_core_out_data_xp
!===============================================================================

subroutine pmf_core_out_data_xp(xp)

    use pmf_dat

    implicit none
    real(PMFDP),intent(out)    :: xp(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    real(PMFDP)                :: xconv
    ! --------------------------------------------------------------------------

    xconv = 1.0d0 / LengthConv

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        xp(:,ridx) = CrdP(:,i)*xconv
    end do

end subroutine pmf_core_out_data_xp

!===============================================================================
! Subroutine:  pmf_core_out_data_vp
!===============================================================================

subroutine pmf_core_out_data_vp(vp)

    use pmf_dat

    implicit none
    real(PMFDP),intent(out)    :: vp(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    real(PMFDP)                :: vconv
    ! --------------------------------------------------------------------------

    vconv = 1.0d0 / VelocityConv

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        vp(:,ridx) = VelP(:,i)*vconv
    end do

end subroutine pmf_core_out_data_vp

!===============================================================================
! Subroutine:  pmf_core_out_data_xpvp
!===============================================================================

subroutine pmf_core_out_data_xpvp(xp,vp)

    use pmf_dat

    implicit none
    real(PMFDP),intent(out)    :: xp(:,:)
    real(PMFDP),intent(out)    :: vp(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    real(PMFDP)                :: xconv,vconv
    ! --------------------------------------------------------------------------

    xconv = 1.0d0 / LengthConv
    vconv = 1.0d0 / VelocityConv

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        xp(:,ridx) = CrdP(:,i)*xconv
        vp(:,ridx) = VelP(:,i)*vconv
    end do

end subroutine pmf_core_out_data_xpvp

!===============================================================================
! Subroutine:  pmf_core_in_data_xbar
!===============================================================================

subroutine pmf_core_in_data_xbar(xbar)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: xbar(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        CrdBar(:,i) = xbar(:,ridx)*LengthConv
    end do

end subroutine pmf_core_in_data_xbar

!===============================================================================
! Subroutine:  pmf_core_in_data_vbar
!===============================================================================

subroutine pmf_core_in_data_vbar(vbar)

    use pmf_dat

    implicit none
    real(PMFDP),intent(in)     :: vbar(:,:)
    ! -----------------------------------------------
    integer                    :: i,ridx
    ! --------------------------------------------------------------------------

    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        VelBar(:,i) = vbar(:,ridx)*VelocityConv
    end do

end subroutine pmf_core_in_data_vbar

!===============================================================================
! Function:  pmf_core_get_num_of_constraints
!===============================================================================

integer function pmf_core_get_num_of_constraints()

    use pmf_dat
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    pmf_core_get_num_of_constraints = 0

    if ( .not. cst_enabled) return

    pmf_core_get_num_of_constraints = NumOfCONs

    return

end function pmf_core_get_num_of_constraints


!===============================================================================

end module pmf_core


