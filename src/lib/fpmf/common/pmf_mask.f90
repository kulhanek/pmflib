!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
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

module pmf_mask

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================
! CPMF/FPMF interface
!===============================================================================

interface

    ! init topology of the simulated systems
    subroutine cpmf_mask_topo_init(ret_st,natoms,nresidues,has_box,cbox_x,cbox_y,cbox_z)
        integer         :: ret_st
        integer         :: natoms
        integer         :: nresidues
        integer         :: has_box
        real(8)         :: cbox_x
        real(8)         :: cbox_y
        real(8)         :: cbox_z
    end subroutine cpmf_mask_topo_init

    ! load mask topology
    subroutine cpmf_mask_topo_load(ret_st,name)
        integer         :: ret_st
        character(*)    :: name
    end subroutine cpmf_mask_topo_load

    ! set residues
    subroutine cpmf_mask_set_topo_residue(ret_st,idx,name,first_atom)
        integer         :: ret_st
        integer         :: idx
        character(*)    :: name
        integer         :: first_atom
    end subroutine cpmf_mask_set_topo_residue

    ! set atoms
    subroutine cpmf_mask_set_topo_atom(ret_st,idx,name,atype)
        integer         :: ret_st
        integer         :: idx
        character(*)    :: name
        character(*)    :: atype
    end subroutine cpmf_mask_set_topo_atom

    ! set atoms
    subroutine cpmf_mask_set_topo_atom_mcrd(ret_st,idx,lmass,x,y,z)
        integer         :: ret_st
        integer         :: idx
        real(8)         :: lmass
        real(8)         :: x
        real(8)         :: y
        real(8)         :: z
    end subroutine cpmf_mask_set_topo_atom_mcrd

    ! finalize topology
    subroutine cpmf_mask_topo_finalize(ret_st)
        integer         :: ret_st
    end subroutine cpmf_mask_topo_finalize

    ! get topology atom
    subroutine cpmf_mask_get_topo_atom(ret_st,idx,name,resid,resname,x,y,z,lmass,stype)
        integer         :: ret_st
        integer         :: idx
        character(*)    :: name
        real(8)         :: lmass
        integer         :: resid
        character(*)    :: resname
        real(8)         :: x
        real(8)         :: y
        real(8)         :: z
        real(8)         :: mass
        character(*)    :: stype
    end subroutine cpmf_mask_get_topo_atom

    ! get number of atoms in mask
    subroutine cpmf_mask_natoms_in_mask(ret_st,mask,natoms)
        integer         :: ret_st
        character(*)    :: mask
        integer         :: natoms
    end subroutine cpmf_mask_natoms_in_mask

    ! set mask
    subroutine cpmf_mask_set_mask(ret_st,mask)
        integer         :: ret_st
        character(*)    :: mask
        integer         :: natoms
    end subroutine cpmf_mask_set_mask

    ! atom selection status
    subroutine cpmf_mask_is_atom_selected(ret_st,indx,sel)
        integer         :: ret_st
        integer         :: indx
        integer         :: sel
    end subroutine cpmf_mask_is_atom_selected

    ! release memory used by mask subsystem
    subroutine cpmf_mask_clear()
    end subroutine cpmf_mask_clear

    ! print errors
    subroutine cpmf_mask_print_errors()
    end subroutine cpmf_mask_print_errors

end interface

contains

!===============================================================================
! Subroutine:  pmf_mask_topo_init
!===============================================================================

subroutine pmf_mask_topo_init(natoms,ares,has_box,cbox_x,cbox_y,cbox_z)

    use pmf_utils

    implicit none
    integer        :: natoms
    integer        :: ares
    integer        :: has_box
    real(8)        :: cbox_x
    real(8)        :: cbox_y
    real(8)        :: cbox_z
    ! -----------------------------------------------
    integer        :: ret_status
    ! --------------------------------------------------------------------------

    call cpmf_mask_topo_init(ret_status,natoms,ares,has_box,cbox_x,cbox_y,cbox_z)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to init topology of mask subsystem!')
    end if

end subroutine pmf_mask_topo_init

!===============================================================================
! Subroutine:  pmf_mask_topo_load
!===============================================================================

subroutine pmf_mask_topo_load(name)

    use pmf_utils

    implicit none
    character(*)    :: name
    ! -----------------------------------------------
    integer         :: ret_st
    ! --------------------------------------------------------------------------

    call cpmf_mask_topo_load(ret_st,name)

    if( ret_st .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to load topology '//trim(name)//' for mask subsystem!')
    end if

end subroutine pmf_mask_topo_load

!===============================================================================
! Subroutine:  pmf_mask_set_topo_residue
!===============================================================================

subroutine pmf_mask_set_topo_residue(idx,name,first_atom)

    use pmf_utils

    implicit none
    integer            :: idx
    character(*)       :: name
    integer            :: first_atom
    ! -----------------------------------------------
    integer            :: ret_status
    ! --------------------------------------------------------------------------

    call cpmf_mask_set_topo_residue(ret_status,idx,name,first_atom)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to set residue of mask subsystem!')
    end if

end subroutine pmf_mask_set_topo_residue

!===============================================================================
! Subroutine:  pmf_mask_set_topo_atom
!===============================================================================

subroutine pmf_mask_set_topo_atom(idx,name,atype)

    use pmf_utils

    implicit none
    integer         :: idx
    character(*)    :: name
    character(*)    :: atype
    ! -----------------------------------------------
    integer            :: ret_status
    ! --------------------------------------------------------------------------

    call cpmf_mask_set_topo_atom(ret_status,idx,name,atype)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to set atom of mask subsystem!')
    end if

end subroutine pmf_mask_set_topo_atom

!===============================================================================
! Subroutine:  pmf_mask_set_topo_atom
!===============================================================================

subroutine pmf_mask_set_topo_atom_mcrd(idx,m,x,y,z)

    use pmf_utils
    use pmf_dat

    implicit none
    integer        :: idx
    real(PMFDP)    :: m
    real(PMFDP)    :: x
    real(PMFDP)    :: y
    real(PMFDP)    :: z
    ! -----------------------------------------------
    integer        :: ret_status
    real(8)        :: tmpmass
    real(8)        :: tmpx
    real(8)        :: tmpy
    real(8)        :: tmpz
    ! --------------------------------------------------------------------------

    tmpmass = m*MassConv
    tmpx = x*LengthConv
    tmpy = y*LengthConv
    tmpz = z*LengthConv

    call cpmf_mask_set_topo_atom_mcrd(ret_status,idx,tmpmass,tmpx,tmpy,tmpz)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to set atom of mask subsystem!')
    end if

end subroutine pmf_mask_set_topo_atom_mcrd

!===============================================================================
! Subroutine:  pmf_mask_topo_finalize
!===============================================================================

subroutine pmf_mask_topo_finalize()

    use pmf_utils

    implicit none
    integer            :: ret_status
    ! --------------------------------------------------------------------------

    call cpmf_mask_topo_finalize(ret_status)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to finalize topology of mask subsystem!')
    end if

end subroutine pmf_mask_topo_finalize

!===============================================================================
! Subroutine:  pmf_mask_topo_print
!===============================================================================

subroutine pmf_mask_topo_print()

    use pmf_utils
    use pmf_dat

    implicit none
    integer        :: id,resid,ret_id
    real(PMFDP)    :: x,y,z,lmass
    character(4)   :: name,resname,stype
    ! --------------------------------------------------------------------------

    write(PMF_OUT,10)
    write(PMF_OUT,20)

    do id=1,fnatoms
        call cpmf_mask_get_topo_atom(ret_id,id,name,resid,resname,x,y,z,lmass,stype)
        write(PMF_OUT,30) id,name,resid,resname,x,y,z,lmass,stype
    end do

    return

10 format("#    ID    Name  ResID  Res     X [A]      Y [A]      Z [A]    Mass  Type")
20 format("# -------- ---- ------- ---- ---------- ---------- ---------- ------ ----")
30 format(2X,I8,1X,A4,1X,I7,1X,A4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F6.2,1X,A4)

end subroutine pmf_mask_topo_print

!===============================================================================
! Function:  pmf_mask_natoms_in_mask
!===============================================================================

integer function pmf_mask_natoms_in_mask(mask)

    use pmf_utils

    implicit none
    character(*)       :: mask
    ! -----------------------------------------------
    integer            :: ret_status
    ! --------------------------------------------------------------------------

    call cpmf_mask_natoms_in_mask(ret_status,mask,pmf_mask_natoms_in_mask)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
             '[PMFLIB] Unable to count number of atoms in the mask '''//trim(mask)//'''!')
    end if

end function pmf_mask_natoms_in_mask

!===============================================================================
! Subroutine:  pmf_mask_set_mask
!===============================================================================

subroutine pmf_mask_set_mask(mask)

    use pmf_utils

    implicit none
    character(*)    :: mask
    ! -----------------------------------------------
    integer            :: ret_status
    ! --------------------------------------------------------------------------

    call cpmf_mask_set_mask(ret_status,mask)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to set mask!')
    end if

end subroutine pmf_mask_set_mask

!===============================================================================
! Subroutine:  pmf_mask_is_atom_selected
!===============================================================================

subroutine pmf_mask_is_atom_selected(indx,sel)

    use pmf_utils

    implicit none
    integer            :: indx
    integer            :: sel
    ! -----------------------------------------------
    integer            :: ret_status
    ! --------------------------------------------------------------------------

    call cpmf_mask_is_atom_selected(ret_status,indx,sel)

    if( ret_status .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
             '[PMFLIB] Unable to get atom selection status!')
    end if

end subroutine pmf_mask_is_atom_selected

!===============================================================================
! Subroutine:  pmf_mask_is_atom_selected
!===============================================================================

subroutine pmf_mask_print()

    use pmf_utils
    use pmf_dat

    implicit none
    integer         :: ret_status
    integer         :: i
    integer         :: sel
    integer         :: resid
    real(PMFDP)     :: x,y,z,lmass
    character(4)    :: name,resname,stype
    ! --------------------------------------------------------------------------

    write(PMF_OUT,10)
    write(PMF_OUT,20)

    do i=1,fnatoms
        call cpmf_mask_is_atom_selected(ret_status,i,sel)

        if( ret_status .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                 '[PMFLIB] Unable to get atom selection status!')
        end if

        if( sel .eq. 1 ) then
            call cpmf_mask_get_topo_atom(ret_status,i,name,resid,resname,x,y,z,lmass,stype)
            if( ret_status .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                     '[PMFLIB] Unable to get topology atom status!')
            end if
            write(PMF_OUT,30) i,name,resid,resname,x,y,z,lmass,stype
        end if

    end do

10 format("      #    ID    Name  ResID  Res     X [A]      Y [A]      Z [A]    Mass  Type")
20 format("      # -------- ---- ------- ---- ---------- ---------- ---------- ------ ----")
30 format(6X,2X,I8,1X,A4,1X,I7,1X,A4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F6.2,1X,A4)

end subroutine pmf_mask_print

!===============================================================================

end module pmf_mask
