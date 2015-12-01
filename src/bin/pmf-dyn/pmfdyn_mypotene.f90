! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License along
!     with this program; if not, write to the Free Software Foundation, Inc.,
!     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ==============================================================================

module pmfdyn_mypotene

contains

! ------------------------------------------------------------------------------

subroutine pot_ext_energy_default(x,d,ene)

    use smf_profiling
    use smf_xyzfile_type
    use smf_xyzfile
    use pmfdyn_system_dat
    use pmf_utils
    ! use ifport

    implicit none
    real(PMFDP)            :: x(:,:)
    real(PMFDP)            :: d(:,:)
    real(PMFDP)            :: ene
    ! -----------------------------------------------
    type(XYZFILE_TYPE)     :: io_xyz
    real(PMFDP)            :: edum
    ! --------------------------------------------------------------------------

    call start_timer(EXTERNAL_TIMER)
    ! -----------------------------------------------
    ! allocate arrays
    call allocate_xyz(io_xyz,natoms)

    ! open xyz file
    call open_xyz(IO_FRC,'input.xyz',io_xyz,'UNKNOWN')

    ! copy symbols and coordinates
    io_xyz%comment = 'cvs'
    io_xyz%symbols(:) = symbols(:)
    io_xyz%cvs(:,:) = x(:,:)

    ! write coordinates
    call write_xyz(IO_FRC,io_xyz)

    ! close input cordinates
    call close_xyz(IO_FRC,io_xyz)

    ! -----------------------------------------------
    !iresult = system('./calc-eg')
    call system('./calc-eg')

    ! if( iresult .ne. 0 ) then
    !    call pmf_utils_exit(PMF_OUT,1,'Unable to calculate external forces!')
    ! end if

    ! -----------------------------------------------
    ! init xyz file structure
    call init_xyz(io_xyz)

    ! open xyz file
    call open_xyz(IO_FRC,'result.xyz',io_xyz,'OLD')

    ! read gradients
    call read_xyz(IO_FRC,io_xyz)

    ! close gradient file
    call close_xyz(IO_FRC,io_xyz)

    ! decode result data
    d(:,:) = io_xyz%cvs(:,:)
    read(io_xyz%comment,*) edum
    ene = edum

    call stop_timer(EXTERNAL_TIMER)

    return

end subroutine pot_ext_energy_default

! ------------------------------------------------------------------------------

end module pmfdyn_mypotene

!===============================================================================

