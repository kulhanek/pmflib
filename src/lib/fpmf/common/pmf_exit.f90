!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module pmf_exit

implicit none
contains

!===============================================================================
! Subroutine:   pmf_exit_mdloop(exitcode)
!===============================================================================

subroutine pmf_exit_mdloop(unitnum,errcode,message)

    use pmf_dat
    use pmf_utils

    implicit none
    integer                :: unitnum
    integer                :: errcode
    character(*),optional  :: message
    ! --------------------------------------------------------------------------

    ! fexit_mdloop must be broadcasted to other MPI processes
    ! most likely in force_mpi subroutines
    fexit_mdloop = errcode

    if( .not. fcanexmdloop ) then
        if( present(message) ) then
            call pmf_utils_exit(unitnum,fexit_mdloop,message)
        else
            call pmf_utils_exit(unitnum,fexit_mdloop,'Client was terminated as it does not know how to terminate MD loop')
        end if
    else
        if( present(message) ) then
            write(unitnum,'(/,A)') '>>> MD LOOP TERMINATED: ' // trim(message)
        end if
    end if

end subroutine pmf_exit_mdloop

!===============================================================================

end module pmf_exit
