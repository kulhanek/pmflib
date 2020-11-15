!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abp_restart

implicit none
contains

!===============================================================================
! Subroutine:  abp_restart_read
!===============================================================================

subroutine abp_restart_read

    use pmf_dat
    use pmf_utils
    use abp_dat
    use abp_accumulator

    implicit none
    ! --------------------------------------------------------------------------

    ! test if restart file exists
    if( frestart .and. .not. pmf_utils_fexist(fabprst) ) then
        frestart = .false.
        write(ABP_OUT,10) trim(fabprst)
    end if

    if( frestart ) then
        write(ABP_OUT,20)
        ! open restart file ----------------------------------------------------
        call pmf_utils_open(ABP_RST,fabprst,'O')

        call abp_accumulator_read(ABP_RST)

        close(ABP_RST)
    else
        write(ABP_OUT,30)
    end if

    return

    10 format('# WARNING: frestart = on, but file (',A,') does not exist! => frestart = off')
    20 format('# RST: frestart = on')
    30 format('# RST: frestart = off')

end subroutine abp_restart_read

!===============================================================================
! Subroutine:  abp_restart_update
!===============================================================================

subroutine abp_restart_update

    use pmf_dat
    use pmf_utils
    use abp_accumulator
    use abp_dat

    implicit none
    !---------------------------------------------------------------------------

    if( frstupdate .le. 0 ) return ! trajectory is written only of frstupdate > 0

    if( mod(fstep,frstupdate) .ne. 0 ) return

    call pmf_utils_open(ABP_RST,fabprst,'U')

    call abp_accumulator_write(ABP_RST)

    close(ABP_RST)

    return

end subroutine abp_restart_update

!===============================================================================
! Subroutine:  abp_restart_write
!===============================================================================

subroutine abp_restart_write

    use pmf_dat
    use pmf_utils
    use abp_accumulator

    implicit none
    !---------------------------------------------------------------------------

    call pmf_utils_open(ABP_RST,fabprst,'U')

    call abp_accumulator_write(ABP_RST)

    close(ABP_RST)

    return

end subroutine abp_restart_write

!===============================================================================

end module abp_restart

