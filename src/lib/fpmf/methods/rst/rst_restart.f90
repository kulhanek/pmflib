!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module rst_restart

implicit none
contains

!===============================================================================
! Subroutine:  rst_restart_read
!===============================================================================

subroutine rst_restart_read

    use pmf_dat
    use pmf_utils
    use rst_dat
    use rst_accumulator

    implicit none
    ! --------------------------------------------------------------------------

    ! test if restart file exists
    if( frestart .and. .not. pmf_utils_fexist(frsthist) ) then
        frestart = .false.
        write(RST_OUT,10) trim(frsthist)
    end if

    if( frestart ) then
        write(RST_OUT,20)
        ! open restart file -----------------------------------------------------------
        call pmf_utils_open(RST_RST,frsthist,'O')

        call rst_accumulator_read(RST_RST)

        close(RST_RST)
    else
        write(RST_OUT,30)
    end if

    return

    10 format('# WARNING: frestart = on, but file (',A,') does not exist! => frestart = off')
    20 format('# RST: frestart = on')
    30 format('# RST: frestart = off')

end subroutine rst_restart_read

!===============================================================================
! Subroutine:  rst_restart_update
!===============================================================================

subroutine rst_restart_update

    use pmf_dat
    use pmf_utils
    use rst_accumulator
    use rst_dat

    implicit none
    !---------------------------------------------------------------------------

    ! only if we would like to update restart file or if restart is explicitly required
    if( (fhistupdate .eq. 0) .and. (frestart .eqv. .false.) ) return

    if( fhistupdate .le. 0 ) return ! trajectory is written only of fhistupdate > 0

    if( mod(fstep,fhistupdate) .ne. 0 ) return

    call pmf_utils_open(RST_RST,frsthist,'U')

    call rst_accumulator_write(RST_RST)

    close(RST_RST)

    return

end subroutine rst_restart_update

!===============================================================================
! Subroutine:  rst_restart_write
!===============================================================================

subroutine rst_restart_write

    use pmf_dat
    use rst_dat
    use pmf_utils
    use rst_accumulator

    implicit none
    !---------------------------------------------------------------------------

    ! only if we would like to update restart file or if restart is explicitly required
    if( (fhistupdate .eq. 0) .and. (frestart .eqv. .false.) ) return

    call pmf_utils_open(RST_RST,frsthist,'U')

    call rst_accumulator_write(RST_RST)

    close(RST_RST)

    return

end subroutine rst_restart_write

!===============================================================================

end module rst_restart

