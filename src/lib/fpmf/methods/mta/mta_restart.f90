!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mta_restart

implicit none
contains

!===============================================================================
! Subroutine:  mta_restart_read
!===============================================================================

subroutine mta_restart_read

    use pmf_dat
    use pmf_utils
    use mta_dat
    use mta_accu

    implicit none
    ! --------------------------------------------------------------------------

    ! test if restart file exists
    if( frestart .and. .not. pmf_utils_fexist(fmtarst) ) then
        frestart = .false.
        write(MTA_OUT,10) trim(fmtarst)
    end if

    if( frestart ) then
        write(MTA_OUT,20)
        ! open restart file ----------------------------------------------------
        call pmf_utils_open(MTA_RST,fmtarst,'O')

        call mta_accu_read(MTA_RST)

        close(MTA_RST)
    else
        write(MTA_OUT,30)
    end if

    return

 10 format('# WARNING: frestart = on, but file (',A,') does not exist! => frestart = off')
 20 format('# RST: frestart = on')
 30 format('# RST: frestart = off')

end subroutine mta_restart_read

!===============================================================================
! Subroutine:  mta_restart_update
!===============================================================================

subroutine mta_restart_update

    use pmf_dat
    use pmf_utils
    use mta_accu
    use mta_dat

    implicit none
    !---------------------------------------------------------------------------

    if( frstupdate .le. 0 ) return ! trajectory is written only of frstupdate > 0

    if( mod(fstep,frstupdate) .ne. 0 ) return

    call pmf_utils_open(MTA_RST,fmtarst,'U')
    call mta_accu_write(MTA_RST)
    close(MTA_RST)

    write(MTA_OUT,10) fstep, insidesamples, outsidesamples

    return

 10 format('# [ACCU] Total steps     = ',I12,' Inside samples  = ',I12,' Outside samples = ',I12  )

end subroutine mta_restart_update

!===============================================================================
! Subroutine:  mta_restart_write
!===============================================================================

subroutine mta_restart_write()

    use pmf_dat
    use pmf_utils
    use mta_accu

    implicit none
    !---------------------------------------------------------------------------

    call pmf_utils_open(MTA_RST,fmtarst,'U')

    call mta_accu_write(MTA_RST)

    close(MTA_RST)

    return

end subroutine mta_restart_write

!===============================================================================

end module mta_restart

