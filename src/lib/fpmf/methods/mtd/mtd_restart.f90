!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module mtd_restart

implicit none
contains

!===============================================================================
! Subroutine:  mtd_restart_read
!===============================================================================

subroutine mtd_restart_read

    use pmf_dat
    use pmf_utils
    use mtd_dat
    use mtd_accu

    implicit none
    ! --------------------------------------------------------------------------

    ! test if restart file exists
    if( frestart .and. .not. pmf_utils_fexist(fmtdrst) ) then
        frestart = .false.
        write(MTD_OUT,10) trim(fmtdrst)
    end if

    if( frestart ) then
        write(MTD_OUT,20)
        ! open restart file ----------------------------------------------------
        call pmf_utils_open(MTD_RST,fmtdrst,'O')

        call mtd_accu_read(MTD_RST)

        close(MTD_RST)
    else
        write(MTD_OUT,30)
    end if

    return

 10 format('# WARNING: frestart = on, but file (',A,') does not exist! => frestart = off')
 20 format('# RST: frestart = on')
 30 format('# RST: frestart = off')

end subroutine mtd_restart_read

!===============================================================================
! Subroutine:  mtd_restart_update
!===============================================================================

subroutine mtd_restart_update

    use pmf_dat
    use pmf_utils
    use mtd_accu
    use mtd_dat

    implicit none
    !---------------------------------------------------------------------------

    if( frstupdate .le. 0 ) return ! trajectory is written only of frstupdate > 0

    if( mod(fstep,frstupdate) .ne. 0 ) return

    call pmf_utils_open(MTD_RST,fmtdrst,'U')
    call mtd_accu_write(MTD_RST)
    close(MTD_RST)

    write(MTD_OUT,10) fstep, insidesamples, outsidesamples

    return

 10 format('# [ACCU] Total steps     = ',I12,' Inside gaussians  = ',I12,' Outside gaussians = ',I12  )

end subroutine mtd_restart_update

!===============================================================================
! Subroutine:  mtd_restart_write
!===============================================================================

subroutine mtd_restart_write

    use pmf_dat
    use pmf_utils
    use mtd_accu

    implicit none
    !---------------------------------------------------------------------------

    call pmf_utils_open(MTD_RST,fmtdrst,'U')

    call mtd_accu_write(MTD_RST)

    close(MTD_RST)

    return

end subroutine mtd_restart_write

!===============================================================================

end module mtd_restart

