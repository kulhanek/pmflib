!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cst_restart

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_restart_read
!===============================================================================

subroutine cst_restart_read

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_accu

    implicit none
    ! --------------------------------------------------------------------------

    if( frestart .and. (.not. fallconstant) ) then
        frestart = .false.
        write(CST_OUT,5) trim(fcstrst)
    end if

    ! test if restart file exists
    if( frestart .and. .not. pmf_utils_fexist(fcstrst) ) then
        frestart = .false.
        write(CST_OUT,10) trim(fcstrst)
    end if

    if( frestart ) then
        write(CST_OUT,20)
        ! open restart file ----------------------------------------------------
        call pmf_utils_open(CST_RST,fcstrst,'O')
        call cst_accu_read(CST_RST)
        close(CST_RST)
    else
        write(CST_OUT,30)
    end if

    return

  5 format('# WARNING: frestart = on, not all CVs are constant  => frestart = off')
 10 format('# WARNING: frestart = on, but file (',A,') does not exist! => frestart = off')
 20 format('# RST: frestart = on')
 30 format('# RST: frestart = off')

end subroutine cst_restart_read

!===============================================================================
! Subroutine:  cst_restart_update
!===============================================================================

subroutine cst_restart_update

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_accu

    implicit none
    !---------------------------------------------------------------------------

    if( .not. fallconstant ) return     ! all must be constant
    if( frstupdate .le. 0 ) return      ! restart is updated only of frstupdate > 0

    if( mod(fstep,frstupdate) .ne. 0 ) return

    call pmf_utils_open(CST_RST,fcstrst,'U')
    call cst_accu_write(CST_RST)
    close(CST_RST)

    write(CST_OUT,10) fstep

    return

 10 format('# [ACCU] Total steps     = ',I12)

end subroutine cst_restart_update

!===============================================================================
! Subroutine:  cst_restart_write
!===============================================================================

subroutine cst_restart_write

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_accu

    implicit none
    !---------------------------------------------------------------------------

    if( .not. fallconstant ) return

    call pmf_utils_open(CST_RST,fcstrst,'U')
    call cst_accu_write(CST_RST)
    close(CST_RST)

    return

end subroutine cst_restart_write

!===============================================================================

end module cst_restart

