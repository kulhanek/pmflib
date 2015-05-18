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

module abf_restart

implicit none
contains

!===============================================================================
! Subroutine:  abf_restart_read
!===============================================================================

subroutine abf_restart_read

    use pmf_dat
    use pmf_utils
    use abf_dat
    use abf_accumulator

    implicit none
    ! --------------------------------------------------------------------------

    ! read mask
    if( fmask_mode .eq. 1 ) then
        write(ABF_OUT,5) trim(fabfmask)
        call pmf_utils_open(ABF_RST,fabfmask,'O')
        call abf_accumulator_read_mask(ABF_RST)
        close(ABF_RST)
    end if

    ! test if restart file exists
    if( frestart .and. .not. pmf_utils_fexist(fabfrst) ) then
        frestart = .false.
        write(ABF_OUT,10) trim(fabfrst)
    end if

    if( frestart ) then
        write(ABF_OUT,20)
        ! open restart file ----------------------------------------------------
        call pmf_utils_open(ABF_RST,fabfrst,'O')

        call abf_accumulator_read(ABF_RST)

        close(ABF_RST)
    else
        write(ABF_OUT,30)
    end if

    ! apply mask to read restart file
    if( (fmask_mode .eq. 1) .and. (.not. fserver_enabled) ) then
        write(ABF_OUT,40)
        call abf_accumulator_apply_mask()
    end if

    return

  5 format('# MASK: reading mask: ',A)
 10 format('# WARNING: frestart = on, but file (',A,') does not exist! => frestart = off')
 20 format('# RST: frestart = on')
 30 format('# RST: frestart = off')
 40 format('# MASK: applying mask to current ABF accumulator')

end subroutine abf_restart_read

!===============================================================================
! Subroutine:  abf_restart_update
!===============================================================================

subroutine abf_restart_update

    use pmf_dat
    use pmf_utils
    use abf_accumulator
    use abf_dat

    implicit none
    !---------------------------------------------------------------------------

    if( frstupdate .le. 0 ) return ! trajectory is written only of frstupdate > 0

    if( mod(fstep,frstupdate) .ne. 0 ) return

    call pmf_utils_open(ABF_RST,fabfrst,'U')
    call abf_accumulator_write(ABF_RST)
    close(ABF_RST)

    write(ABF_OUT,10) fstep, insidesamples, outsidesamples

    return

 10 format('# [ACCU] Total steps     = ',I12,' Inside samples  = ',I12,' Outside samples = ',I12  )

end subroutine abf_restart_update

!===============================================================================
! Subroutine:  abf_restart_write
!===============================================================================

subroutine abf_restart_write

    use pmf_dat
    use pmf_utils
    use abf_accumulator

    implicit none
    !---------------------------------------------------------------------------

    call pmf_utils_open(ABF_RST,fabfrst,'U')

    call abf_accumulator_write(ABF_RST)

    close(ABF_RST)

    return

end subroutine abf_restart_write

!===============================================================================

end module abf_restart

