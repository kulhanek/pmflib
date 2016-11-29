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

module rst_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mon_init_method
!===============================================================================

subroutine rst_init_method

    use rst_output
    use rst_restart
    use rst_accumulator

    implicit none
    ! --------------------------------------------------------------------------

    call rst_init_core
    call rst_accumulator_init
    call rst_init_print_header
    call rst_output_open
    call rst_restart_read
    call rst_output_write_header

end subroutine rst_init_method

!===============================================================================
! Subroutine:  rst_init_dat
!===============================================================================

subroutine rst_init_dat

    use rst_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode        = 0         ! 0 - disable RST, 1 - enabled RST
    fsample      = 500       ! output sample period in steps
    fplevel      = 0         ! print level
    frestart     = .false.   ! restart job with previous data
    fhistupdate   = 0        ! how often is restart file written
    fhistclear     = 0       ! after first 'fhistclear' steps histogram
                             ! will be reset (default 0)
    fsamplefreq   = 1        ! how often take samples
    faccumulation = 0        ! number of accumulated data

    NumOfRSTItems = 0
    insidesamples = 0
    outsidesamples = 0

end subroutine rst_init_dat

!===============================================================================
! Subroutine:  rst_init_print_header
!===============================================================================

subroutine rst_init_print_header

    use prmfile
    use pmf_dat
    use rst_dat
    use rst_restraints

    implicit none
    integer :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' ***************************** RESTRAINED DYNAMICS **************************** '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restrained Dynamics Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Restrained dynamics mode (fmode)        : ', fmode
    write(PMF_OUT,130)  ' Number of restraints                    : ', NumOfRSTItems
    write(PMF_OUT,125)  ' Restraint definition file (frstdef)     : ', trim(frstdef)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (frstout)                   : ', trim(frstout)
    write(PMF_OUT,130)  ' Sample period (fsample)                 : ', fsample
    write(PMF_OUT,130)  ' Print level (fplevel)                   : ', fplevel

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Histogram options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Histogram file (frsthist)               : ', trim(frsthist)
    write(PMF_OUT,125)  ' Histogram restart enabled (frestart)    : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Histogram file update (fhistupdate)     : ', fhistupdate
    write(PMF_OUT,130)  ' Histogram sampling restart (fhistclear) : ', fhistclear
    write(PMF_OUT,120)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of restraints'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfRSTItems
     write(PMF_OUT,140) i
        call rst_restraints_rst_info(RSTCVList(i))
        write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

    return

120 format(A)
125 format(A,A)
130 format(A,I6)
135 format(A,E12.5)

140 format(' == Collective variable #',I4.4)

end subroutine rst_init_print_header

!===============================================================================
! Subroutine:  rst_init_core
!===============================================================================

subroutine rst_init_core

 use pmf_utils
 use pmf_dat
 use rst_dat

 implicit none
 integer      :: alloc_failed

 ! ------------------------------------------------------------------------------

 ! allocate arrays for accumulating data

 allocate(svalues(NumOfRSTItems), &
          s2values(NumOfRSTItems), &
          stat=alloc_failed)

 if( alloc_failed .ne. 0 ) then
    call pmf_utils_exit(PMF_OUT,1,'[RST] Unable to allocate memory for arrays of accumutalion!')
 end if

end subroutine rst_init_core

!===============================================================================

end module rst_init
