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

module mtd_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_init_method
!===============================================================================

subroutine mtd_init_method

 use mtd_output
 use mtd_restart
 use mtd_client

 implicit none
 ! -----------------------------------------------------------------------------

 call mtd_init_print_header
 call mtd_init_core
 call mtd_output_open
 call mtd_restart_open_and_read
 call mtd_client_register
 call mtd_client_get_initial_data
 call mtd_output_write_header

end subroutine mtd_init_method

!===============================================================================
! Subroutine:  mtd_init_dat
!===============================================================================

subroutine mtd_init_dat

 use mtd_dat

 implicit none
 ! -----------------------------------------------------------------------------

! control section --------------------------------------------------------------
 fmode            = 0         ! mode of mtddynamics
 fsample          = 500       ! output print sampling
 fmetastep        = 1000      ! step size of mtddynamics
 fheight          = 0.10d0    ! height of gaussian
 frestart         = .false.   ! restart job with previous data
 fextout          = 0         ! control extended output (cvs,hills)
 fbuffersize      = 1000      ! buffer size

! server part ------------------------------------------------------------------
 fserver_enabled        = .false.   ! is mtddyn-server enabled?
 fserver_key_enabled    = .false.
 fclient_id             = 0         ! client id
 fserverkey             = ''        ! mtddyn-server name
 fserver                = ''        ! mtddyn-server name
 fpassword              = ''        ! mtddyn-server password

 NumOfMTDCVs     = 0 

end subroutine mtd_init_dat

!===============================================================================
! Subroutine:  mtd_init_print_header
!===============================================================================

subroutine mtd_init_print_header

 use mtd_dat
 use pmf_dat
 use pmf_constants
 use mtd_cvs
 use prmfile

 implicit none
 integer        :: i
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,120) 
 write(PMF_OUT,120)  '================================================================================'
 write(PMF_OUT,120)  ' ******************************* METADYNAMICS ********************************* '
 write(PMF_OUT,120)  '================================================================================' 
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' Metadynamics Mode'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,130)  ' Metadynamics mode (fmode)               : ', fmode
 write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfMTDCVs
 write(PMF_OUT,125)  ' Coordinate definition file (fmtddef)    : ', trim(fmtddef)
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' Metadynamics Controls:'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,130)  ' Metadynamics step size (fmetastep)      : ', fmetastep
 write(PMF_OUT,135)  ' Height (fheight) [kcal/mol]             : ', fheight
 write(PMF_OUT,130)  ' Buffer size (fbuffersize)               : ', fbuffersize
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' Restart options:'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,125)  ' Restart file (fmtdrst)                  : ', trim(fmtdrst)
 write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' Output options:'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,125)  ' Output file (fmtdout)                   : ', trim(fmtdout)
 write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' Extended output options:'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,130)  ' Extended output mode (fextout)          : ', fextout
 write(PMF_OUT,125)  ' CVs file (fmtdcvs)                      : ', trim(fmtdcvs)
 write(PMF_OUT,125)  ' Hills file (fmtdhills)                  : ', trim(fmtdhills)
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' MTD server options:'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,125)  ' Server communication is                 : ', prmfile_onoff(fserver_enabled)
 if( fserver_enabled ) then
 if( fserver_key_enabled ) then
 write(PMF_OUT,125)  ' Server Key file name (fserverkey)       : ', trim(fserverkey)
 else
 write(PMF_OUT,125)  ' Server URL (fserver)                    : ', trim(fserver)
 write(PMF_OUT,125)  ' Server password (fpassword)             : ', trim(fpassword)
 end if
 else
 write(PMF_OUT,125)  ' Server URL (fserver)                    : ', 'none'
 write(PMF_OUT,125)  ' Server password (fpassword)             : ', 'none'
 end if
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' List of mtddynamics coordinates'
 write(PMF_OUT,120)  ' -------------------------------------------------------'
 write(PMF_OUT,120)

 do i=1,NumOfMTDCVs
    write(PMF_OUT,140) i
    call mtd_cvs_cv_info(MTDCVList(i))
    write(PMF_OUT,120)
 enddo

 write(PMF_OUT,120)  '================================================================================'

 return

120 format(A)
125 format(A,A)
130 format(A,I6)
135 format(A,E12.5)

140 format(' == Collective variable #',I4.4)

end subroutine mtd_init_print_header

!===============================================================================
! Subroutine:  mtd_init_core
!===============================================================================

subroutine mtd_init_core

 use mtd_dat
 use mtd_history

 implicit none
 ! -----------------------------------------------------------------------------

 ! create first buffer ---------------------------
 if( .not. fserver_enabled ) then
    ! in the case of multiple-walker approach, initial allocation is done
    ! in mtd_get_initial_data
    call mtd_history_allocate_buffer(hill_history,fbuffersize)
 end if

 return

end subroutine mtd_init_core

!===============================================================================

end module mtd_init
