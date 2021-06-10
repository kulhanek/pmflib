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
 use mtd_grid

 implicit none
 ! -----------------------------------------------------------------------------

 call mtd_init_print_header
 call mtd_init_core
 call mtd_grid_init
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
 fplevel          = 0         ! output print level
 fpsample         = 500       ! output print sampling
 fmetastep        = 1000      ! step size of mtddynamics
 fheight          = 0.10d0    ! height of gaussian
 fmetatemp        = 0         ! temperature for well-tempered mtd
 fmetavary        = 0         ! what to vary when well-tempered mtd is on
                              ! 0 - nothing, 1 - height, 2 - step
 frestart         = .false.   ! restart job with previous data
 fextout          = 0         ! control extended output (cvs,hills)
 fbuffersize      = 1000      ! buffer size
 fscaling         = 0         ! scaling type: 0 - none, 1 - step signal, 2 - ramp signal
 fdelaying        = 0         ! delaying time in md steps
 ftransition      = 0         ! transition time in md steps 

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
 use mtd_cvs_mod
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
 write(PMF_OUT,135)  ' Meta temperature (fmetatemp)            : ', fmetatemp
 write(PMF_OUT,130)  ' Varying quantity (fmetavary)            : ', fmetavary
 select case(fmetavary)
     case(0)
        write(PMF_OUT,120)  '  == varying quantity is not applied'
     case(1)
        write(PMF_OUT,120)  '  == varying quantity is fheight'
     case(2)
        write(PMF_OUT,120)  '  == varying quantity is fmetastep'
 end select
 write(PMF_OUT,130)  ' Scaling type (fscaling)                 : ', fscaling
 select case(fscaling)
     case(0)
        write(PMF_OUT,120)  '  == scaling is not applied'
     case(1)
        write(PMF_OUT,120)  '  == scaling type is step signal'
     case(2)
        write(PMF_OUT,120)  '  == scaling type is ramp signal'
 end select
 write(PMF_OUT,130)  ' Transition step size (ftransition)      : ', ftransition
 write(PMF_OUT,130)  ' Delaying step size (fdelaying)          : ', fdelaying
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
 write(PMF_OUT,130)  ' Output sampling (fpsample)              : ', fpsample
 write(PMF_OUT,130)  ' Print level (fplevel)                   : ', fplevel
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
 write(PMF_OUT,120)  ' List of metadynamics coordinates'
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

 use pmf_dat
 use pmf_utils
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
