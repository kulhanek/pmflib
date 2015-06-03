!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module remd_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  remd_init_method
!===============================================================================

subroutine remd_init_method

 use remd_output
 use remd_client
 use pmf_dat
 use remd_dat

 implicit none
 ! -----------------------------------------------------------------------------

 OldBathTemp = ftemp

 call remd_init_print_header
 call remd_output_open
 call remd_client_register
 call remd_client_get_initial_data
 call remd_output_write_header

end subroutine remd_init_method

!===============================================================================
! Subroutine:  remd_init_dat
!===============================================================================

subroutine remd_init_dat

 use remd_dat

 implicit none
 ! -----------------------------------------------------------------------------

 fmode              = 0
 fsample            = 500
 fserverkey         = ''
 fserver            = ''
 fpassword          = ''

 use_key = .false.
 ReplicaId = 0
 BathId = 0

 REMDFullSwap       = .false.
 REMDPeriod         = 5000
 OldBathTemp        = 300.0
 CurBathTemp        = 300.0

end subroutine remd_init_dat

!===============================================================================
! Subroutine:  remd_init_print_header
!===============================================================================

subroutine remd_init_print_header

 use pmf_constants
 use remd_dat
 use pmf_dat

 implicit none
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,120)
 write(PMF_OUT,120)  '================================================================================'
 write(PMF_OUT,120)  ' ******************** REPLICA EXCHANGE MOLECULAR DYNAMICS ********************* '
 write(PMF_OUT,120)  '================================================================================'
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' REMD Options'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,130)  ' REMD mode (fmode)                       : ', fmode
 write(PMF_OUT,130)  ' REMD print frequency (fsample)          : ', fsample
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' REMD server options:'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 if( use_key ) then
 write(PMF_OUT,125)  ' Server Key file name (fserverkey)       : ', trim(fserverkey)
 else
 write(PMF_OUT,125)  ' Server URL (fserver)                    : ', trim(fserver)
 write(PMF_OUT,125)  ' Server password (fpassword)             : ', trim(fpassword)
 end if
 write(PMF_OUT,120)
 write(PMF_OUT,120)  ' Output options:'
 write(PMF_OUT,120)  ' ------------------------------------------------------'
 write(PMF_OUT,125)  ' Output file (fremdout)                  : ', trim(fremdout)
 write(PMF_OUT,120)
 write(PMF_OUT,120)  '================================================================================'

 return

120 format(A)
125 format(A,A)
130 format(A,I6)
135 format(A,E12.5)

end subroutine remd_init_print_header

!===============================================================================

end module remd_init
