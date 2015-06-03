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

module remd_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  remd_output_open
!===============================================================================

subroutine remd_output_open

 use pmf_utils
 use pmf_dat
 use pmf_constants

 implicit none
 ! -----------------------------------------------------------------------------

 call pmf_utils_open(REMD_OUT,fremdout,'R')

 write(REMD_OUT,10)
 write(REMD_OUT,20)
 write(REMD_OUT,30)

 return

10 format('#===============================================================================')
20 format('# Replica Exchange Molecular Dynamics                                           ')
30 format('#===============================================================================')

end subroutine remd_output_open

!===============================================================================
! Subroutine:  remd_output_write_header
!===============================================================================

subroutine remd_output_write_header

 use pmf_constants
 use pmf_dat
 use abf_dat

 implicit none
 ! -----------------------------------------------------------------------------

 write(REMD_OUT,*)
 write(REMD_OUT,10) '#   NSTEP      Temp         Epot    '
 write(REMD_OUT,10) '#----------- ---------- ------------'

 return

10 format(A)

end subroutine remd_output_write_header

!===============================================================================
! Subroutine:  remd_output_write
!===============================================================================

subroutine remd_output_write(temp)

 use pmf_constants
 use pmf_dat
 use remd_dat

 implicit none
 real(PMFDP)    :: temp
 ! -----------------------------------------------------------------------------

 if( fsample .le. 0 ) return
 if( fstep .le. 1 ) return
 if( mod(fstep,fsample) .ne. 1 ) return

 write(REMD_OUT,10) fstep-1,temp,PotEne

 return

 10 format(2X,I10,1X,F10.3,1X,F12.3)

end subroutine remd_output_write

!===============================================================================
! Subroutine:  remd_output_close
!===============================================================================

subroutine remd_output_close

 use pmf_constants
 use pmf_dat

 implicit none
 ! -----------------------------------------------------------------------------

 close(REMD_OUT)

 return

end subroutine remd_output_close

!===============================================================================

end module remd_output
