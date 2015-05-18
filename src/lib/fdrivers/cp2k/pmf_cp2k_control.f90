!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2010,2011 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
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

module pmf_cp2k_control

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_cp2k_process_control

 use pmf_constants
 use pmf_cp2k_dat
 use prmfile
 use pmf_utils
 use pmf_control
 use pmf_dat

 implicit none
 ! -----------------------------------------------------------------------------

 ! write header
 write(PMF_OUT,*)
 call pmf_utils_heading(PMF_OUT,'Reading control file','-')
 write(PMF_OUT,'(a,a)') 'Control file name : ',trim(ControlFileName)

 ! open file
 if( .not. prmfile_read(ControlPrmfile,ControlFileName) ) then
    call pmf_utils_exit(PMF_OUT,1,'Specified control file cannot be opened!')
 end if

 ! read groups
 call pmf_control_read_pmflib_group(ControlPrmfile)

 ! write footer
 write(PMF_OUT,*)
 call pmf_utils_heading(PMF_OUT,'{END}',':')
 write(PMF_OUT,*)

 ! now we check if everything was understood from control file
 if( prmfile_count_ulines(ControlPrmfile) .gt. 0 ) then
    write(PMF_OUT,'(/,a,/)') '=== [unprocessed options] ========================================'
    call prmfile_dump(ControlPrmfile,PMF_OUT,.true.)
    call pmf_utils_exit(PMF_OUT,1,'Some items in the control file were not understood (see above)!')
 end if

 return

end subroutine pmf_cp2k_process_control

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmf_cp2k_control
