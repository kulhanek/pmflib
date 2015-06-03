!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module gap_cvs

implicit none
contains

!===============================================================================
! Subroutine:  gap_cvs_reset_cv
!===============================================================================

subroutine gap_cvs_reset_cv(gap_item)

    use gap_dat

    implicit none
    type(CVTypeGAP)        :: gap_item
    ! --------------------------------------------------------------------------

    gap_item%cvindx          = 0   ! CV index

end subroutine gap_cvs_reset_cv

!===============================================================================
! Subroutine:  gap_cvs_read_cv
!===============================================================================

subroutine gap_cvs_read_cv(prm_fin,gap_item)

    use prmfile
    use pmf_dat
    use pmf_cvs
    use pmf_utils
    use gap_dat

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    type(CVTypeGAP)                    :: gap_item
    logical                            :: file_exist
    ! --------------------------------------------------------------------------

    ! load rest of definition
    if( .not. prmfile_get_string_by_key(prm_fin,'gpfile',gap_item%gpfile) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: gpfile is not specified!')
    else
        write(PMF_OUT,100) trim(gap_item%gpfile)
    end if
    inquire(file=gap_item%gpfile, exist=file_exist)
    if( .not. file_exist ) then
        write(PMF_OUT,200) trim(gap_item%gpfile)
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: gpfile is not found!')
    end if

    if( .not. prmfile_get_string_by_key(prm_fin,'groupname',gap_item%groupname) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: groupname is not specified!')
    else
        write(PMF_OUT,110) trim(gap_item%groupname)
    end if

    if( .not. prmfile_get_integer_by_key(prm_fin,'indexid',gap_item%indexid) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: indexid is not specified!')
    else
        write(PMF_OUT,120) gap_item%indexid
    end if

    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',gap_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: min_value is not specified!')
    else
        write(PMF_OUT,130) gap_item%min_value, trim(gap_item%cv%get_ulabel())
    end if
    call gap_item%cv%conv_to_ivalue(gap_item%min_value)

    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',gap_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: max_value is not specified!')
    else
        write(PMF_OUT,140) gap_item%max_value, trim(gap_item%cv%get_ulabel())
    end if
    call gap_item%cv%conv_to_ivalue(gap_item%max_value)

100 format('   ** gpfile             : ',A)
110 format('   ** groupname          : ',A)
120 format('   ** indexid            : ',I10)
130 format('   ** min_value          : ',E16.7,' [',A,']')
140 format('   ** max_value          : ',E16.7,' [',A,']')

200 format('>>> ERROR: file: ',A,' does not exist!')

end subroutine gap_cvs_read_cv

!===============================================================================
! Subroutine:  gap_cvs_cv_info
!===============================================================================

subroutine gap_cvs_cv_info(gap_item)

    use pmf_dat
    use pmf_cvs
    use gap_dat

    implicit none
    type(CVTypeGAP)    :: gap_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,140) trim(gap_item%cv%name)
    write(PMF_OUT,150) trim(gap_item%cv%ctype)
    write(PMF_OUT,160) gap_item%cv%get_rvalue(CVContext%CVsValues(gap_item%cvindx)), &
                    trim(gap_item%cv%get_ulabel())
    write(PMF_OUT,170) trim(gap_item%gpfile)
    write(PMF_OUT,180) gap_item%cv%get_rvalue(gap_item%min_value), &
                    trim(gap_item%cv%get_ulabel()) 
    write(PMF_OUT,190) gap_item%cv%get_rvalue(gap_item%max_value), &
                    trim(gap_item%cv%get_ulabel())

    return

140 format('    ** Name              : ',a)
150 format('    ** Type              : ',a)
160 format('    ** Current value     : ',E16.7,' [',A,']')
170 format('    ** GP file           : ',a)
180 format('    ** Min value         : ',E16.7,' [',A,']')
190 format('    ** Max value         : ',E16.7,' [',A,']')

end subroutine gap_cvs_cv_info

!===============================================================================

end module gap_cvs
