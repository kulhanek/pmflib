!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module abf_cvs

implicit none
contains

!===============================================================================
! Subroutine:  abf_cvs_reset_cv
!===============================================================================

subroutine abf_cvs_reset_cv(abf_item)

    use abf_dat

    implicit none
    type(CVTypeABF) :: abf_item
    ! --------------------------------------------------------------------------

    abf_item%cvindx         = 0     ! CV index
    abf_item%cv             => null()
    abf_item%set            = 0
    abf_item%min_value      = 0.0   ! left range
    abf_item%max_value      = 0.0   ! right range
    abf_item%nbins          = 0     ! number of bins

end subroutine abf_cvs_reset_cv

!===============================================================================
! Subroutine:  abf_cvs_read_cv
!===============================================================================

subroutine abf_cvs_read_cv(prm_fin,abf_item)

    use prmfile
    use abf_dat
    use pmf_cvs
    use pmf_unit
    use pmf_utils
    use pmf_paths
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeABF)                     :: abf_item
    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( abf_item%cv%pathidx .gt. 0 ) then
        if( PathList(abf_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in a driven mode!')
        end if
    end if

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',abf_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) abf_item%min_value, trim(abf_item%cv%get_ulabel())
    call abf_item%cv%conv_to_ivalue(abf_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',abf_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) abf_item%max_value, trim(abf_item%cv%get_ulabel())
    call abf_item%cv%conv_to_ivalue(abf_item%max_value)

    if( abf_item%max_value .le. abf_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',abf_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    end if
    if( abf_item%nbins .lt. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins has to be greater then zero!')
    end if
    write(PMF_OUT,125) abf_item%nbins

    select case(feimode)
        case(0,1)
            ! nothing to be here
        !------------------------------
        case(2)
            if( .not. prmfile_get_integer_by_key(prm_fin,'wfac',abf_item%wfac) ) then
                call pmf_utils_exit(PMF_OUT,1,'wfac is not specified!')
            end if
            call pmf_ctrl_check_integer('ABF','wfac',abf_item%wfac,0,CND_GT)
            write(PMF_OUT,170) abf_item%wfac
        !------------------------------
        case default
            ! nothing to be here
    end select


    return

110 format('    ** Min value         : ',F16.7,' [',A,']')
120 format('    ** Max value         : ',F16.7,' [',A,']')
125 format('    ** Number of bins    : ',I8)
170 format('    ** Smoothing W-fac   : ',I8)

end subroutine abf_cvs_read_cv

!===============================================================================
! Subroutine:  abf_cvs_cv_info
!===============================================================================

subroutine abf_cvs_cv_info(abf_item)

    use abf_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit

    implicit none
    type(CVTypeABF) :: abf_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,145) trim(abf_item%cv%name)
    write(PMF_OUT,146) trim(abf_item%cv%ctype)
    write(PMF_OUT,150) abf_item%cv%get_rvalue(CVContext%CVsValues(abf_item%cvindx)), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,155) abf_item%cv%get_rvalue(abf_item%min_value), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,160) abf_item%cv%get_rvalue(abf_item%max_value), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,165) abf_item%nbins

    if( feimode .eq. 2 ) then
    write(PMF_OUT,170) abf_item%wfac
    end if

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)
170 format('    ** GKS Wfactor       : ',F16.7)

end subroutine abf_cvs_cv_info

!===============================================================================

end module abf_cvs
