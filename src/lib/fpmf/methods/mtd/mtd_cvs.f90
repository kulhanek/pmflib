!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

! ifort version 17.0.1
! mtd_init.f90(117): error #6406: Conflicting attributes or multiple declaration of name.   [MTD_CVS]
! renamed mtd_cvs -> mtd_cvs_mod

module mtd_cvs_mod

implicit none
contains

!===============================================================================
! Subroutine:  mtd_cvs_reset_cv
!===============================================================================

subroutine mtd_cvs_reset_cv(mtd_item)

    use mtd_dat

    implicit none
    type(CVTypeMTD)        :: mtd_item
    ! --------------------------------------------------------------------------

    mtd_item%cvindx         = 0     ! CV index
    mtd_item%width          = 0.0   ! width
    mtd_item%min_value      = 0.0   ! left range
    mtd_item%max_value      = 0.0   ! right range
    mtd_item%nbins          = 0     ! number of bins

end subroutine mtd_cvs_reset_cv

!===============================================================================
! Subroutine:  mtd_cvs_read_cv
!===============================================================================

subroutine mtd_cvs_read_cv(prm_fin,mtd_item)

    use prmfile
    use mtd_dat
    use pmf_dat
    use pmf_cvs
    use pmf_utils
    use pmf_paths
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    type(CVTypeMTD)                    :: mtd_item
    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( mtd_item%cv%pathidx .gt. 0 ) then
        if( PathList(mtd_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in a driven mode!')
        end if
    end if

    ! load rest of definition
    if( .not. prmfile_get_real8_by_key(prm_fin,'width',mtd_item%width) ) then
        call pmf_utils_exit(PMF_OUT,1,'width is not specified!')
    end if
    write(PMF_OUT,100) mtd_item%width, trim(mtd_item%cv%get_ulabel())
    call pmf_ctrl_check_real8_wunit('MTD','alpha',mtd_item%cv%unit,mtd_item%width,0.0d0,CND_GT,'f12.2')
    call mtd_item%cv%conv_to_ivalue(mtd_item%width)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',mtd_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) mtd_item%min_value, trim(mtd_item%cv%get_ulabel())
    call mtd_item%cv%conv_to_ivalue(mtd_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_deposit',mtd_item%min_deposit) ) then
        mtd_item%min_deposit = mtd_item%min_value    ! already in internal units
        write(PMF_OUT,115) mtd_item%cv%get_rvalue(mtd_item%min_deposit), trim(mtd_item%cv%get_ulabel())
        if( mtd_item%min_deposit .lt. mtd_item%min_value ) then
            call pmf_utils_exit(PMF_OUT,1,'min_deposit >= min_value')
        end if
    else
        write(PMF_OUT,115) mtd_item%min_deposit, trim(mtd_item%cv%get_ulabel())
        call mtd_item%cv%conv_to_ivalue(mtd_item%min_deposit)
    end if

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',mtd_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) mtd_item%max_value, trim(mtd_item%cv%get_ulabel())
    call mtd_item%cv%conv_to_ivalue(mtd_item%max_value)

    if( mtd_item%max_value .le. mtd_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_deposit',mtd_item%max_deposit) ) then
        mtd_item%max_deposit = mtd_item%max_value    ! already in internal units
        write(PMF_OUT,125) mtd_item%cv%get_rvalue(mtd_item%max_deposit), trim(mtd_item%cv%get_ulabel())
        if( mtd_item%max_deposit .gt. mtd_item%max_value ) then
            call pmf_utils_exit(PMF_OUT,1,'max_deposit <= max_value')
        end if
    else
        write(PMF_OUT,125) mtd_item%max_deposit, trim(mtd_item%cv%get_ulabel())
        call mtd_item%cv%conv_to_ivalue(mtd_item%max_deposit)
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',mtd_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    end if
    write(PMF_OUT,130) mtd_item%nbins

    return

100 format('    ** Width             : ',F16.6,' [',A,']')
110 format('    ** Min value         : ',F16.6,' [',A,']')
115 format('    ** Min deposit value : ',F16.6,' [',A,']')
120 format('    ** Max value         : ',F16.6,' [',A,']')
125 format('    ** Max deposit value : ',F16.6,' [',A,']')
130 format('    ** Number of bins    : ',I10)

end subroutine mtd_cvs_read_cv

!===============================================================================
! Subroutine:  mtd_cvs_cv_info
!===============================================================================

subroutine mtd_cvs_cv_info(mtd_item)

    use mtd_dat
    use pmf_dat
    use pmf_cvs

    implicit none
    type(CVTypeMTD)    :: mtd_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,145) trim(mtd_item%cv%name)
    write(PMF_OUT,146) trim(mtd_item%cv%ctype)
    write(PMF_OUT,150) mtd_item%cv%get_rvalue(CVContext%CVsValues(mtd_item%cvindx)), &
                    trim(mtd_item%cv%get_ulabel())
    write(PMF_OUT,152) mtd_item%cv%get_rvalue(mtd_item%width), &
                    trim(mtd_item%cv%get_ulabel())
    write(PMF_OUT,155) mtd_item%cv%get_rvalue(mtd_item%min_value), &
                    trim(mtd_item%cv%get_ulabel())
    write(PMF_OUT,156) mtd_item%cv%get_rvalue(mtd_item%min_deposit), &
                    trim(mtd_item%cv%get_ulabel())
    write(PMF_OUT,157) mtd_item%cv%get_rvalue(mtd_item%max_deposit), &
                    trim(mtd_item%cv%get_ulabel())
    write(PMF_OUT,160) mtd_item%cv%get_rvalue(mtd_item%max_value), &
                    trim(mtd_item%cv%get_ulabel())
    write(PMF_OUT,165) mtd_item%nbins

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
152 format('    ** Width             : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
156 format('    ** Min deposit value : ',E16.7,' [',A,']')
157 format('    ** Max deposit value : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I8)

end subroutine mtd_cvs_cv_info

!===============================================================================

end module mtd_cvs_mod
