!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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

    abf_item%cvindx         = 0         ! CV index
    abf_item%cv             => null()
    abf_item%set            = 0
    abf_item%min_value      = 0.0       ! left range
    abf_item%max_value      = 0.0       ! right range
    abf_item%nbins          = 0         ! number of bins
    abf_item%wfac           = 1.0d0
    abf_item%buffer         = 0.0d0
    abf_item%energy         = 0.0d0
    abf_item%shake          = .false.

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
    ! -----------------------------------------------
    type(UnitType)                      :: forceunit
    character(1)                        :: buffer
    ! --------------------------------------------------------------------------

! used CV cannot be controlled by the path subsystem
    if( abf_item%cv%pathidx .gt. 0 ) then
        if( PathList(abf_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in a driven mode!')
        end if
    end if

    if( prmfile_get_logical_by_key(prm_fin,'shake',abf_item%shake) ) then
        write(PMF_OUT,300) prmfile_onoff(abf_item%shake)
    else
        write(PMF_OUT,300) prmfile_onoff(abf_item%shake)
    end if

    if( abf_item%shake ) return

! main CV setup
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

    ! ========================
    abf_item%buffer = 0.0
    if( prmfile_get_real8_by_key(prm_fin,'buffer',abf_item%buffer) ) then
        write(PMF_OUT,180) abf_item%buffer, trim(abf_item%cv%get_ulabel())
        call abf_item%cv%conv_to_ivalue(abf_item%buffer)

        if( abf_item%buffer .lt. 0.0d0 ) then
            call pmf_utils_exit(PMF_OUT,1,'buffer has to be greater or equal to zero!')
        end if
    end if

    ! ========================
    if( feimode .eq. 2 ) then
        if( .not. prmfile_get_real8_by_key(prm_fin,'wfac',abf_item%wfac) ) then
            call pmf_utils_exit(PMF_OUT,1,'wfac is not specified!')
        end if
        call pmf_ctrl_check_real8('ABF','wfac',abf_item%wfac,0.0d0,CND_GT,'f12.2')
        write(PMF_OUT,170) abf_item%wfac
    end if

    if( prmfile_get_string_by_key(prm_fin,'value',buffer) ) then
        if( trim(buffer) .eq. '@' ) then
            abf_item%set_value = .true.
            write(PMF_OUT,201) '  @initial value'
        end if
    end if

    ! ========================
    if( fusmode ) then
    ! prepare force unit
    forceunit       = pmf_unit_div_units( EnergyUnit, pmf_unit_power_unit(abf_item%cv%unit,2) )

    if( abf_item%set_value .eqv. .false. ) then
        if( prmfile_get_real8_by_key(prm_fin,'value',abf_item%target_value) ) then
            call abf_item%cv%conv_to_ivalue(abf_item%target_value)
            write(PMF_OUT,200) abf_item%cv%get_rvalue(abf_item%target_value), trim(abf_item%cv%get_ulabel())
        else
            call pmf_utils_exit(PMF_OUT,1,'Target value is not specified (required by US-ABF)!')
        end if
    end if

    ! force constant =========
    if( .not. prmfile_get_real8_by_key(prm_fin,'force_constant',abf_item%force_constant) ) then
        call pmf_utils_exit(PMF_OUT,1,'force_constant is not specified!')
    end if
    call pmf_unit_conv_to_ivalue(forceunit,abf_item%force_constant)
    write(PMF_OUT,210) pmf_unit_get_rvalue(forceunit,abf_item%force_constant), trim(pmf_unit_label(forceunit))
    end if



    return

110 format('    ** Min value         : ',F16.7,' [',A,']')
120 format('    ** Max value         : ',F16.7,' [',A,']')
125 format('    ** Number of bins    : ',I8)
170 format('    ** KS W-factor       : ',F16.7)
180 format('    ** Buffer width      : ',E16.7,' [',A,']')

200 format('    ** Value             : ',F16.7,' [',A,']')
201 format('    ** Value             : ',A)
210 format('    ** Force constant    : ',F16.7,' [',A,']')

300 format('    ** SHAKE             : ',A)

end subroutine abf_cvs_read_cv

!===============================================================================
! Subroutine:  abf_cvs_init_values
!===============================================================================

subroutine abf_cvs_init_values()

    use abf_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit

    implicit none
    integer             :: i,gi0
    ! --------------------------------------------------------------------------

    do i=1,NumOfABFCVs
        if( ABFCVList(i)%set_value ) then
            ABFCVList(i)%target_value = CVContext%CVsValues(ABFCVList(i)%cvindx)
        end if
    end do

    ! align bias if requested
    ! FIXME
!    if( falignbias ) then
!        do i=1,NumOfABFCVs
!            cvhist(i,1) = ABFCVList(i)%target_value
!        end do
!        gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvhist(:,1))
!        if( gi0 .ne. 0 ) then
!            do i=1,NumOfABFCVs
!                ABFCVList(i)%target_value = abfaccu%binpos(i,gi0)
!            end do
!        end if
!    end if

end subroutine abf_cvs_init_values

!===============================================================================
! Subroutine:  abf_cvs_cv_info
!===============================================================================

subroutine abf_cvs_cv_info(abf_item)

    use abf_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit
    use prmfile

    implicit none
    type(CVTypeABF) :: abf_item
    ! -----------------------------------------------
    type(UnitType)                      :: forceunit
    ! --------------------------------------------------------------------------

    write(PMF_OUT,145) trim(abf_item%cv%name)
    write(PMF_OUT,146) trim(abf_item%cv%ctype)
    write(PMF_OUT,150) abf_item%cv%get_rvalue(CVContext%CVsValues(abf_item%cvindx)), &
                    trim(abf_item%cv%get_ulabel())

    write(PMF_OUT,300) prmfile_onoff(abf_item%shake)
    if( abf_item%shake ) return

    write(PMF_OUT,155) abf_item%cv%get_rvalue(abf_item%min_value), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,160) abf_item%cv%get_rvalue(abf_item%max_value), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,165) abf_item%nbins
    write(PMF_OUT,180) abf_item%cv%get_rvalue(abf_item%buffer), &
                    trim(abf_item%cv%get_ulabel())

    ! ========================
    if( feimode .eq. 2 ) then
    write(PMF_OUT,170) abf_item%wfac
    end if

    ! ========================
    if( fusmode ) then
    ! prepare force unit
    forceunit       = pmf_unit_div_units( EnergyUnit, pmf_unit_power_unit(abf_item%cv%unit,2) )
    write(PMF_OUT,210) abf_item%cv%get_rvalue(abf_item%target_value), &
                        trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,220) pmf_unit_get_rvalue(forceunit,abf_item%force_constant), &
                        trim(pmf_unit_label(forceunit))
    end if

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)
170 format('    ** KS W-factor       : ',F16.7)
180 format('    ** Buffer width      : ',E16.7,' [',A,']')

210 format('    ** Target value      : ',E16.7,' [',A,']')
220 format('    ** Force constant    : ',E16.7,' [',A,']')

300 format('    ** SHAKE             : ',A)

end subroutine abf_cvs_cv_info

!===============================================================================

end module abf_cvs
