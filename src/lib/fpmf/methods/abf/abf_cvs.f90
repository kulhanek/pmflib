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

    abf_item%cvindx         = 0   ! CV index
    abf_item%cv             => null()
    abf_item%set            = 0
    abf_item%min_value      = 0.0 ! left range
    abf_item%max_value      = 0.0 ! right range
    abf_item%nbins          = 0.0 ! number of bins
    abf_item%maxforce       = 0.0 ! max force to be applied
    abf_item%switch         = 0.0 ! switch to zero width
    abf_item%fgplen         = 0.0  ! characteristic lengh-scale
    abf_item%fgpsigmaoffset = 0.0  ! force sigma offset
    abf_item%fgpsigmafac    = 1.0  ! force sigma factor

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

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeABF)                     :: abf_item
    ! -----------------------------------------------
    type(UnitType)                      :: forceunit
    type(UnitType)                      :: forceunit_recip
    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( abf_item%cv%pathidx .gt. 0 ) then
        if( PathList(abf_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in driven mode!')
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
    write(PMF_OUT,125) abf_item%nbins

    ! ============================================
    ! prepare force unit
    forceunit = pmf_unit_div_units( EnergyUnit, abf_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    ! ========================
    ! optional parameter
    if( prmfile_get_real8_by_key(prm_fin,'maxforce',abf_item%maxforce) ) then
        abf_item%maxforce = pmf_unit_get_ivalue(forceunit,abf_item%maxforce)
    end if

    if( abf_item%maxforce .eq. 0.d0 ) then
        write(PMF_OUT,135)
    else
        write(PMF_OUT,130) pmf_unit_get_rvalue(forceunit,abf_item%maxforce),trim(pmf_unit_label(forceunit))
    end if

    if( abf_item%maxforce .lt. 0.0d0 ) then
        call pmf_utils_exit(PMF_OUT,1,'maxforce must be zero or greater than zero!')
    end if

    ! ========================
    ! optional parameter
    if( prmfile_get_real8_by_key(prm_fin,'switch',abf_item%switch) ) then
        call abf_item%cv%conv_to_ivalue(abf_item%switch)
    end if

    if( abf_item%switch .eq. 0.d0 ) then
        write(PMF_OUT,145)
    else
        write(PMF_OUT,140) abf_item%cv%get_rvalue(abf_item%switch),trim(abf_item%cv%get_ulabel())
    end if

    if( abf_item%switch .gt. 0.5d0*(abf_item%max_value - abf_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'switch must be lower than half of CV interval!')
    end if

    ! ============================================
    if( feimode .eq. 3 ) then

        ! ========================
        if( prmfile_get_real8_by_key(prm_fin,'gplen',abf_item%fgplen) ) then
            abf_item%maxforce = pmf_unit_get_ivalue(forceunit,abf_item%maxforce)
        else
            abf_item%fgplen = 4*(abf_item%max_value - abf_item%min_value) &
                              / real(abf_item%nbins)
        end if
        write(PMF_OUT,150) abf_item%cv%get_rvalue(abf_item%fgplen),trim(abf_item%cv%get_ulabel())


        ! ========================
        if( prmfile_get_real8_by_key(prm_fin,'soffset',abf_item%fgpsigmaoffset) ) then
            abf_item%fgpsigmaoffset = pmf_unit_get_ivalue(forceunit,abf_item%fgpsigmaoffset)
        end if
        write(PMF_OUT,160) pmf_unit_get_rvalue(forceunit,abf_item%fgpsigmaoffset), &
                           trim(pmf_unit_label(forceunit))

        ! ========================
        if( prmfile_get_real8_by_key(prm_fin,'sfac',abf_item%fgpsigmafac) ) then
            abf_item%fgpsigmaoffset = pmf_unit_get_ivalue(forceunit_recip,abf_item%fgpsigmafac)
        end if
        write(PMF_OUT,170) pmf_unit_get_rvalue(forceunit_recip,abf_item%fgpsigmafac), &
                           trim(pmf_unit_label(forceunit_recip))

    end if

    return

110 format('   ** Min value          : ',F16.7,' [',A,']')
120 format('   ** Max value          : ',F16.7,' [',A,']')
125 format('   ** Number of bins     : ',I8)
130 format('   ** Max allowed force  : ',F16.7,' [',A,']')
135 format('   ** Max allowed force  : -disabled-')
140 format('   ** Switch width       : ',F16.7,' [',A,']')
145 format('   ** Switch width       : -disabled-')
150 format('   ** Character length   : ',F16.7,' [',A,']')
160 format('   ** Sigma offset       : ',F16.7,' [',A,']')
170 format('   ** Sigma factor       : ',F16.7,' [',A,']')

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
    ! -----------------------------------------------
    type(UnitType)  :: forceunit
    type(UnitType)  :: forceunit_recip
    ! --------------------------------------------------------------------------

    ! prepare force unit
    forceunit       = pmf_unit_div_units( EnergyUnit, abf_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    write(PMF_OUT,145) trim(abf_item%cv%name)
    write(PMF_OUT,146) trim(abf_item%cv%ctype)
    write(PMF_OUT,150) abf_item%cv%get_rvalue(CVContext%CVsValues(abf_item%cvindx)), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,155) abf_item%cv%get_rvalue(abf_item%min_value), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,160) abf_item%cv%get_rvalue(abf_item%max_value), &
                    trim(abf_item%cv%get_ulabel())
    write(PMF_OUT,165) abf_item%nbins

    if( abf_item%maxforce .eq. 0.d0 ) then
        write(PMF_OUT,135)
    else
        write(PMF_OUT,130) pmf_unit_get_rvalue(forceunit,abf_item%maxforce), &
                           trim(pmf_unit_label(forceunit))
    end if

    if( abf_item%switch .eq. 0.d0 ) then
        write(PMF_OUT,142)
    else
        write(PMF_OUT,140) abf_item%cv%get_rvalue(abf_item%switch), &
                           trim(abf_item%cv%get_ulabel())
    end if

    ! gaussian process setup
    if( feimode .eq. 3 ) then
        write(PMF_OUT,170) abf_item%cv%get_rvalue(abf_item%fgplen), &
                           trim(abf_item%cv%get_ulabel())
        write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,abf_item%fgpsigmaoffset), &
                           trim(pmf_unit_label(forceunit))
        write(PMF_OUT,190) pmf_unit_get_rvalue(forceunit_recip,abf_item%fgpsigmafac), &
                           trim(pmf_unit_label(forceunit_recip))
    end if

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)
130 format('    ** Max allowed force : ',E16.7,' [',A,']')
135 format('    ** Max allowed force : -disabled-')
140 format('    ** Switch width      : ',E16.7,' [',A,']')
142 format('    ** Switch width      : -disabled-')
170 format('    ** Character length  : ',E16.7,' [',A,']')
180 format('    ** Sigma offset      : ',E16.7,' [',A,']')
190 format('    ** Sigma factor      : ',E16.7,' [',A,']')

end subroutine abf_cvs_cv_info

!===============================================================================

end module abf_cvs
