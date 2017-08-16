!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

    mtd_item%cvindx          = 0   ! CV index
    mtd_item%width           = 0.0 ! width
    mtd_item%min_value       = 0.0 ! left range
    mtd_item%max_value       = 0.0 ! right range
    mtd_item%nbins           = 0.0 ! number of bins
    mtd_item%max_dist        = 0.0 ! maximum distance
    mtd_item%gptheta        = 0.0  ! GP gptheta
    mtd_item%gpperiodicity  = 0.0  ! GP periodicities
    mtd_item%gpfixgrid      = 0    ! GP fixgrid

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

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    type(CVTypeMTD)                    :: mtd_item
    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( mtd_item%cv%pathidx .gt. 0 ) then
        if( PathList(mtd_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in driven mode!')
        end if
    end if

    ! load rest of definition
    if( .not. prmfile_get_real8_by_key(prm_fin,'width',mtd_item%width) ) then
        call pmf_utils_exit(PMF_OUT,1,'width is not specified!')
    else
        write(PMF_OUT,100) mtd_item%width, trim(mtd_item%cv%get_ulabel())
    end if
    call mtd_item%cv%conv_to_ivalue(mtd_item%width)

    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',mtd_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    else
        write(PMF_OUT,110) mtd_item%min_value, trim(mtd_item%cv%get_ulabel())
    end if
    call mtd_item%cv%conv_to_ivalue(mtd_item%min_value)

    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',mtd_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    else
        write(PMF_OUT,120) mtd_item%max_value, trim(mtd_item%cv%get_ulabel())
    end if
    call mtd_item%cv%conv_to_ivalue(mtd_item%max_value)

    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',mtd_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    else
        write(PMF_OUT,125) mtd_item%nbins
    end if

    if( .not. prmfile_get_real8_by_key(prm_fin,'max_dist',mtd_item%max_dist) ) then
        if( mtd_item%cv%is_periodic_cv() ) then
            mtd_item%max_dist = mtd_item%cv%get_period_cv_value()
        else
            mtd_item%max_dist = mtd_item%max_value - mtd_item%min_value
        end if
    else
        write(PMF_OUT,130) mtd_item%max_dist, trim(mtd_item%cv%get_ulabel())
    end if

    if( fmode .eq. 3) then
         if( .not. prmfile_get_real8_by_key(prm_fin,'gptheta',mtd_item%gptheta) ) then
             call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: gptheta is not specified!')
         else
             write(PMF_OUT,140) mtd_item%gptheta, trim(mtd_item%cv%get_ulabel())
             call mtd_item%cv%conv_to_ivalue(mtd_item%gptheta)
         end if
         if( mtd_item%gptheta .le. 0 ) then
             call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: gptheta must be greater than 0!')
         end if

         if( mtd_item%cv%is_periodic_cv() ) then
             mtd_item%gpperiodicity = mtd_item%cv%get_period_cv_value()
         end if

         if(prmfile_get_real8_by_key(prm_fin,'gpperiodicity',mtd_item%gpperiodicity) ) then
             write(PMF_OUT,160) mtd_item%gpperiodicity, trim(mtd_item%cv%get_ulabel())
         else
             if( mtd_item%cv%is_periodic_cv() ) then
                 mtd_item%gpperiodicity = mtd_item%cv%get_period_cv_value()
             end if
             write(PMF_OUT,170) mtd_item%gpperiodicity, trim(mtd_item%cv%get_ulabel())
         end if
         if( mtd_item%gpperiodicity .lt. 0 ) then
             call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: gpperiodicity must be greater or equal than 0!')
         end if

         if( fgpsparsification .eq. 2 ) then
             if( .not. prmfile_get_integer_by_key(prm_fin,'gpnumgrid',mtd_item%gpnumgrid) ) then
                 call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: gpnumgrid is not specified!')
             else
                 write(PMF_OUT,180) mtd_item%gpnumgrid
             end if

             if( prmfile_get_integer_by_key(prm_fin,'gpfixgrid',mtd_item%gpfixgrid) ) then
                 write(PMF_OUT,190) mtd_item%gpfixgrid
             else
                 write(PMF_OUT,200) mtd_item%gpfixgrid
             end if
             if( (mtd_item%gpfixgrid .ne. 0) .and. (mtd_item%gpfixgrid .ne. 1) ) then
                 call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: gpfixgrid iseither 0 or 1!')
             end if
         end if
    end if

    return

100 format('   ** width              : ',F16.6,' [',A,']')
110 format('   ** min value          : ',F16.6,' [',A,']')
120 format('   ** max value          : ',F16.6,' [',A,']')
125 format('   ** number of bins     : ',I10)
130 format('   ** maximum distance   : ',F16.6,' [',A,']')
140 format('   ** GP gptheta value     : ',F16.6,' [',A,']')
160 format('   ** GP gpperiodicity is  : ',F16.6,' [',A,']')
170 format('   ** GP gpperiodicity is  : ',F16.6,' [',A,']','  (default)')
180 format('   ** GP gpnumgrid is  : ',I9)
190 format('   ** GP gpfixgrid is  : ',I9)
200 format('   ** GP gpfixgrid is  : ',I9,'  (default)')

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
    write(PMF_OUT,160) mtd_item%cv%get_rvalue(mtd_item%max_value), &
                    trim(mtd_item%cv%get_ulabel())
    write(PMF_OUT,165) mtd_item%nbins
    write(PMF_OUT,170) mtd_item%max_dist, &
                    trim(mtd_item%cv%get_ulabel())

    if( fmode .eq. 3 ) then
         write(PMF_OUT,180) mtd_item%gptheta, trim(mtd_item%cv%get_ulabel())
         write(PMF_OUT,200) mtd_item%gpperiodicity, trim(mtd_item%cv%get_ulabel())
         if ( fgpsparsification .eq. 2 ) then
              write(PMF_OUT,210) mtd_item%gpnumgrid
              write(PMF_OUT,220) mtd_item%gpfixgrid
         end if
    end if

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
152 format('    ** Width             : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I8)
170 format('    ** Maximum distance  : ',E16.7,' [',A,']')
180 format('    ** GP theta value    : ',F16.6,' [',A,']')
200 format('    ** GP periodicity is : ',F16.6,' [',A,']')
210 format('    ** GP number of grids: ',I8)
220 format('    ** GP fixgrid        : ',I8)


end subroutine mtd_cvs_cv_info

!===============================================================================

end module mtd_cvs_mod
