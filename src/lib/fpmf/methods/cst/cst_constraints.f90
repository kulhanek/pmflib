!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cst_constraints

use pmf_cvs
use cst_dat

implicit none
contains

!===============================================================================
! Subroutine:  cst_constraints_reset_con
!===============================================================================

subroutine cst_constraints_reset_con(cst_item)

    implicit none
    type(CVTypeBM)       :: cst_item
    ! -----------------------------------------------------------------------------

    cst_item%cvindx         = 0       ! CV index
    cst_item%mode           = 'C'
    cst_item%cv             => null()

    cst_item%startvalue     = 0.0d0   ! start value
    cst_item%stopvalue      = 0.0d0   ! stop value
    cst_item%value          = 0.0d0   ! current value in time t

    cst_item%deviation      = 0.0d0   ! deviation between real and value
    cst_item%sdevtot        = 0.0d0   ! total sum of deviation squares
    cst_item%value_set      = .false. ! initial value is user provided

    cst_item%min_value      = 0.0   ! left range
    cst_item%max_value      = 0.0   ! right range
    cst_item%nbins          = 0     ! number of bins

    cst_item%ibin           = 0

end subroutine cst_constraints_reset_con

!===============================================================================
! Subroutine:  cst_constraints_read_con
!===============================================================================

subroutine cst_constraints_read_con(prm_fin,cst_item)

    use prmfile
    use pmf_utils
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeBM)                      :: cst_item
    ! --------------------------------------------
    integer                             :: ibin
    ! --------------------------------------------------------------------------

    if( cst_item%cv%pathidx .gt. 0 ) then
        if( PathList(cst_item%cv%pathidx)%path%driven_mode ) then
            write(PMF_OUT,96)
            return
        end if
    end if

    if( freadranges ) then
        ! ========================
        if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',cst_item%min_value) ) then
            call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
        end if
        write(PMF_OUT,210) cst_item%min_value, trim(cst_item%cv%get_ulabel())
        call cst_item%cv%conv_to_ivalue(cst_item%min_value)

        ! ========================
        if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',cst_item%max_value) ) then
            call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
        end if
        write(PMF_OUT,220) cst_item%max_value, trim(cst_item%cv%get_ulabel())
        call cst_item%cv%conv_to_ivalue(cst_item%max_value)

        if( cst_item%max_value .le. cst_item%min_value ) then
            call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
        end if

        ! ========================
        if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',cst_item%nbins) ) then
            call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
        end if
        if( cst_item%nbins .lt. 1 ) then
            call pmf_utils_exit(PMF_OUT,1,'nbins has to be greater then zero!')
        end if
        write(PMF_OUT,225) cst_item%nbins
    end if

    if( prmfile_get_real8_by_key(prm_fin,'value',cst_item%value) ) then
        write(PMF_OUT,90) cst_item%value, trim(cst_item%cv%get_ulabel())
        call cst_item%cv%conv_to_ivalue(cst_item%value)
        cst_item%value_set = .true.
    else if( prmfile_get_integer_by_key(prm_fin,'value_at_bin',ibin) ) then
        if( (ibin .lt. 1) .or. (ibin .gt. cst_item%nbins) ) then
            call pmf_utils_exit(PMF_OUT,1,'value_at_bin out-of-range!')
        end if
        cst_item%ibin  = ibin
        cst_item%value = (real(ibin,PMFDP)-0.5d0)*(cst_item%max_value -  cst_item%min_value)/real(cst_item%nbins) &
                       + cst_item%min_value
        write(PMF_OUT,91) cst_item%cv%get_rvalue(cst_item%value), trim(cst_item%cv%get_ulabel()), ibin
        cst_item%value_set = .true.
    else
        write(PMF_OUT,95)
        cst_item%value_set = .false.
    end if

    ! modes -----------------------------------------------------------------------
    ! constant mode by default
    cst_item%mode = 'C'

    if( prmfile_get_real8_by_key(prm_fin,'change_to',cst_item%stopvalue) ) then
        cst_item%mode = 'V'
        write(PMF_OUT,100) cst_item%stopvalue, trim(cst_item%cv%get_ulabel())
        call cst_item%cv%conv_to_ivalue(cst_item%stopvalue)
    else if( prmfile_get_integer_by_key(prm_fin,'change_to_bin',ibin) ) then
        cst_item%mode = 'V'
        cst_item%stopvalue = (real(ibin,PMFDP)-0.5d0)*(cst_item%max_value -  cst_item%min_value)/real(cst_item%nbins) &
                           + cst_item%min_value
        write(PMF_OUT,101) cst_item%cv%get_rvalue(cst_item%stopvalue), trim(cst_item%cv%get_ulabel()), ibin
    end if

    if( prmfile_get_real8_by_key(prm_fin,'increment',cst_item%stopvalue) ) then
        if( cst_item%mode .eq. 'V' ) then
            call pmf_utils_exit(PMF_OUT,1,'change_to/change_to_bin and increment keywords cannot be used together!')
        end if
        cst_item%mode = 'I'
        write(PMF_OUT,110) cst_item%stopvalue, trim(cst_item%cv%get_ulabel())
        call cst_item%cv%conv_to_ivalue(cst_item%stopvalue)
    else if( prmfile_get_integer_by_key(prm_fin,'increment_by_bins',ibin) ) then
        if( cst_item%mode .eq. 'V' ) then
            call pmf_utils_exit(PMF_OUT,1,'change_to/change_to_bin and increment_by_bins keywords cannot be used together!')
        end if
        cst_item%stopvalue = real(ibin,PMFDP)*(cst_item%max_value -  cst_item%min_value)/real(cst_item%nbins)
        cst_item%mode = 'I'
        write(PMF_OUT,111) cst_item%cv%get_rvalue(cst_item%stopvalue), trim(cst_item%cv%get_ulabel()), ibin
    end if

    return

210 format('   ** Min value          : ',F16.7,' [',A,']')
220 format('   ** Max value          : ',F16.7,' [',A,']')
225 format('   ** Number of bins     : ',I8)

 90 format('   ** Constrained value  :',E16.7,' [',A,'] (user specified)')
 91 format('   ** Constrained value  :',E16.7,' [',A,'] (user specified at bin: ',I2,')')
 95 format('   ** Constrained value  : value from input coordinates or CST restart file')
 96 format('   ** Constrained value  : controlled by path subsystem')
100 format('   ** Change to value    :',E16.7,' [',A,']')
101 format('   ** Change to value    :',E16.7,' [',A,'] at bin: ', I2)
110 format('   ** Increment value    :',E16.7,' [',A,']')
111 format('   ** Increment value    :',E16.7,' [',A,'] by ', I2, ' bins')

end subroutine cst_constraints_read_con

!===============================================================================
! Subroutine:  cst_constraints_cst_info
!===============================================================================

subroutine cst_constraints_cst_info(cst_item)

    use pmf_paths

    implicit none
    type(CVTypeBM)     :: cst_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,70) trim(cst_item%cv%name)
    write(PMF_OUT,80) trim(cst_item%cv%ctype)

    if( cst_item%mode .eq. 'P' ) then
        write(PMF_OUT,105) trim(PathList(cst_item%cv%pathidx)%path%name)
    end if

    write(PMF_OUT,90)  cst_item%cv%get_rvalue(cst_item%value), &
                       trim(cst_item%cv%get_ulabel())
    write(PMF_OUT,100) cst_item%cv%get_rvalue(CVContext%CVsValues(cst_item%cvindx)), &
                       trim(cst_item%cv%get_ulabel())

    select case(cst_item%mode)
    case ('I','V')
        write(PMF_OUT,110) cst_item%cv%get_rvalue(cst_item%startvalue), &
                           trim(cst_item%cv%get_ulabel())
        write(PMF_OUT,120) cst_item%cv%get_rvalue(cst_item%stopvalue - cst_item%startvalue), &
                           trim(cst_item%cv%get_ulabel())
        write(PMF_OUT,130) cst_item%cv%get_rvalue(cst_item%stopvalue), &
                           trim(cst_item%cv%get_ulabel())
    case ('S')
        write(PMF_OUT,110) cst_item%cv%get_rvalue(cst_item%startvalue), &
                           trim(cst_item%cv%get_ulabel())
    end select

    if( freadranges ) then
        write(PMF_OUT,255) cst_item%cv%get_rvalue(cst_item%min_value), &
                        trim(cst_item%cv%get_ulabel())
        write(PMF_OUT,260) cst_item%cv%get_rvalue(cst_item%max_value), &
                        trim(cst_item%cv%get_ulabel())
        write(PMF_OUT,265) cst_item%nbins
    end if

    return

 70 format('    ** Name              : ',a)
 80 format('    ** Type              : ',a)
 90 format('    ** Constrained value : ',E16.7,' [',A,']')
100 format('    ** Current value     : ',E16.7,' [',A,']')

105 format('    ** Managed by path   : ',A)

110 format('    ** Start value       : ',E16.7,' [',A,']')
120 format('    ** Increment value   : ',E16.7,' [',A,']')
130 format('    ** Final value       : ',E16.7,' [',A,']')

255 format('    ** Min value         : ',E16.7,' [',A,']')
260 format('    ** Max value         : ',E16.7,' [',A,']')
265 format('    ** Number of bins    : ',I9)

end subroutine cst_constraints_cst_info

!===============================================================================
! Subroutine:  cst_constraints_init_all
!===============================================================================

subroutine cst_constraints_init_all

    implicit none
    integer            :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfCONs
        call cst_constraints_cst_init(CONList(i))
    end do

end subroutine cst_constraints_init_all

!===============================================================================
! Subroutine:  cst_constraints_cst_init
!===============================================================================

subroutine cst_constraints_cst_init(cst_item)

    use pmf_paths

    implicit none
    type(CVTypeBM)     :: cst_item
    ! --------------------------------------------------------------------------

    ! get initial value only for constraints specified by user
    if( cst_item%value_set .neqv. .true. ) then
        cst_item%value = CVContext%CVsValues(cst_item%cvindx)
    end if

    if( cst_item%mode .eq. 'S' ) then
        cst_item%value = cst_item%control_values(1)
    end if

    cst_item%startvalue = cst_item%value
    cst_item%deviation  = 0.0
    cst_item%sdevtot    = 0.0

    ! correct increment and stop value -----------------------------------------------
    if( cst_item%mode .eq. 'I' ) then
        cst_item%stopvalue = cst_item%startvalue + cst_item%stopvalue
    end if

    if( cst_item%cv%pathidx .gt. 0 ) then
        if( PathList(cst_item%cv%pathidx)%path%driven_mode ) then
            cst_item%mode = 'P'
        end if
    end if

    return

end subroutine cst_constraints_cst_init

!===============================================================================
! Subroutine:  cst_constraints_increment
!===============================================================================

subroutine cst_constraints_increment

    implicit none
    integer            :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfCONs
        call cst_constraints_cst_increment(CONList(i))
    end do

end subroutine cst_constraints_increment

!===============================================================================
! Subroutine:  cst_constraints_increment_all
!===============================================================================

subroutine cst_constraints_cst_increment(cst_item)

    use pmf_paths

    implicit none
    type(CVTypeBM)     :: cst_item
    ! --------------------------------------------------------------------------

    ! value changes pro-grammatically by managed path
    if( cst_item%mode .eq. 'P' ) then
        cst_item%value = pmf_paths_get_rpos(cst_item%cvindx)
        return
    end if

    if( cst_item%mode .eq. 'C' ) return  ! no change of constraint

    if( (fstep .lt. 0) .or. (fstep .gt. fnstlim) ) return

    if( cst_item%mode .eq. 'S' ) then
        ! set corresponding value
        cst_item%value = cst_item%control_values(fstep)
    else
        ! incrementation is in linear mode
        cst_item%value = cst_item%startvalue + (cst_item%stopvalue - cst_item%startvalue)*fstep/fnstlim
    end if

    return

end subroutine cst_constraints_cst_increment

!===============================================================================
! Subroutine:  cst_constraints_calc_fdxp
!===============================================================================

subroutine cst_constraints_calc_fdxp

    use pmf_timers
    use pmf_cvs

    implicit none
    integer              :: i,ci
    ! --------------------------------------------------------------------------

    ! we need to reset this since gradients are added to it
    CVContextP%CVsValues(:) = 0.0d0
    CVContextP%CVsDrvs(:,:,:) = 0.0d0

    do i=1,NumOfCVs
        CVList(i)%cv%processed = .false.
    end do

    ! timer cannot be in cst_constraints_calc_fdxp_cvitem
    ! which is called recursively
    call pmf_timers_start_timer(PMFLIB_CVS_TIMER)

    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        call cst_constraints_calc_fdxp_cvitem(ci)
        cv(i) = get_deviation(CONList(i)%cv,CONList(i)%value,CVContextP%CVsValues(ci))
    end do

    call pmf_timers_stop_timer(PMFLIB_CVS_TIMER)

end subroutine cst_constraints_calc_fdxp

!===============================================================================
! Subroutine:  cst_constraints_calc_fdxp_cvitem
!===============================================================================

recursive subroutine cst_constraints_calc_fdxp_cvitem(ci)

    implicit none
    integer             :: ci
    ! --------------------------------------------
    integer             :: i
    ! --------------------------------------------------------------------------

    ! already processed CV
    if( CVList(ci)%cv%processed ) return

    ! is it algebraic CV?
    if( .not. CVList(ci)%cv%isalgebraic ) then
        call CVList(ci)%cv%calculate_cv(CrdP,CVContextP)
        CVList(ci)%cv%processed = .true.
    else
        ! first dependent CVs
        if( associated(CVList(ci)%cv%algebraicidxs) ) then
            do i=1,size(CVList(ci)%cv%algebraicidxs)
                call cst_constraints_calc_fdxp_cvitem(CVList(ci)%cv%algebraicidxs(i))
            end do
        end if
        ! then the CV
        call CVList(ci)%cv%calculate_cv(CrdP,CVContextP)
        CVList(ci)%cv%processed = .true.
    end if

end subroutine cst_constraints_calc_fdxp_cvitem

!===============================================================================

!===============================================================================
! Subroutine:  cst_constraints_read_control_file
!===============================================================================

subroutine cst_constraints_read_control_file(cst_item)

    use pmf_utils

    implicit none
    type(CVTypeBM)     :: cst_item
    ! ----------------------------------------------
    integer            :: alloc_failed, ios
    integer            :: i
    ! --------------------------------------------------------------------------

    ! test if control file exists
    if( .not. pmf_utils_fexist(fcstctr) ) then
        write(PMF_OUT,10) trim(fcstctr)
        call pmf_utils_exit(PMF_OUT,1,'fcstctr file does not exist!')
    end if

    ! open control file
    call pmf_utils_open(CST_CTR,fcstctr,'O')

    ! allocate cst_item%control_values
    allocate(cst_item%control_values(fnstlim), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
         call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to allocate memory for control_values!')
    end if

    ! read data from control_file
    ios = 0
    do i=1,fnstlim
       read(CST_CTR,*,iostat=ios) cst_item%control_values(i)
       if ( ios .ne. 0 ) then
            write(PMF_OUT,20) trim(fcstctr), i-1, fnstlim
            call pmf_utils_exit(PMF_OUT,1,'[CST] There is not enough data in fcstctr file!')
       end if
    end do

    close(CST_CTR)

    return

 10 format('[CST] fcstctr file (',A,') does not exist!')
 20 format('[CST] fcstctr file (',A,') has less number of data (',I8,') than fnstlim (',I8,')!')

end subroutine cst_constraints_read_control_file

!===============================================================================

end module cst_constraints
