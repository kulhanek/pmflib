!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module con_constraints

use pmf_cvs

implicit none
contains

!===============================================================================
! Subroutine:  con_constraints_reset_con
!===============================================================================

subroutine con_constraints_reset_con(con_item)

    use con_dat

    implicit none
    type(CVTypeBM)       :: con_item
    ! -----------------------------------------------------------------------------

    con_item%cvindx         = 0       ! CV index
    con_item%mode           = 'C'
    con_item%cv             => null()

    con_item%startvalue     = 0.0d0   ! start value
    con_item%stopvalue      = 0.0d0   ! stop value
    con_item%value          = 0.0d0   ! current value in time t

    con_item%deviation      = 0.0d0   ! deviation between real and value
    con_item%sdevtot        = 0.0d0   ! total sum of deviation squares
    con_item%value_set      = .false. ! initial value is user provided

end subroutine con_constraints_reset_con

!===============================================================================
! Subroutine:  con_constraints_read_con
!===============================================================================

subroutine con_constraints_read_con(prm_fin,con_item)

    use prmfile
    use pmf_utils
    use con_dat
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeBM)                      :: con_item
    ! --------------------------------------------------------------------------

    if( con_item%cv%pathidx .gt. 0 ) then
        if( PathList(con_item%cv%pathidx)%path%driven_mode ) then
            write(PMF_OUT,96)
            return
        end if
    end if

    if( prmfile_get_real8_by_key(prm_fin,'value',con_item%value) ) then
        write(PMF_OUT,90) con_item%value, trim(con_item%cv%get_ulabel())
        call con_item%cv%conv_to_ivalue(con_item%value)
        con_item%value_set = .true.
    else
        write(PMF_OUT,95)
        con_item%value_set = .false.
    end if

    ! modes -----------------------------------------------------------------------
    ! constant mode by default
    con_item%mode = 'C'

    if( prmfile_get_real8_by_key(prm_fin,'change_to',con_item%stopvalue) ) then
        con_item%mode = 'V'
        write(PMF_OUT,100) con_item%stopvalue, trim(con_item%cv%get_ulabel())
        call con_item%cv%conv_to_ivalue(con_item%stopvalue)
    end if

    if( prmfile_get_real8_by_key(prm_fin,'increment',con_item%stopvalue) ) then
        if( con_item%mode .eq. 'V' ) then
            call pmf_utils_exit(PMF_OUT,1,'change_to and increment keywords cannot be used together!')
        end if
        con_item%mode = 'I'
        write(PMF_OUT,110) con_item%stopvalue, trim(con_item%cv%get_ulabel())
        call con_item%cv%conv_to_ivalue(con_item%stopvalue)
    end if

    return

 90 format('   ** Constrained value  :',E16.7,' [',A,'] (user specified)')
 95 format('   ** Constrained value  : value from input coordinates or CON restart file')
 96 format('   ** Constrained value  : controlled by path subsystem')
100 format('   ** Change to value    :',E16.7' [',A,']')
110 format('   ** Increment value    :',E16.7' [',A,']')

end subroutine con_constraints_read_con

!===============================================================================
! Subroutine:  con_constraints_con_info
!===============================================================================

subroutine con_constraints_con_info(con_item)

    use con_dat
    use pmf_paths

    implicit none
    type(CVTypeBM)     :: con_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,70) trim(con_item%cv%name)
    write(PMF_OUT,80) trim(con_item%cv%ctype)

    if( con_item%mode .eq. 'P' ) then
        write(PMF_OUT,105) trim(PathList(con_item%cv%pathidx)%path%name)
    end if

    write(PMF_OUT,90)  con_item%cv%get_rvalue(con_item%value), &
                       trim(con_item%cv%get_ulabel())
    write(PMF_OUT,100) con_item%cv%get_rvalue(CVContext%CVsValues(con_item%cvindx)), &
                       trim(con_item%cv%get_ulabel())

    select case(con_item%mode)
    case ('I','V')
        write(PMF_OUT,110) con_item%cv%get_rvalue(con_item%startvalue), &
                           trim(con_item%cv%get_ulabel())
        write(PMF_OUT,120) con_item%cv%get_rvalue(con_item%stopvalue - con_item%startvalue), &
                           trim(con_item%cv%get_ulabel())
        write(PMF_OUT,130) con_item%cv%get_rvalue(con_item%stopvalue), &
                           trim(con_item%cv%get_ulabel())
    end select

    return

 70 format('    ** Name              : ',a)
 80 format('    ** Type              : ',a)
 90 format('    ** Constrained value : ',E16.7,' [',A,']')
100 format('    ** Current value     : ',E16.7,' [',A,']')

105 format('    ** Managed by path   : ',A)

110 format('    ** Start value       : ',E16.7,' [',A,']')
120 format('    ** Increment value   : ',E16.7,' [',A,']')
130 format('    ** Final value       : ',E16.7,' [',A,']')

end subroutine con_constraints_con_info

!===============================================================================
! Subroutine:  con_constraints_init_all
!===============================================================================

subroutine con_constraints_init_all

    use con_dat

    implicit none
    integer            :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfCONs
        call con_constraints_con_init(CONList(i))
    end do

end subroutine con_constraints_init_all

!===============================================================================
! Subroutine:  con_constraints_con_init
!===============================================================================

subroutine con_constraints_con_init(con_item)

    use pmf_dat
    use con_dat
    use pmf_paths

    implicit none
    type(CVTypeBM)     :: con_item
    ! --------------------------------------------------------------------------

    ! get initial value only for constraints specified by user
    if( con_item%value_set .neqv. .true. ) then
        con_item%value = CVContext%CVsValues(con_item%cvindx)
    end if

    con_item%startvalue = con_item%value
    con_item%deviation  = 0.0
    con_item%sdevtot    = 0.0

    ! correct increment and stop value -----------------------------------------------
    if( con_item%mode .eq. 'I' ) then
        con_item%stopvalue = con_item%startvalue + con_item%stopvalue
    end if

    if( con_item%cv%pathidx .gt. 0 ) then
        if( PathList(con_item%cv%pathidx)%path%driven_mode ) then
            con_item%mode = 'P'
        end if
    end if

    return

end subroutine con_constraints_con_init

!===============================================================================
! Subroutine:  con_constraints_increment
!===============================================================================

subroutine con_constraints_increment

    use con_dat

    implicit none
    integer            :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfCONs
        call con_constraints_con_increment(CONList(i))
    end do

end subroutine con_constraints_increment

!===============================================================================
! Subroutine:  con_constraints_increment_all
!===============================================================================

subroutine con_constraints_con_increment(con_item)

    use con_dat
    use pmf_dat
    use pmf_paths

    implicit none
    type(CVTypeBM)     :: con_item
    ! --------------------------------------------------------------------------

    ! value changes programmatically by managed path
    if( con_item%mode .eq. 'P' ) then
        con_item%value = pmf_paths_get_rpos(con_item%cvindx)
        return
    end if

    if( con_item%mode .eq. 'C' ) return  ! no change of constraint

    if( (fstep .lt. 0) .or. (fstep .gt. fnstlim) ) return

    ! incrementation is in linear mode
    con_item%value = con_item%startvalue + (con_item%stopvalue - con_item%startvalue)*fstep/fnstlim

    ! image incremented value
    con_item%value = get_imaged_value(con_item%cv,con_item%value)

    return

end subroutine con_constraints_con_increment

!===============================================================================
! Subroutine:  con_constraints_calc_fdxp
!===============================================================================

subroutine con_constraints_calc_fdxp

    use pmf_dat
    use con_dat
    use pmf_timers

    implicit none
    integer              :: i,ci
    ! --------------------------------------------------------------------------

    ! we need to reset this since gradients are added to it
    CVContextP%CVsValues(:) = 0.0d0
    CVContextP%CVsDrvs(:,:,:) = 0.0d0

    do i=1,NumOfCVs
        CVList(i)%cv%processed = .false.
    end do

    ! timer cannot be in con_constraints_calc_fdxp_cvitem
    ! which is called recursively
    call pmf_timers_start_timer(PMFLIB_CVS_TIMER)

    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        call con_constraints_calc_fdxp_cvitem(ci)
        cv(i) = get_deviation(CONList(i)%cv,CONList(i)%value,CVContextP%CVsValues(ci))
    end do

    call pmf_timers_stop_timer(PMFLIB_CVS_TIMER)

end subroutine con_constraints_calc_fdxp

!===============================================================================
! Subroutine:  con_constraints_calc_fdxp
!===============================================================================

recursive subroutine con_constraints_calc_fdxp_cvitem(ci)

    use pmf_dat

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
                call con_constraints_calc_fdxp_cvitem(CVList(ci)%cv%algebraicidxs(i))
            end do
        end if
        ! then the CV
        call CVList(ci)%cv%calculate_cv(CrdP,CVContextP)
        CVList(ci)%cv%processed = .true.
    end if

end subroutine con_constraints_calc_fdxp_cvitem

!===============================================================================

end module con_constraints
