!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

module rst_restraints

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  rst_restraints_reset_rst
!===============================================================================

subroutine rst_restraints_reset_rst(rst_item)

    use rst_dat

    implicit none
    type(CVTypeUM)       :: rst_item
    ! --------------------------------------------------------------------------

    rst_item%cvindx          = 0    ! CV index
    rst_item%mode            = ' '  ! mode

    rst_item%startvalue      = 0.0 ! start value
    rst_item%stopvalue       = 0.0 ! stop value
    rst_item%target_value    = 0.0 ! required value of restraint

    rst_item%left_value      = 0.0 ! left value for wall restraint
    rst_item%right_value     = 0.0 ! right value for wall restraint

    rst_item%force_constant  = 0.0 ! sigma value
    rst_item%deviation       = 0.0 ! deviation between real and actual value
    rst_item%energy          = 0.0 ! restraint energy

    rst_item%set_value       = .false.

end subroutine rst_restraints_reset_rst

!===============================================================================
! Subroutine:  rst_restraints_read_rst
!===============================================================================

subroutine rst_restraints_read_rst(prm_fin,rst_item)

    use prmfile
    use rst_dat
    use pmf_cvs
    use pmf_dat
    use pmf_unit
    use pmf_utils
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    type(CVTypeUM)                     :: rst_item
    ! -----------------------------------------------
    logical                            :: value_set
    logical                            :: left_value_set
    logical                            :: right_value_set
    logical                            :: ctr_set
    character(1)                       :: buffer
    type(UnitType)                     :: forceunit
    ! --------------------------------------------------------------------------

    ! prepare force unit
    forceunit = pmf_unit_div_units( EnergyUnit, pmf_unit_power_unit(rst_item%cv%unit,2) )

    if( rst_item%cv%pathidx .gt. 0 ) then
        if( PathList(rst_item%cv%pathidx)%path%driven_mode ) then        
            write(PMF_OUT,101) '  controlled by path subsystem'
            rst_item%mode = 'P'
            ! force_constant ==============================================================
            if( .not. prmfile_get_real8_by_key(prm_fin,'force_constant',rst_item%force_constant) ) then
                call pmf_utils_exit(PMF_OUT,1,'force_constant is not specified!')
            end if
            call pmf_unit_conv_to_ivalue(forceunit,rst_item%force_constant)
            write(PMF_OUT,110) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), trim(pmf_unit_label(forceunit))

            ! read histogram setup if necessary
            call rst_restraints_read_rst_hist(prm_fin,rst_item)

            return
        end if
    end if

    value_set = .false.
    left_value_set = .false.
    right_value_set = .false.
    ctr_set = .false.

    if( prmfile_get_real8_by_key(prm_fin,'value',rst_item%startvalue) ) then
        value_set = .true.
        call rst_item%cv%conv_to_ivalue(rst_item%startvalue)
        write(PMF_OUT,100) rst_item%cv%get_rvalue(rst_item%startvalue), trim(rst_item%cv%get_ulabel())
    end if

    if( prmfile_get_string_by_key(prm_fin,'value',buffer) ) then
        if( trim(buffer) .eq. '@' ) then
            value_set = .true.
            rst_item%set_value = .true.
            write(PMF_OUT,101) '  @initial value'
        end if
    end if

    if( prmfile_get_real8_by_key(prm_fin,'left_value',rst_item%left_value) ) then
        left_value_set = .true.
        call rst_item%cv%conv_to_ivalue(rst_item%left_value)
        write(PMF_OUT,102) rst_item%cv%get_rvalue(rst_item%left_value), trim(rst_item%cv%get_ulabel())
    end if

    if( prmfile_get_real8_by_key(prm_fin,'right_value',rst_item%right_value) ) then
        right_value_set = .true.
        call rst_item%cv%conv_to_ivalue(rst_item%right_value)
        write(PMF_OUT,104) rst_item%cv%get_rvalue(rst_item%right_value), trim(rst_item%cv%get_ulabel())
    end if

    if( prmfile_get_string_by_key(prm_fin,'control_file',frstctr) ) then
        ctr_set = .true.
        write(PMF_OUT,105) trim(frstctr)
        call rst_restraints_read_control_file(rst_item)
    end if

    if( .not. (value_set .or. left_value_set .or. right_value_set .or. ctr_set) ) then
        call pmf_utils_exit(PMF_OUT,1,'Neither value nor left_value and right_value nor control_file keys are specified!')
    end if

    if( value_set .and. (left_value_set .or. right_value_set) ) then
        call pmf_utils_exit(PMF_OUT,1,'Either value or left_value and right_value keys can be specified!')
    end if

    if( value_set .and. ctr_set ) then
        call pmf_utils_exit(PMF_OUT,1, &
        '>>> ERROR: Either value or left_value and right_value or control_file keys can be specified!')
    end if

    if( ctr_set .and. (left_value_set .or. right_value_set) ) then
        call pmf_utils_exit(PMF_OUT,1,&
        '>>> ERROR: Either value or left_value and right_value or control_file keys can be specified!')
    end if

    if( (.not. value_set) .and. (.not. ctr_set) .and. (left_value_set .or. right_value_set) ) then
        if( .not. (left_value_set .and. right_value_set) ) then
            call pmf_utils_exit(PMF_OUT,1,'Either left_value or right_value key is missing!')
        else
            ! wall restraint mode
            rst_item%mode = 'W'
        end if
    else if( value_set ) then
        ! constant mode
        rst_item%mode = 'C'
    else if( ctr_set ) then
        ! steering mode
        rst_item%mode = 'S'
    end if

    ! force_constant ==============================================================
    if( .not. prmfile_get_real8_by_key(prm_fin,'force_constant',rst_item%force_constant) ) then
        call pmf_utils_exit(PMF_OUT,1,'force_constant is not specified!')
    end if
    call pmf_unit_conv_to_ivalue(forceunit,rst_item%force_constant)
    write(PMF_OUT,110) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), trim(pmf_unit_label(forceunit))

    ! modes =======================================================================

    if( rst_item%mode .ne. 'W' ) then
        if( prmfile_get_real8_by_key(prm_fin,'change_to',rst_item%stopvalue) ) then
            if( rst_item%mode .eq. 'S' ) then
                call pmf_utils_exit(PMF_OUT,1,'change_to and control_file keywords cannot be used together!')
            end if
            rst_item%mode = 'V'
            call rst_item%cv%conv_to_ivalue(rst_item%stopvalue)
            write(PMF_OUT,120) rst_item%cv%get_rvalue(rst_item%stopvalue), trim(rst_item%cv%get_ulabel())
        end if

        if( prmfile_get_real8_by_key(prm_fin,'increment',rst_item%stopvalue) ) then
            if( rst_item%mode .eq. 'V' ) then
                call pmf_utils_exit(PMF_OUT,1,'change_to and increment keywords cannot be used together!')
            end if
            if( rst_item%mode .eq. 'S' ) then
                call pmf_utils_exit(PMF_OUT,1,'increment and control_file keywords cannot be used together!')
            end if
            rst_item%mode = 'I'
            call rst_item%cv%conv_to_ivalue(rst_item%stopvalue)
            write(PMF_OUT,130) rst_item%cv%get_rvalue(rst_item%stopvalue), trim(rst_item%cv%get_ulabel())
        end if
    else if ( rst_item%mode .ne. 'S' ) then
        if( prmfile_get_real8_by_key(prm_fin,'change_to',rst_item%stopvalue) ) then
            call pmf_utils_exit(PMF_OUT,1,'change_to key cannot be specified in wall restraint mode!')
        end if

        if( prmfile_get_real8_by_key(prm_fin,'increment',rst_item%stopvalue) ) then
            call pmf_utils_exit(PMF_OUT,1,'increment key cannot be specified in wall restraint mode!')
        end if
    end if

    rst_item%target_value = rst_item%startvalue

    ! read histogram setup if necessary
    call rst_restraints_read_rst_hist(prm_fin,rst_item)

    return

100 format('   ** Value          :',F16.6,' [',A,']')
101 format('   ** Value          :',A)
102 format('   ** Left value <=  :',F16.6,' [',A,']')
104 format('   ** Right value >= :',F16.6,' [',A,']')
105 format('   ** Control file   :',A)
110 format('   ** Force constant :',F16.6,' [',A,']')
120 format('   ** Change to      :',F16.6,' [',A,']')
130 format('   ** Increment by   :',F16.6,' [',A,']')

end subroutine rst_restraints_read_rst

!===============================================================================
! Subroutine:  rst_restraints_read_rst_hist
!===============================================================================

subroutine rst_restraints_read_rst_hist(prm_fin,rst_item)

    use prmfile
    use rst_dat
    use pmf_cvs
    use pmf_dat
    use pmf_unit
    use pmf_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    type(CVTypeUM)                     :: rst_item
    ! --------------------------------------------------------------------------

    ! only if we would like to update restart file and of if restart is explicitly required
    if( (fhistupdate .eq. 0) .and. (frestart .eqv. .false.) ) return

    ! histogram ----------------------------------------------------------------
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',rst_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified but histogram accumulation is required!')
    end if
    write(PMF_OUT,210) rst_item%min_value, trim(rst_item%cv%get_ulabel())
    call rst_item%cv%conv_to_ivalue(rst_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',rst_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified but histogram accumulation is required!')
    end if
    write(PMF_OUT,220) rst_item%max_value, trim(rst_item%cv%get_ulabel())
    call rst_item%cv%conv_to_ivalue(rst_item%max_value)

    if( rst_item%max_value .le. rst_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',rst_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified but histogram accumulation is required!')
    end if
    write(PMF_OUT,225) rst_item%nbins

210 format('   ** Min value      : ',F16.6,' [',A,']')
220 format('   ** Max value      : ',F16.6,' [',A,']')
225 format('   ** Number of bins : ',I8)

end subroutine rst_restraints_read_rst_hist

!===============================================================================
! Subroutine:  rst_restraints_rst_info
!===============================================================================

subroutine rst_restraints_rst_info(rst_item)

    use rst_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit
    use pmf_paths

    implicit none
    type(CVTypeUM)     :: rst_item
    ! -----------------------------------------------
    real(PMFDP)        :: increment
    type(UnitType)     :: forceunit
    ! -----------------------------------------------------------------------------

    ! prepare force unit
    forceunit = pmf_unit_div_units( EnergyUnit, pmf_unit_power_unit(rst_item%cv%unit,2) )

    write(PMF_OUT,100) trim(rst_item%cv%name)
    write(PMF_OUT,105) trim(rst_item%cv%ctype)

    if( rst_item%cv%pathidx .gt. 0 ) then
        if( PathList(rst_item%cv%pathidx)%path%driven_mode ) then
            rst_item%set_value = .false.
            rst_item%mode = 'P'
            rst_item%target_value = pmf_paths_get_rpos(rst_item%cvindx)
            write(PMF_OUT,115) trim(PathList(rst_item%cv%pathidx)%path%name)
        end if
    end if

    if( rst_item%set_value ) then
        rst_item%startvalue = CVContext%CVsValues(rst_item%cvindx)
        rst_item%target_value = rst_item%startvalue
    end if

    rst_item%deviation = rst_item%cv%get_deviation(CVContext%CVsValues(rst_item%cvindx),rst_item%target_value)

    select case(rst_item%mode)
        case('P')
            write(PMF_OUT,120) 'P','managed by path'
            write(PMF_OUT,130) rst_item%cv%get_rvalue(rst_item%target_value), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,110) rst_item%cv%get_rvalue(CVContext%CVsValues(rst_item%cvindx)), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,170) rst_item%cv%get_rvalue(rst_item%deviation), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), &
                                trim(pmf_unit_label(forceunit))
        case('C')
            write(PMF_OUT,120) 'C','constant value mode'
            write(PMF_OUT,130) rst_item%cv%get_rvalue(rst_item%target_value), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,110) rst_item%cv%get_rvalue(CVContext%CVsValues(rst_item%cvindx)), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,170) rst_item%cv%get_rvalue(rst_item%deviation), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), &
                                trim(pmf_unit_label(forceunit))
        case('W')
            write(PMF_OUT,120) 'W','wall restraint mode'
            write(PMF_OUT,140) rst_item%cv%get_rvalue(rst_item%left_value), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,110) rst_item%cv%get_rvalue(CVContext%CVsValues(rst_item%cvindx)), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,150) rst_item%cv%get_rvalue(rst_item%right_value), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), &
                                trim(pmf_unit_label(forceunit))
        case('I')
            write(PMF_OUT,120) 'I','incremental mode'
            write(PMF_OUT,160) rst_item%cv%get_rvalue(rst_item%startvalue), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,110) rst_item%cv%get_rvalue(CVContext%CVsValues(rst_item%cvindx)), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,170) rst_item%cv%get_rvalue(rst_item%deviation), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), &
                                trim(pmf_unit_label(forceunit))
            rst_item%stopvalue = rst_item%stopvalue + rst_item%startvalue
            increment = rst_item%cv%get_deviation(rst_item%stopvalue,rst_item%startvalue)
            write(PMF_OUT,190) rst_item%cv%get_rvalue(increment), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,200) rst_item%cv%get_rvalue(rst_item%stopvalue), &
                                trim(rst_item%cv%get_ulabel())
        case('V')
            write(PMF_OUT,120) 'V','change to value'
            write(PMF_OUT,160) rst_item%cv%get_rvalue(rst_item%startvalue), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,110) rst_item%cv%get_rvalue(CVContext%CVsValues(rst_item%cvindx)), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,170) rst_item%cv%get_rvalue(rst_item%deviation), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), &
                                trim(pmf_unit_label(forceunit))
            increment = rst_item%cv%get_deviation(rst_item%stopvalue,rst_item%startvalue)
            write(PMF_OUT,190) rst_item%cv%get_rvalue(increment), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,200) rst_item%cv%get_rvalue(rst_item%stopvalue), &
                                trim(rst_item%cv%get_ulabel())
        case('S')
            write(PMF_OUT,120) 'S','control mode'
            write(PMF_OUT,160) rst_item%cv%get_rvalue(rst_item%startvalue), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,110) rst_item%cv%get_rvalue(CVContext%CVsValues(rst_item%cvindx)), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,170) rst_item%cv%get_rvalue(rst_item%deviation), &
                                trim(rst_item%cv%get_ulabel())
            write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,rst_item%force_constant), &
                                trim(pmf_unit_label(forceunit))
    end select

    if( (fhistupdate .gt. 0) .or. (frestart .eqv. .true.) ) then
        write(PMF_OUT,355) rst_item%cv%get_rvalue(rst_item%min_value), &
                        trim(rst_item%cv%get_ulabel())
        write(PMF_OUT,360) rst_item%cv%get_rvalue(rst_item%max_value), &
                        trim(rst_item%cv%get_ulabel())
        write(PMF_OUT,365) rst_item%nbins
    end if

    return

100 format('    ** Name             : ',A)
105 format('    ** Type             : ',A)
110 format('    ** Current value    : ',E16.9,' [',A,']')

115 format('    ** Managed by path   : ',A)

120 format('    ** Restraining mode : ',A2,' (',A,')')
130 format('    ** Target value     : ',E16.9,' [',A,']')
140 format('    ** Left value <=    : ',E16.9,' [',A,']')
150 format('    ** Right value =>   : ',E16.9,' [',A,']')
160 format('    ** Start value      : ',E16.9,' [',A,']')
170 format('    ** Deviation        : ',E16.9,' [',A,']')
180 format('    ** Force constant   : ',E16.9,' [',A,']')
190 format('    ** Increment value  : ',E16.9,' [',A,']')
200 format('    ** Final value      : ',E16.9,' [',A,']')

355 format('    ** Min value        : ',E16.9,' [',A,']')
360 format('    ** Max value        : ',E16.9,' [',A,']')
365 format('    ** Number of bins   : ',I9)

end subroutine rst_restraints_rst_info

!===============================================================================
! Subroutine:  rst_restraints_increment
!===============================================================================

subroutine rst_restraints_increment(rst_item)

    use rst_dat
    use pmf_dat
    use pmf_cvs
    use pmf_paths

    implicit none
    type(CVTypeUM)     :: rst_item
    ! ----------------------------------------------
    real(PMFDP)        :: difference
    ! --------------------------------------------------------------------------

    select case(rst_item%mode)
        case('P')
            ! value changes programmatically by managed path
            rst_item%target_value = pmf_paths_get_rpos(rst_item%cvindx)
            return
        case('C')
            return  ! no change of constant restraint
        case('W')
            return  ! no change of wall restraint
        case('S')
            ! set corresponding value
            rst_item%target_value = rst_item%control_values(fstep)
            return
        case('I','V')
            ! incrementation is in linear mode
            difference = rst_item%cv%get_deviation(rst_item%stopvalue,rst_item%startvalue)
            rst_item%target_value = rst_item%startvalue + difference * fstep / fnstlim
        case default
    end select

end subroutine rst_restraints_increment

!===============================================================================
! Subroutine:  rst_restraints_read_control_file
!===============================================================================

subroutine rst_restraints_read_control_file(rst_item)

    use pmf_dat
    use pmf_utils
    use rst_dat

    implicit none
    type(CVTypeUM)     :: rst_item
    ! ----------------------------------------------
    integer            :: alloc_failed, ios
    integer            :: i
    ! --------------------------------------------------------------------------

    ! test if control file exists
    if( .not. pmf_utils_fexist(frstctr) ) then
        write(PMF_OUT,10) trim(frstctr)
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: frstctr file does not exist!')
    end if

    ! open control file
    call pmf_utils_open(RST_CTR,frstctr,'O')

    ! allocate rst_item%control_values
    allocate(rst_item%control_values(fnstlim), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
         call pmf_utils_exit(PMF_OUT,1,'[RST] Unable to allocate memory for control_values!')
    end if

    ! read data from control_file
    ios = 0
    do i=1,fnstlim
       read(RST_CTR,*,iostat=ios) rst_item%control_values(i)
       if ( ios .ne. 0 ) then
            write(PMF_OUT,20) trim(frstctr), i-1, fnstlim
            call pmf_utils_exit(PMF_OUT,1,'[RST] There is not enough data in frstctr file!')
       end if
    end do

    ! specify start value
    rst_item%startvalue = rst_item%control_values(1)

    close(RST_CTR)

    return

 10 format('[RST] frstctr file (',A,') does not exist!')
 20 format('[RST] frstctr file (',A,') has less number of data (',I8,') than fnstlim (',I8,')!')

end subroutine rst_restraints_read_control_file

!===============================================================================

end module rst_restraints
