!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module pmf_control_utils

implicit none
contains

!===============================================================================
! Subroutine:   pmf_ctrl_print_default_logical
!===============================================================================

subroutine pmf_ctrl_print_default_logical(name, value)

    use prmfile
    use pmf_constants

    implicit none
    character(*)                        :: name
    logical                             :: value
    ! --------------------------------------------
    character(len=37)                   :: buf2
    ! --------------------------------------------------------------------------

    buf2 = name

    ! setup format string
    write(PMF_OUT,15) adjustl(buf2), prmfile_onoff(value)

15 format(A37,1X,'=',1X,A12,19X,'(default)')

end subroutine pmf_ctrl_print_default_logical

!===============================================================================
! Subroutine:   pmf_ctrl_print_default_real8
!===============================================================================

subroutine pmf_ctrl_print_default_real8( name, value, fmt)

    use prmfile
    use pmf_constants

    implicit none
    character(*)                        :: name
    real(PMFDP)                         :: value
    character(*)                        :: fmt
    ! --------------------------------------------
    character(len=37)                   :: buf2
    character(len=30)                   :: buf3
    ! --------------------------------------------------------------------------

    buf2 = name
    write(buf3,'('//trim(fmt)//')') value

    ! setup format string
    write(PMF_OUT,15) adjustl(buf2), buf3

15 format(A37,1X,'=',1X,A12,19X,'(default)')

end subroutine pmf_ctrl_print_default_real8

!===============================================================================
! Subroutine:   pmf_ctrl_print_default_integer
!===============================================================================

subroutine pmf_ctrl_print_default_integer(name, value, fmt)

    use prmfile
    use pmf_constants

    implicit none
    character(*)                        :: name
    integer                             :: value
    character(*)                        :: fmt
    ! --------------------------------------------
    character(len=37)                   :: buf2
    character(len=30)                   :: buf3
    ! --------------------------------------------------------------------------

    buf2 = name
    write(buf3,'('//trim(fmt)//')') value

    ! setup format string
    write(PMF_OUT,15) adjustl(buf2), buf3

15 format(A37,1X,'=',1X,A12,19X,'(default)')

end subroutine pmf_ctrl_print_default_integer

!===============================================================================
! Subroutine:  pmf_ctrl_print_default_stritem
!===============================================================================

subroutine pmf_ctrl_print_default_stritem(sname,svalue)

    use pmf_constants

    implicit none
    character(*)                        :: sname
    character(*)                        :: svalue
    character(len=27)                   :: buff
    ! --------------------------------------------------------------------------

    write(buff,5) trim(sname)
    write(PMF_OUT,15) adjustl(buff),trim(svalue)

 5 format(A27)
15 format(A27,' = ',A40,' (default)')

end subroutine pmf_ctrl_print_default_stritem

!===============================================================================
! Subroutine:  pmf_ctrl_read_stritem
!===============================================================================

subroutine pmf_ctrl_read_stritem(prm_fin,sname,svalue)

    use prmfile
    use pmf_constants

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(*)                        :: sname
    character(*)                        :: svalue
    character(len=27)                   :: buff
    ! --------------------------------------------------------------------------

    write(buff,5) trim(sname)
    if( prmfile_get_string_by_key(prm_fin,sname,svalue) ) then
        write(PMF_OUT,10) adjustl(buff),trim(svalue)
    else
        write(PMF_OUT,15) adjustl(buff),trim(svalue)
    end if

 5 format(A27)
10 format(A27,' = ',A40)
15 format(A27,' = ',A40,' (default)')

end subroutine pmf_ctrl_read_stritem

!===============================================================================
! Subroutine:   pmf_ctrl_read_logical
!===============================================================================

subroutine pmf_ctrl_read_logical(prm_fin, name, value)

    use prmfile
    use pmf_constants

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(*)                        :: name
    logical                             :: value
    ! --------------------------------------------
    logical                             :: defval
    character(len=37)                   :: buf2
    ! --------------------------------------------------------------------------

    ! read value
    defval = prmfile_get_logical_by_key(prm_fin,name,value)

    buf2 = name

    ! setup format string
    if( defval ) then
        write(PMF_OUT,10) adjustl(buf2), prmfile_onoff(value)
    else
        write(PMF_OUT,15) adjustl(buf2), prmfile_onoff(value)
    end if

10 format(A37,1X,'=',1X,A12)
15 format(A37,1X,'=',1X,A12,19X,'(default)')

end subroutine pmf_ctrl_read_logical

!===============================================================================
! Subroutine:   pmf_ctrl_read_real8
!===============================================================================

subroutine pmf_ctrl_read_real8(prm_fin, name, value, fmt)

    use prmfile
    use pmf_constants

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(*)                        :: name
    real(PMFDP)                         :: value
    character(*)                        :: fmt
    ! --------------------------------------------
    logical                             :: defval
    character(len=37)                   :: buf2
    character(len=30)                   :: buf3
    ! --------------------------------------------------------------------------

    ! read value
    defval = prmfile_get_real8_by_key(prm_fin,name,value)

    buf2 = name
    write(buf3,'('//trim(fmt)//')') value

    ! setup format string
    if( defval ) then
        write(PMF_OUT,10) adjustl(buf2), buf3
    else
        write(PMF_OUT,15) adjustl(buf2), buf3
    end if

10 format(A37,1X,'=',1X,A12)
15 format(A37,1X,'=',1X,A12,19X,'(default)')

end subroutine pmf_ctrl_read_real8

!===============================================================================
! Subroutine:   pmf_ctrl_check_real8
!===============================================================================

subroutine pmf_ctrl_check_real8(sec,name,value,testv,cond,fmt)

    use prmfile
    use pmf_constants
    use pmf_unit
    use pmf_utils

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    real(PMFDP)                         :: value
    real(PMFDP)                         :: testv
    integer                             :: cond
    character(*)                        :: fmt
    ! -------------------------------------------
    character(len=PMF_MAX_PATH)         :: buffer
    character(len=37)                   :: buf1,buf2
    ! --------------------------------------------------------------------------

    write(buf1,'('//trim(fmt)//')') value
    write(buf2,'('//trim(fmt)//')') testv

    select case(cond)
        case(CND_GT)
            if( value .le. testv ) then
                write(buffer,10) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_GE)
            if( value .lt. testv ) then
                write(buffer,20) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
         case(CND_LT)
            if( value .ge. testv ) then
                write(buffer,30) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_LE)
            if( value .gt. testv ) then
                write(buffer,40) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_EQ)
            if( value .ne. testv ) then
                write(buffer,50) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_NE)
            if( value .eq. testv ) then
                write(buffer,60) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
    end select

 10 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be grater than ',A)
 20 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be grater than or equal to ',A)
 30 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be lesser than ',A)
 40 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be lesser than or equal to ',A)
 50 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to equal to ',A)
 60 format('[',A,']',1X,A,1X,'parameter: ',A,' value does not have to be ',A)

end subroutine pmf_ctrl_check_real8

!===============================================================================
! Subroutine:   pmf_ctrl_read_integer
!===============================================================================

subroutine pmf_ctrl_read_integer(prm_fin, name, value, fmt)

    use prmfile
    use pmf_constants

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(*)                        :: name
    integer                             :: value
    character(*)                        :: fmt
    ! --------------------------------------------
    logical                             :: defval
    character(len=37)                   :: buf2
    character(len=30)                   :: buf3
    ! --------------------------------------------------------------------------

    ! read value
    defval = prmfile_get_integer_by_key(prm_fin,name,value)

    buf2 = name
    write(buf3,'('//trim(fmt)//')') value

    ! setup format string
    if( defval ) then
        write(PMF_OUT,10) adjustl(buf2), buf3
    else
        write(PMF_OUT,15) adjustl(buf2), buf3
    end if

10 format(A37,1X,'=',1X,A12)
15 format(A37,1X,'=',1X,A12,19X,'(default)')

end subroutine pmf_ctrl_read_integer

!===============================================================================
! Subroutine:   pmf_ctrl_check_integer_in_range
!===============================================================================

subroutine pmf_ctrl_check_integer_in_range(sec,name,value,minv,maxv)

    use prmfile
    use pmf_constants
    use pmf_utils

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    integer                             :: value
    integer                             :: minv
    integer                             :: maxv
    ! -------------------------------------------
    character(len=PMF_MAX_PATH)         :: buffer
    character(len=37)                   :: buf1,buf2,buf3
    ! --------------------------------------------------------------------------

    write(buf1,'(I12)') value
    write(buf2,'(I12)') minv
    write(buf3,'(I12)') maxv

    if( (value .lt. minv) .or. (value .gt. maxv) ) then
        write(buffer,10) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2)),trim(adjustl(buf3))
        call pmf_utils_exit(PMF_OUT,1,trim(buffer))
    end if

 10 format('[',A,']',1X,A,1X,'parameter: ',A,' value is out-of-allowed range: ',A,'-',A)

end subroutine pmf_ctrl_check_integer_in_range

!===============================================================================
! Subroutine:   pmf_ctrl_check_integer
!===============================================================================

subroutine pmf_ctrl_check_integer(sec,name,value,testv,cond)

    use prmfile
    use pmf_constants
    use pmf_utils

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    integer                             :: value
    integer                             :: testv
    integer                             :: cond
    ! -------------------------------------------
    character(len=PMF_MAX_PATH)         :: buffer
    character(len=37)                   :: buf1,buf2
    ! --------------------------------------------------------------------------

    write(buf1,'(I12)') value
    write(buf2,'(I12)') testv

    select case(cond)
        case(CND_GT)
            if( value .le. testv ) then
                write(buffer,10) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_GE)
            if( value .lt. testv ) then
                write(buffer,20) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
         case(CND_LT)
            if( value .ge. testv ) then
                write(buffer,30) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_LE)
            if( value .gt. testv ) then
                write(buffer,40) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_EQ)
            if( value .ne. testv ) then
                write(buffer,50) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_NE)
            if( value .eq. testv ) then
                write(buffer,60) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
    end select

 10 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be grater than ',A)
 20 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be grater than or equal to ',A)
 30 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be lesser than ',A)
 40 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be lesser than or equal to ',A)
 50 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to equal to ',A)
 60 format('[',A,']',1X,A,1X,'parameter: ',A,' value does not have to be ',A)

end subroutine pmf_ctrl_check_integer

!===============================================================================
! Subroutine:   pmf_ctrl_read_real8_wunit
!===============================================================================

subroutine pmf_ctrl_read_real8_wunit(prm_fin, name, iunit, value, fmt)

    use prmfile
    use pmf_constants
    use pmf_unit
    use pmf_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(*)                        :: name
    type(UnitType)                      :: iunit
    real(PMFDP)                         :: value
    character(*)                        :: fmt
    ! --------------------------------------------
    logical                             :: defval
    character(len=37)                   :: buf2
    character(len=12)                   :: buf3
    character(len=17)                   :: buf4
    ! --------------------------------------------------------------------------

    ! read value
    defval = prmfile_get_real8_by_key(prm_fin,name,value)

    call pmf_unit_conv_to_ivalue(iunit,value)

    buf2 = name
    write(buf3,'('//trim(fmt)//')') pmf_unit_get_rvalue(iunit,value)
    write(buf4,20) trim(pmf_unit_label(iunit))

    ! setup format string
    if( defval ) then
        write(PMF_OUT,10) adjustl(buf2), adjustr(buf3), adjustl(buf4)
    else
        write(PMF_OUT,15) adjustl(buf2), adjustr(buf3), adjustl(buf4)
    end if

10 format(A37,1X,'=',1X,A12,1X,A17)
15 format(A37,1X,'=',1X,A12,1X,A17,1X,'(default)')
20 format('[',A,']')

end subroutine pmf_ctrl_read_real8_wunit

!===============================================================================
! Subroutine:   pmf_ctrl_check_real8_wunit
!===============================================================================

subroutine pmf_ctrl_check_real8_wunit(sec,name, iunit,value,testv,cond,fmt)

    use prmfile
    use pmf_constants
    use pmf_unit
    use pmf_utils

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    type(UnitType)                      :: iunit
    real(PMFDP)                         :: value
    real(PMFDP)                         :: testv
    integer                             :: cond
    character(*)                        :: fmt
    ! -------------------------------------------
    character(len=PMF_MAX_PATH)         :: buffer
    real(PMFDP)                         :: ravlue,rtestv
    character(len=37)                   :: buf1,buf2
    ! --------------------------------------------------------------------------

    ravlue = pmf_unit_get_rvalue(iunit,value)
    rtestv = pmf_unit_get_rvalue(iunit,testv)

    write(buf1,'('//trim(fmt)//')') pmf_unit_get_rvalue(iunit,ravlue)
    write(buf2,'('//trim(fmt)//')') pmf_unit_get_rvalue(iunit,rtestv)

    select case(cond)
        case(CND_GT)
            if( value .le. testv ) then
                write(buffer,10) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2)),trim(pmf_unit_label(iunit))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_GE)
            if( value .lt. testv ) then
                write(buffer,20) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2)),trim(pmf_unit_label(iunit))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
         case(CND_LT)
            if( value .ge. testv ) then
                write(buffer,30) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2)),trim(pmf_unit_label(iunit))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_LE)
            if( value .gt. testv ) then
                write(buffer,40) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2)),trim(pmf_unit_label(iunit))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_EQ)
            if( value .ne. testv ) then
                write(buffer,50) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2)),trim(pmf_unit_label(iunit))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
        case(CND_NE)
            if( value .eq. testv ) then
                write(buffer,60) trim(sec),trim(name),trim(adjustl(buf1)),trim(adjustl(buf2)),trim(pmf_unit_label(iunit))
                call pmf_utils_exit(PMF_OUT,1,trim(buffer))
            end if
    end select

 10 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be grater than ',A              ,1X,'[',A,']')
 20 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be grater than or equal to ',A  ,1X,'[',A,']')
 30 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be lesser than ',A              ,1X,'[',A,']')
 40 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to be lesser than or equal to ',A  ,1X,'[',A,']')
 50 format('[',A,']',1X,A,1X,'parameter: ',A,' value has to equal to ',A                    ,1X,'[',A,']')
 60 format('[',A,']',1X,A,1X,'parameter: ',A,' value does not have to be ',A                ,1X,'[',A,']')

end subroutine pmf_ctrl_check_real8_wunit

!===============================================================================

end module pmf_control_utils

