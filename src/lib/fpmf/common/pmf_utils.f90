!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_utils

use pmf_sizes

implicit none

interface
    subroutine cpmf_print_errors()
    end subroutine cpmf_print_errors
end interface

contains

!===============================================================================
! Subroutine:  pmf_utils_open
!===============================================================================

subroutine pmf_utils_open(unitnum, filename, mystatus)

    implicit none
    integer           unitnum      ! logical unit number
    character(*)      filename     ! file name
    character(1)      mystatus     ! N - new, U - unknown, O - old, R - replace
    ! -----------------------------------------------
    integer           ierr
    character(7)      ustatus       ! status keyword
    character(7)      uposition
    ! --------------------------------------------------------------------------

    if( mystatus .eq. 'N' ) then
        ustatus = 'NEW'
        uposition = 'REWIND'
    else if( mystatus .eq. 'O' ) then
        ustatus = 'OLD'
        uposition = 'REWIND'
    else if( mystatus .eq. 'U' ) then
        ustatus = 'UNKNOWN'
        uposition = 'REWIND'
    else if( mystatus .eq. 'R' ) then
        ustatus = 'REPLACE'
        uposition = 'REWIND'
    else if( mystatus .eq. 'A' ) then
        ustatus = 'UNKNOWN'
        uposition = 'APPEND'
    else
        call pmf_utils_exit(6, 1,'Incorrect file status!')
    endif

    open(unit = unitnum, file = filename, status = ustatus, &
    position = uposition, form = 'FORMATTED', iostat = ierr)

    if( ierr .ne. 0 ) then
        write(6, '(/,a,a)') 'Unable to open file: ',filename
        call pmf_utils_exit(6, 1)
    endif

    return

end subroutine pmf_utils_open

!===============================================================================
! Subroutine:  pmf_utils_fexist
!===============================================================================

logical function pmf_utils_fexist(filename)

    use pmf_constants

    implicit none
    character(*)      filename     ! file name
    ! -----------------------------------------------
    integer           ierr
    ! --------------------------------------------------------------------------

    open(unit = PMF_TEST, file = filename, status = 'OLD', form = 'FORMATTED', iostat = ierr)

    pmf_utils_fexist = ierr .eq. 0

    if( pmf_utils_fexist ) close(PMF_TEST)

end function pmf_utils_fexist

!===============================================================================
! Subroutine:   pmf_utils_read_ctrl_real8
!===============================================================================

subroutine pmf_utils_read_ctrl_real8(prm_fin, name, value, fmt)

    use prmfile
    use pmf_constants

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(*)                        :: name
    real(PMFDP)                         :: value
    character(*)                        :: fmt
    ! --------------------------------------------
    logical                             :: defval
    character(len=80)                   :: buf1
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

10 format(A37,1X,'=',1X,A30)
15 format(A37,1X,'=',1X,A30,1X,'(default)')

end subroutine pmf_utils_read_ctrl_real8

!===============================================================================
! Subroutine:   pmf_utils_check_real8_in_range
!===============================================================================

subroutine pmf_utils_check_real8_in_range(sec,name,value,minv,maxv)

    use prmfile
    use pmf_constants

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    real(PMFDP)                         :: value
    real(PMFDP)                         :: minv
    real(PMFDP)                         :: maxv
    ! --------------------------------------------------------------------------

end subroutine pmf_utils_check_real8_in_range

!===============================================================================
! Subroutine:   pmf_utils_check_real8
!===============================================================================

subroutine pmf_utils_check_real8(sec,name,value,testv,cond)

    use prmfile
    use pmf_constants

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    real(PMFDP)                         :: value
    real(PMFDP)                         :: testv
    real(PMFDP)                         :: cond
    ! --------------------------------------------------------------------------

end subroutine pmf_utils_check_real8

!===============================================================================
! Subroutine:   pmf_utils_read_ctrl_integer
!===============================================================================

subroutine pmf_utils_read_ctrl_integer(prm_fin, name, value, fmt)

    use prmfile
    use pmf_constants

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(*)                        :: name
    integer                             :: value
    character(*)                        :: fmt
    ! --------------------------------------------
    logical                             :: defval
    character(len=80)                   :: buf1
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

10 format(A37,1X,'=',1X,A30)
15 format(A37,1X,'=',1X,A30,1X,'(default)')

end subroutine pmf_utils_read_ctrl_integer

!===============================================================================
! Subroutine:   pmf_utils_check_integer_in_range
!===============================================================================

subroutine pmf_utils_check_integer_in_range(sec,name,value,minv,maxv)

    use prmfile
    use pmf_constants

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    integer                             :: value
    integer                             :: minv
    integer                             :: maxv
    ! --------------------------------------------------------------------------

end subroutine pmf_utils_check_integer_in_range

!===============================================================================
! Subroutine:   pmf_utils_check_integer
!===============================================================================

subroutine pmf_utils_check_integer(sec,name,value,testv,cond)

    use prmfile
    use pmf_constants

    implicit none
    character(*)                        :: sec
    character(*)                        :: name
    integer                             :: value
    integer                             :: testv
    integer                             :: cond
    ! --------------------------------------------------------------------------

end subroutine pmf_utils_check_integer

!===============================================================================
! Subroutine:   pmf_utils_exit
!===============================================================================

subroutine pmf_utils_exit(unitnum, errcode, message)

    implicit none
    integer                :: unitnum
    integer                :: errcode
    character(*),optional  :: message
    ! --------------------------------------------------------------------------

    if(present(message)) then
        write(unitnum,'(/,A)') '>>> ERROR: ' // message
    end if

    call cpmf_print_errors

    write(unitnum,'(/,A)') '>>> ERROR: Some fatal error occured in PMFLib!'
    write(unitnum,'(A)')   '           Look above for detailed message (if any).'
    write(unitnum,'(A,/)') '           Program execution is terminated.'

    if( errcode .eq. 0 ) then
        stop 0
    else
        stop 1
    end if

end subroutine pmf_utils_exit

!===============================================================================
! Subroutine:   pmf_utils_heading
!===============================================================================

subroutine pmf_utils_heading(iunit, msg, fill)

    implicit none
    integer        :: iunit
    character(*)   :: msg
    character      :: fill
    ! -----------------------------------------------
    integer        :: n,m, i
    ! -------------------------------------------------------------------------

    if( mod(len(msg),2) .eq. 0 ) then
        n = (80 - len(msg) - 2) / 2
        m = n
    else
        n = (80 - len(msg) - 2) / 2
        m = n+1
    end if

    write(iunit,'(80a)') (fill, i=1,n), ' ', msg, ' ', (fill, i=1,m)

end subroutine pmf_utils_heading

!===============================================================================
! helper functions
!===============================================================================

real(PMFDP) function kdelta(i,j)

    implicit none
    integer             :: i
    integer             :: j
    ! --------------------------------------------------------------------------

    if( i .eq. j ) then
        kdelta = 1.0d0
    else
        kdelta = 0.0d0
    end if

end function kdelta

!===============================================================================

real(PMFDP) function vdot(x,y)

    implicit none
    real(PMFDP)     :: x(3)
    real(PMFDP)     :: y(3)
    ! --------------------------------------------------------------------------

    vdot = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)

end function vdot

!===============================================================================

function vcross(x,y)

    implicit none
    real(PMFDP)     :: x(3)
    real(PMFDP)     :: y(3)
    real(PMFDP)     :: vcross(3)
    ! --------------------------------------------------------------------------

    ! i j k i j
    ! 1 2 3 1 2
    ! 1 2 3 1 2

    ! cross-product
    vcross(1) = x(2)*y(3) - x(3)*y(2)
    vcross(2) = x(3)*y(1) - x(1)*y(3)
    vcross(3) = x(1)*y(2) - x(2)*y(1)

end function vcross

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_utils_header(progname)

    use pmf_ver
    use pmf_constants

    implicit none
    character(*)   :: progname
    ! -----------------------------------------------
    integer        :: datum(8)
    ! --------------------------------------------------------------------------

    ! get current date and time
    call date_and_time(values=datum)

    ! write header
    write (PMF_OUT,*)
    write (PMF_OUT,10)
    call pmf_utils_heading(PMF_OUT,'*** '//trim(progname)//' ***',' ')
    write (PMF_OUT,20)
    write (PMF_OUT,50) trim(PMFLIBVER)
    write (PMF_OUT,60)
    write (PMF_OUT,70) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)
    write (PMF_OUT,80)

    return

 10 format('|==============================================================================|')
 20 format('|------------------------------------------------------------------------------|')
 50 format('| Version : ',A66,' |')
 60 format('|------------------------------------------------------------------------------|')
 70 format('| Current date ',i4,'-',i2.2,'-',i2.2,' and time ',i2,':',i2.2,':',i2.2,'                                    |')
 80 format('|==============================================================================|')

end subroutine pmf_utils_header

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_utils_footer(progname)

    use pmf_constants

    implicit none
    character(*)   :: progname
    ! -----------------------------------------------
    integer        :: i
    ! --------------------------------------------------------------------------

    write (PMF_OUT,*)
    write (PMF_OUT,'(80a)') ('=',i=1,80)
    call pmf_utils_heading(PMF_OUT,trim(progname) // ' terminated normally.',' ')
    write (PMF_OUT,'(80a)') ('=',i=1,80)
    write (PMF_OUT,*)

    return

end subroutine pmf_utils_footer

!===============================================================================

end module pmf_utils


