! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License along
!     with this program; if not, write to the Free Software Foundation, Inc.,
!     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ==============================================================================

program test_stat

    use pmf_sizes
    use pmf_constants
    use pmf_utils

    implicit none
    integer,parameter           ::  RPREC=16

    integer                     ::  command_argument_count
    character(PMF_MAX_PATH)     ::  SizeString
    character(PMF_MAX_PATH)     ::  InputFile
    real(RPREC),allocatable     ::  x(:)
    real(RPREC),allocatable     ::  y(:)
    integer                     ::  DataSize,i
    ! --------------------------------------------------------------------------

    call pmf_utils_header('Test Data Statistics')

! check number of arguments ---------------------
    if( command_argument_count() .eq. 0 ) then
        call print_usage
        write(PMF_OUT,*)
        stop
    end if

    if( command_argument_count() .ne. 2 ) then
        call print_usage
        call pmf_utils_exit(PMF_OUT,1,'Incorrect number of arguments was specified (two expected)!')
    end if

! get arguments and print -----------------------
    call get_command_argument(1, SizeString)
    call get_command_argument(2, InputFile)

! allocate and read data
    read(SizeString,*) DataSize
    allocate(x(DataSize),y(DataSize))
    call pmf_utils_open(CST_INP,InputFile,'O')
    do i=1,DataSize
        read(CST_INP,*) x(i), y(i)
    end do
    close(CST_INP)

    write(PMF_OUT,*)
    write(PMF_OUT,10) DataSize
    write(PMF_OUT,20) RPREC

! do statistics
    call naive_stat
    call online_stat

! print header ----------------------------------
    call pmf_utils_footer('Test Data Statistics')

    stop

10 format('# Data size: ',I10)
20 format('# Precision: ',I10,' bytes')

!===============================================================================
contains

!===============================================================================
! Subroutine: print_usage
!===============================================================================

subroutine print_usage

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a,/)') '=== [usage] ===================================================================='

    write(PMF_OUT,10)

    return

10 format('    pmf-test-stat <size> <xy-data')

end subroutine print_usage

!===============================================================================
! Subroutine: naive_stat
!===============================================================================

subroutine naive_stat

    implicit none
    real(RPREC)     :: sum_x
    real(RPREC)     :: sum_y
    real(RPREC)     :: sum_xy
    real(RPREC)     :: sum_x2
    real(RPREC)     :: sum_y2
    real(RPREC)     :: n,dx,dy,mx,my
    real(RPREC)     :: sum_dxy
    real(RPREC)     :: sum_dxy2,mdxy,sum_sxy
    integer         :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,5)

    sum_x   = 0.0d0
    sum_y   = 0.0d0
    sum_xy  = 0.0d0
    sum_x2  = 0.0d0
    sum_y2  = 0.0d0
    n       = 0.0d0

    do i=1,DataSize
        sum_x   = sum_x  + x(i)
        sum_y   = sum_y  + y(i)
        sum_xy  = sum_xy + x(i)*y(i)
        sum_x2  = sum_x2 + x(i)*x(i)
        sum_y2  = sum_y2 + y(i)*y(i)
        n = n + 1
    end do

    mx = sum_x/n
    my = sum_y/n

    write(PMF_OUT,10) mx
    write(PMF_OUT,15) sqrt((sum_x2-sum_x*sum_x/n)/n)
    write(PMF_OUT,20) my
    write(PMF_OUT,25) sqrt((sum_y2-sum_y*sum_y/n)/n)
    write(PMF_OUT,30) (sum_xy - sum_x*sum_y/n)/n

    write(PMF_OUT,*)
    write(PMF_OUT,105)

    sum_dxy     = 0.0d0
    sum_dxy2    = 0.0d0
    n           = 0.0d0
    do i=1,DataSize
        dx          = x(i) - mx
        dy          = y(i) - my
        sum_dxy     = sum_dxy + dx*dy
        sum_dxy2    = sum_dxy2 + (dx*dy)**2
        n = n + 1
    end do

    write(PMF_OUT,130) sum_dxy/n
    write(PMF_OUT,135) sqrt((sum_dxy2-sum_dxy*sum_dxy/n)/n)

    write(PMF_OUT,*)
    write(PMF_OUT,205)

    mdxy        = sum_dxy/n
    sum_sxy     = 0.0d0
    n           = 0.0d0
    do i=1,DataSize
        dx          = x(i) - mx
        dy          = y(i) - my
        sum_sxy     = sum_sxy + (dx*dy - mdxy)**2
        n = n + 1
    end do

    write(PMF_OUT,230) mdxy
    write(PMF_OUT,235) sqrt(sum_sxy/n)

  5 format('== Naive algorithms')
 10 format('<X>         = ',F26.16)
 15 format('s(X)        = ',F26.16)
 20 format('<Y>         = ',F26.16)
 25 format('s(Y)        = ',F26.16)
 30 format('cov(XY)     = ',F26.16)

105 format('== Two-pass algorithms')
130 format('cov(XY)     = ',F26.16)
135 format('s(cov(XY))  = ',F26.16)

205 format('== Two-pass+two-pass algorithms')
230 format('cov(XY)     = ',F26.16)
235 format('s(cov(XY))  = ',F26.16)

end subroutine naive_stat

!===============================================================================
! Subroutine: online_stat
!===============================================================================

subroutine online_stat

    implicit none
    real(RPREC)     :: mx
    real(RPREC)     :: my
    real(RPREC)     :: m2x
    real(RPREC)     :: m2y
    real(RPREC)     :: cxy
    real(RPREC)     :: n,dx,dy,dx2,dy2
    integer         :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,5)

    mx      = 0.0d0
    my      = 0.0d0
    m2x     = 0.0d0
    m2y     = 0.0d0
    cxy     = 0.0d0
    n       = 0.0d0

    do i=1,DataSize
        n       = n + 1
        dx      = x(i) - mx
        mx      = mx  + dx/n
        dx2     = x(i) - mx
        m2x     = m2x + dx * dx2

        dy      = y(i) - my
        my      = my  + dy/n
        dy2     = y(i) - my
        m2y     = m2y + dy * dy2

        cxy     = cxy + dx * dy2
    end do

    write(PMF_OUT,10) mx
    write(PMF_OUT,15) sqrt(m2x/n)
    write(PMF_OUT,20) my
    write(PMF_OUT,25) sqrt(m2y/n)
    write(PMF_OUT,30) cxy/n
    write(PMF_OUT,35) sqrt(m2x/n) * sqrt(m2y/n)

  5 format('== Online algorithms')
 10 format('<X>         = ',F26.16)
 15 format('s(X)        = ',F26.16)
 20 format('<Y>         = ',F26.16)
 25 format('s(Y)        = ',F26.16)
 30 format('cov(XY)     = ',F26.16)
 35 format('s(cov(XY))  = ',F26.16,' approx = s(X)*s(Y)')

end subroutine online_stat

!===============================================================================

end program test_stat


