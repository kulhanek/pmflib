!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

program pmf_gprocess

    use pmf_sizes
    use pmf_constants
    use pmf_utils

    character(len=PMF_MAX_PATH)  :: inname       ! input file name
    character(len=PMF_MAX_PATH)  :: outname      ! output file name
    integer                 :: size         ! number of data
    integer                 :: gpsize       ! number of valid data
    real(PMFDP),allocatable :: x(:)         ! x - values
    real(PMFDP),allocatable :: y(:)         ! y - values
    real(PMFDP),allocatable :: s(:)         ! standard deviations
    real(PMFDP),allocatable :: nsamples(:)  ! number of samples
    real(PMFDP),allocatable :: iy(:)        ! interpolated values
    real(PMFDP),allocatable :: xcov(:,:)
    real(PMFDP),allocatable :: kstar(:)
    real(PMFDP),allocatable :: values(:)
    real(PMFDP),allocatable :: alpha(:)
    integer,allocatable     :: indx(:)
    integer                 :: gpmin        ! minimum number of samples per bin
    real(PMFDP)             :: gplen        ! characteristic length-scale
    real(PMFDP)             :: gpsigmafac
    real(PMFDP)             :: gpsigmaoffset
    ! --------------------------------------------
    integer                 :: i
    ! --------------------------------------------------------------------------

    call pmf_utils_header('Gaussian process interpolation')

    ! test number of input arguments
    if( command_argument_count() .ne. 2 ) then
        call print_usage()
        call pmf_utils_exit(PMF_OUT,1,'Incorrect number of arguments was specified (two expected)!')
    end if

    call get_command_argument(1, inname)
    call get_command_argument(2, outname)

    ! load imput data
    call load_data()

    ! do gaussian process
    call gaussian_process

    ! calculate interpolated values
    iy = 0.0d0
    do i=1,size
       ! if( nsamples(i) .gt. gpmin ) then
            iy(i) = gprocess_interpolation(x(i))
       ! end if
    end do

    ! write result
    call save_data()

    call pmf_utils_footer('Gaussian process interpolation')

contains

!===============================================================================
! subroutine:  print_usage
!===============================================================================

subroutine print_usage()

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,'(/,a,/)') '=== [usage] ===================================================================='
    write(PMF_OUT,10)
    write(PMF_OUT,*)

    return

10 format('    pmf-gprocess <input> <output>')

end subroutine print_usage

!===============================================================================
! subroutine:  load_data
!===============================================================================

subroutine load_data()

    integer     :: i, alloc_failed
    ! --------------------------------------------------------------------------

    write(PMF_OUT,90) trim(inname)

    call pmf_utils_open(ABF_INP,inname,'O')

    ! read header
    read(ABF_INP,*,end=10,err=10) size
    read(ABF_INP,*,end=20,err=20) gpmin, gplen, gpsigmafac, gpsigmaoffset

    write(PMF_OUT,100) size
    write(PMF_OUT,110) gpmin
    write(PMF_OUT,120) gplen
    write(PMF_OUT,130) gpsigmafac
    write(PMF_OUT,140) gpsigmaoffset

    allocate(x(size),y(size),s(size),nsamples(size),iy(size), &
             xcov(size,size),kstar(size),values(size),alpha(size), &
             indx(size), &
             stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'Unable to allocate memory for the gaussian process!')
    end if

    do i=1,size
        read(ABF_INP,*,end=20,err=20) x(i), y(i), s(i), nsamples(i)
    end do

    close(ABF_INP)

    return

10  call pmf_utils_exit(PMF_OUT,1,'Unable to read input data header!')
15  call pmf_utils_exit(PMF_OUT,1,'Unable to read the gaussian process setup!')
20  call pmf_utils_exit(PMF_OUT,1,'Unable to read input data!')

 90 format('Input file name             = ',A)
100 format('Number of data points       = ',I6)
110 format('Minimum number of samples   = ',I6)
120 format('Characteristic length-scale = ',E13.6)
130 format('Sigma factor                = ',E13.6)
140 format('Sigma offset                = ',E13.6)

end subroutine load_data

!===============================================================================
! subroutine:  save_data
!===============================================================================

subroutine save_data()

    integer     :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,5) trim(outname)

    call pmf_utils_open(ABF_OUT,outname,'R')

    do i=1,size
        write(ABF_OUT,10,err=20) x(i),y(i),iy(i)
    end do

    close(ABF_OUT)

    return

 5 format('Output file name            = ',A)
10 format(E16.6,1X,E16.6,1X,E16.6)
20 call pmf_utils_exit(PMF_OUT,1,'Unable to save interpolated data!')

end subroutine save_data

!===============================================================================
! subroutine:  gaussian_process
!===============================================================================

subroutine gaussian_process()

    implicit none
    integer         :: i, j, ri, rj, info
    real(PMFDP)     :: xi, xj, ei, yi
    ! --------------------------------------------------------------------------

    ! determine size of problem

    gpsize=0
    do i=1,size
        if( nsamples(i) .gt. gpmin ) then
            gpsize = gpsize + 1
        end if
    end do

    write(PMF_OUT,100) gpsize

    if( gpsize .lt. 2 ) then
        call pmf_utils_exit(PMF_OUT,1,'Gaussian process requires at least two valid data points!')
    end if

    ! construct xcov
    xcov = 0
    ri=1
    do i=1,size
        if( nsamples(i) .lt. gpmin ) cycle
        xi = x(i)
        yi = y(i)
        ei = s(i)
        rj=1
        do j=1,size
            if( nsamples(j) .lt. gpmin ) cycle
            xj = x(j)
            xcov(ri,rj) = exp(-(xi-xj)**2/(2.0*gplen**2))
            if( ri .eq. rj )  xcov(ri,rj) = xcov(ri,rj) + ((ei + gpsigmaoffset) * gpsigmafac)**2
            rj = rj + 1
        end do
        values(ri) = yi
        ri = ri + 1
    end do

    ! LU decomposition
    call dgetrf(gpsize,gpsize,xcov,size,indx,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'LU decomposition failed in gaussian_process!')
    end if

    ! inversion
    call dgetri(gpsize,xcov,size,indx,alpha,gpsize,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Matrix inversion failed in gaussian_process!')
    end if

    ! matrix vector mult
    alpha = 0.0d0

    ! multiply by targets
    call dgemv('N',gpsize,gpsize,1.0d0,xcov,size,values,1,0.0d0,alpha,1)

100 format('Number of used data points  = ',I6)

end subroutine gaussian_process

!===============================================================================
! function:  gprocess_interpolation
!===============================================================================

real(PMFDP) function gprocess_interpolation(xv)

    implicit none
    real(PMFDP)     :: xv
    ! --------------------------------------------
    integer         :: i, ri
    real(PMFDP)     :: xi
    ! --------------------------------------------------------------------------

    ri=1
    do i=1,size
        if( nsamples(i) .lt. gpmin ) cycle
        xi = x(i)
        kstar(ri) = exp(-(xi-xv)**2/(2.0*gplen**2))
        ri = ri + 1
    end do

    gprocess_interpolation = 0.0
    do i=1,gpsize
        gprocess_interpolation = gprocess_interpolation + alpha(i)*kstar(i)
    end do

end function gprocess_interpolation

!===============================================================================

end program pmf_gprocess
