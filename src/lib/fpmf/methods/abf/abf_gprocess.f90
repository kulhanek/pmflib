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

module abf_gprocess

    use pmf_sizes
    use pmf_constants
    use pmf_utils

contains

!===============================================================================
! subroutine:  save_data
!===============================================================================

subroutine abf_gprocess_write_output()

    use abf_dat

    integer     :: i
    ! --------------------------------------------------------------------------

    if( fgpprint_period .le. 0 ) return ! output is written only of fgpprint_period > 0
    if( mod(fstep,fgpprint_period) .ne. 0 ) return

    call pmf_utils_open(ABF_GPOUT,fabfgpout,'R')

    do i=1,gpmaxsize
        write(ABF_GPOUT,10,err=20) xvalues(i),yvalues(i),abf_gprocess_interpol(xvalues(i))
    end do

    close(ABF_GPOUT)

    return

10 format(E16.6,1X,E16.6,1X,E16.6)
20 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to save interpolated data from gaussian process!')

end subroutine abf_gprocess_write_output

!===============================================================================
! subroutine:  abf_gprocess_update_model
!===============================================================================

subroutine abf_gprocess_update_model()

    use abf_dat

    implicit none
    integer         :: i, j, ri, rj, info
    real(PMFDP)     :: xi, xj, ei, yi, gplen, gpsoffset, gpsfac
    ! --------------------------------------------------------------------------

    if( fgpmodel_update .le. 0 ) return ! output is written only of fgpprint_period > 0
    if( mod(fstep,fgpmodel_update) .ne. 0 ) return

    ! determine gpsize of problem
    gpsize=0
    do i=1,gpmaxsize
        if( accumulator%nsamples(i) .gt. fgpmin_samples ) then
            gpsize = gpsize + 1
        end if
    end do

    write(ABF_OUT,100) fstep, gpsize

    if( gpsize .lt. 2 ) then
        return
    end if

    ! update data
    do i=1,gpmaxsize
        yvalues(i) = accumulator%abfforce(1,i) / real(accumulator%nsamples(i))
        ei = accumulator%nsamples(i)*accumulator%abfforce2(1,i) - &
             accumulator%abfforce(1,i)*accumulator%abfforce(1,i)
        if( ei .gt. 0 ) then
            ei = sqrt(ei) / accumulator%nsamples(i);
        else
            ei = 0.0
        end if
        svalues(i) = ei / sqrt( real(accumulator%nsamples(i)) )
    end do

    gplen = ABFCVList(1)%fgplen
    gpsoffset = ABFCVList(1)%fgpsigmaoffset
    gpsfac = ABFCVList(1)%fgpsigmafac

    ! construct xcov
    xcov = 0
    ri=1
    do i=1,gpmaxsize
        if( accumulator%nsamples(i) .lt. fgpmin_samples ) cycle
        xi = xvalues(i)
        yi = yvalues(i)
        ei = svalues(i)
        rj=1
        do j=1,gpmaxsize
            if( accumulator%nsamples(j) .lt. fgpmin_samples ) cycle
            xj = xvalues(j)
            xcov(ri,rj) = exp(-(xi-xj)**2/(2.0*gplen**2))
            if( ri .eq. rj )  xcov(ri,rj) = xcov(ri,rj) + ((ei + gpsoffset)* gpsfac)**2
            rj = rj + 1
        end do
        gpyvalues(ri) = yi
        ri = ri + 1
    end do

    ! LU decomposition
    call dgetrf(gpsize,gpsize,xcov,gpmaxsize,gpindx,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'LU decomposition failed in gaussian_process!')
    end if

    ! inversion
    call dgetri(gpsize,xcov,gpmaxsize,gpindx,alpha,gpsize,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Matrix inversion failed in gaussian_process!')
    end if

    ! matrix vector mult
    alpha = 0.0d0

    ! multiply by targets
    call dgemv('N',gpsize,gpsize,1.0d0,xcov,gpmaxsize,gpyvalues,1,0.0d0,alpha,1)

100 format('# [GP]   Total steps     = ',I12,' Number of data points used = ',I12)

end subroutine abf_gprocess_update_model

!===============================================================================
! function:  abf_gprocess_interpol
!===============================================================================

real(PMFDP) function abf_gprocess_interpol(xv)

    use abf_dat

    implicit none
    real(PMFDP)     :: xv
    ! --------------------------------------------
    integer         :: i, ri
    real(PMFDP)     :: xi,gplen
    ! --------------------------------------------------------------------------

    gplen = ABFCVList(1)%fgplen

    ri=1
    do i=1,gpmaxsize
        if( accumulator%nsamples(i) .lt. fgpmin_samples ) cycle
        xi = xvalues(i)
        kstar(ri) = exp(-(xi-xv)**2/(2.0*gplen**2))
        ri = ri + 1
    end do

    abf_gprocess_interpol = 0.0
    do i=1,gpsize
        abf_gprocess_interpol = abf_gprocess_interpol + alpha(i)*kstar(i)
    end do

end function abf_gprocess_interpol

!===============================================================================

end module abf_gprocess
