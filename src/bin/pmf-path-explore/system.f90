! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_system

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine init_timers

    use smf_profiling_dat
    use smf_profiling
    use pmf_system_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! add standard timers --------------------------------
    INITIALIZATION_TIMER   = add_timer(TOTAL_TIMER,'Initialization')
        DERIVATIVES_TIMER      = add_timer(INITIALIZATION_TIMER,'FES Derivatives')
    CORE_TIMER             = add_timer(TOTAL_TIMER,'Program core')
        FORCES_TIMER           = add_timer(CORE_TIMER,'Forces')
        UPDATE_TIMER           = add_timer(CORE_TIMER,'STM Updates')
        SMOOTH_TIMER           = add_timer(CORE_TIMER,'STM Smoothing')
        REPARAM_TIMER          = add_timer(CORE_TIMER,'STM Reparametrization')
    FINALIZATION_TIMER     = add_timer(TOTAL_TIMER,'Finalization')

end subroutine init_timers

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine calc_fes_deriv

    use pmf_system_dat

    implicit none
    integer            :: i,j,findex
    integer            :: ib(ncvs)
    integer            :: it(ncvs)
    real(PMFDP)        :: ty1,ty2
    ! --------------------------------------------------------------------------

    do i=1,tot_nbins
        ! determine ib
        findex = i - 1
        do j=ncvs,1,-1
            ib(j) = mod(findex,cvs(j)%nbins)
            findex = findex / cvs(j)%nbins
        end do

        ! now determine derivatives in all directions
        do j=1,ncvs
            ! point 1
            it = ib
            it(j) = it(j) - 1
            if( it(j) .lt. 0 ) cycle
            ty1 = get_fes_value(it)

            ! point 2
            it = ib
            it(j) = it(j) + 1
            if( it(j) .ge. cvs(j)%nbins ) cycle
            ty2 = get_fes_value(it)

            dfes(j,i) = (ty2 - ty1)/(2.0d0*cvs(j)%width)
        end do
    end do

end subroutine calc_fes_deriv

!-------------------------------------------------------------------------------

real(PMFDP) function get_fes_value(ib)

    use pmf_system_dat

    implicit none
    integer            :: ib(ncvs)
    ! -----------------------------------------------
    integer            :: i,findex
    ! --------------------------------------------------------------------------

    findex = 0
    do i=1,ncvs
        findex = findex*cvs(i)%nbins + ib(i)
        end do

        findex = findex + 1

        if( (findex .lt. 1) .or. (findex .gt. tot_nbins) ) then
            get_fes_value = 0.0d0
            return
        end if

        get_fes_value = fes(findex)
    return

end function get_fes_value

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine open_trajectory_file

    use pmf_system_dat
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    ! do we print trajectory?
    if( traj_freq .le. 0 ) return

    call pmf_utils_open(IO_TRJ,TrajectoryFile,'R')

    write(IO_TRJ,10) ncvs, nbeads
    call write_path_summary_header(IO_TRJ)

    return

10 format('# STMTRAJ ',I3,1X,I3)

end subroutine open_trajectory_file

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_trajectory_snap

    use pmf_system_dat
    use pmf_core

    implicit none
    ! --------------------------------------------------------------------------

    if( .not. ((traj_freq .gt. 0) .and. (mod(istep,traj_freq) .eq. 0)) ) return

    write(IO_TRJ,10) istep/traj_freq
    call write_path_summary(IO_TRJ)
    write(IO_TRJ,*)

    return

10 format('# STMSNAP ',I8)

end subroutine write_trajectory_snap

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine close_trajectory_file

    use pmf_system_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! do we print trajectory?
    if( traj_freq .le. 0 ) return

    close(IO_TRJ)
    return

end subroutine close_trajectory_file

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_path_summary_header(ioout)

    use pmf_system_dat
    use pmf_sizes

    implicit none
    integer     :: ioout
    ! --------------------------------------------
    integer                         :: i,id
    character(len=PMF_MAX_CV_NAME)  :: buffer
    ! --------------------------------------------------------------------------

    write(ioout,10)
    write(ioout,20) trim(path_item%name)
    write(ioout,30) ncvs
    write(ioout,40) nbeads

! header --------------------
    ! legends
    write(ioout,100,ADVANCE='NO')
    do i=1,ncvs
        write(ioout,101,ADVANCE='NO') i
    end do
    do i=1,ncvs
        write(ioout,102,ADVANCE='NO') i
    end do
    do i=1,ncvs
        write(ioout,103,ADVANCE='NO') i
    end do
    do i=1,ncvs
        write(ioout,104,ADVANCE='NO') i
    end do
    write(ioout,*)

    ! delimiters
    write(ioout,300,ADVANCE='NO')
    do i=1,4*ncvs
        write(ioout,310,ADVANCE='NO')
    end do
    write(ioout,*)

! data ----------------------
    write(ioout,200,ADVANCE='NO')
    do i=1,ncvs
        buffer = CVList(i)%cv%name
        buffer = adjustr(buffer)
        write(ioout,130,ADVANCE='NO') buffer
    end do
    do i=1,ncvs
        buffer = CVList(i)%cv%name
        buffer = adjustr(buffer)
        write(ioout,130,ADVANCE='NO') buffer
    end do
    do i=1,ncvs
        buffer = CVList(i)%cv%name
        buffer = adjustr(buffer)
        write(ioout,130,ADVANCE='NO') buffer
    end do
    do i=1,ncvs
        buffer = CVList(i)%cv%name
        buffer = adjustr(buffer)
        write(ioout,130,ADVANCE='NO') buffer
    end do
    write(ioout,*)

    write(ioout,210,ADVANCE='NO')
    do i=1,ncvs
        buffer = CVList(i)%cv%ctype
        buffer = adjustr(buffer)
        write(ioout,130,ADVANCE='NO') buffer
    end do
    write(ioout,*)

    write(ioout,220,ADVANCE='NO')
    do i=1,ncvs
        write(ioout,110,ADVANCE='NO') path_item%minvalues(i)
    end do
    write(ioout,*)

    write(ioout,230,ADVANCE='NO')
    do i=1,ncvs
        write(ioout,110,ADVANCE='NO') path_item%maxvalues(i)
    end do
    write(ioout,*)

    write(ioout,240,ADVANCE='NO')
    do i=1,ncvs
        if( path_item%maxmoves(i) .le. 0 ) then
            write(ioout,120,ADVANCE='NO')
        else
            write(ioout,110,ADVANCE='NO') path_item%maxmoves(i)
        end if
    end do
    write(ioout,*)

    write(ioout,300,ADVANCE='NO')
    do i=1,4*ncvs
        write(ioout,310,ADVANCE='NO')
    end do
    write(ioout,*)

    write(ioout,305,ADVANCE='NO')
    id = 9
    do i=1,4*ncvs
        write(ioout,320,ADVANCE='NO') id
        id = id + 1
    end do
    write(ioout,*)

    write(ioout,300,ADVANCE='NO')
    do i=1,4*ncvs
        write(ioout,310,ADVANCE='NO')
    end do
    write(ioout,*)

 10 format('# === [PATH] ===================================================================')
 20 format('# Path name       = ',A)
 30 format('# Number of CVs   = ',I6)
 40 format('# Number of beads = ',I6)

100 format('#  ID   Type  ST  alpha  dA/dalpha       A        CID Updates')
101 format('     CV',I2.2,'    ')
102 format('   dA/dCV',I2.2,'  ')
103 format(' dCV',I2.2,'/dalpha')
104 format(' -|F',I2.2,'/dalpha')

110 format(1X,E12.5)
120 format(1X,'     --     ')
130 format(1X,A12)

200 format('#      names                                                 ');
210 format('#      types                                                 ');
220 format('#      min                                                   ');
230 format('#      max                                                   ');
240 format('#      maxmov                                                ')

300 format('# ---- ------ -- ------ ------------ ------------ --- -------')
305 format('#    1      2  3      4            5            6   7       8')
310 format(' ------------')
320 format(1X,I12)

end subroutine write_path_summary_header

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_path_summary(iout)

    use pmf_system_dat

    implicit none
    integer     :: iout
    ! --------------------------------------------
    integer     :: i,b
    REAL(PMFDP) :: tau(ncvs)
    ! --------------------------------------------------------------------------

    do b=1,nbeads
        if( path_item%fixed(b) ) then
            write(iout,10,ADVANCE='NO') b
        else
            write(iout,20,ADVANCE='NO') b
        end if

        call pmf_paths_get_intpoint_der(path_item,path_item%alphas(b),tau)
        write(iout,30,ADVANCE='NO')
        write(iout,40,ADVANCE='NO') path_item%alphas(b)
        write(iout,50,ADVANCE='NO') dot_product(deriv(:,b),tau(:))
        write(iout,50,ADVANCE='NO') values(b)
        write(iout,60,ADVANCE='NO')
        write(iout,70,ADVANCE='NO') 0

        do i=1,ncvs
            write(iout,50,ADVANCE='NO') path_item%points(b,i)
        end do
        do i=1,ncvs
            write(iout,50,ADVANCE='NO') deriv(i,b)
        end do
        do i=1,ncvs
            write(iout,50,ADVANCE='NO') tau(i)
        end do
        write(iout,*)
    end do

 10 format('  ',I4,' P      ')
 20 format('  ',I4,' F      ')
 30 format('UN')
 40 format(1X,F6.3)
 50 format(1X,F12.5)
 60 format('  --')
 70 format(1X,I7)

end subroutine write_path_summary

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_results

    use pmf_paths
    use pmf_system_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! summary file
    call pmf_utils_open(IO_SUM,OutputPathSummaryFile,'R')
    call write_path_summary_header(IO_SUM)
    call write_path_summary(IO_SUM)
    close(IO_SUM)

    ! path
    call pmf_utils_open(IO_PATH,OutputPathFile,'R')
    call pmf_paths_write_path(IO_PATH,path_item)
    close(IO_PATH)

    return

end subroutine write_results

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_stmout_header

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,10)
    write(PMF_OUT,20)

    return

10 format('# Step   Path length  Length change   Max movement  BID  Ave movement ')
20 format('# ---- -------------- -------------- -------------- --- --------------')

end subroutine write_stmout_header

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_stmout_results

    use pmf_system_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( .not. ((print_freq .gt. 0) .and. (mod(istep,print_freq) .eq. 0)) ) return

    write(PMF_OUT,10) istep,old_path_len,path_len_change, &
                      MaxMovement, MaxMovementBead, AveMovement

    return

    10 format(I6,1X,E14.7,1X,E14.7,1X,E14.7,1X,I3,1X,E14.7)

end subroutine write_stmout_results

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

! energy values are centered in bins - thus 0.5d0*cvs(i)%width
! value and derivatives are calculated by n-linear interpolation

subroutine fes_value_and_deriv(x,e,d)

    use pmf_system_dat

    implicit none
    real(PMFDP)        :: x(:)
    real(PMFDP)        :: e
    real(PMFDP)        :: d(:)
    ! -----------------------------------------------
    integer            :: i,j,t
    integer            :: ib(ncvs)
    integer            :: it(ncvs)
    real(PMFDP)        :: tx(ncvs)
    real(PMFDP)        :: ty,fac
    real(PMFDP)        :: td(ncvs)
    ! --------------------------------------------------------------------------

    e = 0.0d0
    d = 0.0d0

    ! determine position of point
    do i=1,ncvs
        ib(i) = floor((x(i) - (cvs(i)%min_value + 0.5d0*cvs(i)%width)) / cvs(i)%width)
        if( (ib(i) .lt. 0) .or. (ib(i) .ge. cvs(i)%nbins) ) then
            ! out of boundary
            return
        end if
        tx(i) = (x(i) - ib(i)*cvs(i)%width - (cvs(i)%min_value + 0.5d0*cvs(i)%width)) / cvs(i)%width
    end do

    ! interpolation cell
    do i=0,2**ncvs-1
        t = 1
        fac = 1.0d0
        do j=1,ncvs
            it(j) = mod((i / t),2)
            t = t * 2
            if( it(j) .eq. 0 ) then
                fac = fac * (1.0d0 - tx(j))
            else
                fac = fac * tx(j)
            end if
        end do
        it(:) = ib(:) + it(:)
        call get_dfes_value(it,ty,td)
        e = e + ty*fac
        d = d + td*fac
    end do

end subroutine fes_value_and_deriv

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine get_dfes_value(ib,e,d)

    use pmf_system_dat

    implicit none
    integer            :: ib(ncvs)
    real(PMFDP)        :: e
    real(PMFDP)        :: d(ncvs)
    ! -----------------------------------------------
    integer            :: i,findex
    ! --------------------------------------------------------------------------

    findex = 0
    do i=1,ncvs
        findex = findex*cvs(i)%nbins + ib(i)
    end do

    findex = findex + 1

    if( (findex .lt. 1) .or. (findex .gt. tot_nbins) ) then
        e = 0.0d0
        d = 0.0d0
        return
    end if

    e = fes(findex)
    d = dfes(:,findex)
    return

end subroutine get_dfes_value

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmf_system
