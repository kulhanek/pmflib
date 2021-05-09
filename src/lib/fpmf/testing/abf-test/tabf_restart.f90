!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module tabf_restart

implicit none
contains

!===============================================================================
! Subroutine:  tabf_restart_update
!===============================================================================

subroutine tabf_restart_update

    use pmf_dat
    use pmf_utils
    use tabf_accumulator
    use tabf_dat

    implicit none
    !---------------------------------------------------------------------------

    if( frstupdate .le. 0 ) return ! trajectory is written only of frstupdate > 0

    if( mod(fstep,frstupdate) .ne. 0 ) return

    call pmf_utils_open(TABF_RST,ftabfrst,'U')
    call tabf_accumulator_write(TABF_RST)
    close(TABF_RST)

    write(TABF_OUT,10) fstep, insidesamples, outsidesamples

    return

 10 format('# [ACCU] Total steps     = ',I12,' Inside samples  = ',I12,' Outside samples = ',I12  )

end subroutine tabf_restart_update

!===============================================================================
! Subroutine:  tabf_restart_write
!===============================================================================

subroutine tabf_restart_write

    use pmf_dat
    use pmf_utils
    use tabf_accumulator

    implicit none
    !---------------------------------------------------------------------------

    call pmf_utils_open(TABF_RST,ftabfrst,'U')

    call tabf_accumulator_write(TABF_RST)

    close(TABF_RST)

    return

end subroutine tabf_restart_write

!===============================================================================

end module tabf_restart

