!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module mon_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mon_init_method
!===============================================================================

subroutine mon_init_method

    use mon_output

    implicit none
    ! --------------------------------------------------------------------------

    call mon_init_print_header
    call mon_output_open
    call mon_output_write_header

end subroutine mon_init_method

!===============================================================================
! Subroutine:  mon_init_dat
!===============================================================================

subroutine mon_init_dat

    use mon_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode              = 0         ! 0 - disable MON, 1 - enabled MON
    fsample            = 500       ! output sample pariod in steps

    NumOfMONItems      = 0         ! number of monitored CVs

end subroutine mon_init_dat

!===============================================================================
! Subroutine:  mon_init_print_header
!===============================================================================

subroutine mon_init_print_header

    use mon_dat
    use pmf_dat
    use pmf_utils
    use pmf_cvs

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' ******************* COLLECTIVE VARIABLE MONITORING METHOD ******************** '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' MONITORING Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' MONITORING mode (fmode)                 : ', fmode
    write(PMF_OUT,130)  ' Number of collective variables          : ', NumOfMONItems
    write(PMF_OUT,125)  ' CV definition file (fmondef)            : ', trim(fmondef)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fmonout)                   : ', trim(fmonout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of MONITORING collective variables'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfMONItems
        write(PMF_OUT,140) i
        write(PMF_OUT,145) trim(MONCVList(i)%cv%name)
        write(PMF_OUT,146) trim(MONCVList(i)%cv%ctype)
        write(PMF_OUT,150) MONCVList(i)%cv%get_rvalue(CVContext%CVsValues(MONCVList(i)%cvindx)), &
                           trim(MONCVList(i)%cv%get_ulabel())
        write(PMF_OUT,120)
    enddo

    write(PMF_OUT,120)  '================================================================================'

    return

120 format(A)
125 format(A,A)
130 format(A,I6)
135 format(A,E12.5)
140 format(' == Collective variable #',I2.2)
145 format('    ** Name          : ',a)
146 format('    ** Type          : ',a)
150 format('    ** Current value : ',E16.7,' [',A,']')

end subroutine mon_init_print_header

!===============================================================================

end module mon_init
