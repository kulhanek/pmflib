!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module gap_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  gap_init_method
!===============================================================================

subroutine gap_init_method

 use gap_output

 implicit none
 ! -----------------------------------------------------------------------------

 call gap_init_print_header
 call gap_init_core
 call gap_output_open
 call gap_output_write_header

end subroutine gap_init_method

!===============================================================================
! Subroutine:  gap_init_dat
!===============================================================================

subroutine gap_init_dat

    use gap_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode              = 0         ! 0 - disable GAP, 1 - enabled GAP
    fsample            = 500       ! output sample pariod in steps
    fplevel            = 0         ! print level

    NumOfGAPCVs      = 0         ! number of CVs used in GAP

    NumOfGAPGroups     = 0         ! number of GAP groups 

end subroutine gap_init_dat

!===============================================================================
! Subroutine:  gap_init_print_header
!===============================================================================

subroutine gap_init_print_header

    use pmf_constants
    use pmf_dat
    use pmf_utils
    use pmf_cvs
    use gap_dat

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' *********************** COLLECTIVE VARIABLE GAP METHOD *********************** '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' GAP Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' GAP mode (fmode)                        : ', fmode
    write(PMF_OUT,130)  ' Number of collective variables          : ', NumOfGAPCVs
    write(PMF_OUT,125)  ' CV definition file (fgapdef)            : ', trim(fgapdef)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fgapout)                   : ', trim(fgapout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,130)  ' Print level (fplevel)                   : ', fplevel
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of collective variables used in GAP'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfGAPCVs
        write(PMF_OUT,140) i
        write(PMF_OUT,145) trim(GAPCVList(i)%cv%name)
        write(PMF_OUT,146) trim(GAPCVList(i)%cv%ctype)
        write(PMF_OUT,150) GAPCVList(i)%cv%get_rvalue(CVContext%CVsValues(GAPCVList(i)%cvindx)), &
                           trim(GAPCVList(i)%cv%get_ulabel())
        write(PMF_OUT,160) trim(GAPCVList(i)%gpfile)
        write(PMF_OUT,170) trim(GAPCVList(i)%groupname)
        write(PMF_OUT,180) GAPCVList(i)%indexid
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
160 format('    ** GP file       : ',a)
170 format('    ** Group name    : ',a)
180 format('    ** Index id      : ',I2)

end subroutine gap_init_print_header

!===============================================================================
! Subroutine:  gap_init_core
!===============================================================================

subroutine gap_init_core

 use pmf_dat
 use pmf_utils
 use gp_dat_mod
 use gp_basic_mod
 use gap_dat

 implicit none
 integer              :: i, j, k
 logical              :: newgroup
 integer              :: alloc_failed
 ! -----------------------------------------------------------------------------

 ! first: check if all specifications are unique
 do i=1,NumOfGAPCVs-1
    do j=i+1,NumOfGAPCVs
       if( GAPCVList(i)%groupname .eq. GAPCVList(j)%groupname ) then
           if( GAPCVList(i)%indexid .eq. GAPCVList(j)%indexid ) then
               write(PMF_OUT,100) trim(GAPCVList(i)%cv%name), trim(GAPCVList(j)%cv%name), &
                                  GAPCVList(i)%groupname, GAPCVList(i)%indexid
               call pmf_utils_exit(PMF_OUT, 1,'[GAP] Same groupname and indexid!')
           end if
       end if
    end do
 end do

 ! second: count the number of GAP groups
 NumOfGAPGroups = 0
 do i=1,NumOfGAPCVs
    newgroup = .true.
    do j=1,i-1
       if( GAPCVList(j)%groupname .eq. GAPCVList(i)%groupname ) then
           newgroup = .false.
           exit
       end if
    end do
    if( newgroup ) then
        NumOfGAPGroups = NumOfGAPGroups + 1
    end if
 end do

 ! third: allocate GAP groups
 allocate(GAPGroupList(NumOfGAPGroups), stat = alloc_failed)
 if( alloc_failed .ne. 0 ) then
     call pmf_utils_exit(PMF_OUT, 1,'[GAP] Unable to allocate memory GAPGroupList!')
 end if

 k=0
 do i=1,NumOfGAPCVs
    newgroup = .true.
    do j=1,k
       if( GAPCVList(i)%groupname .eq. GAPGroupList(j)%groupname ) then
           newgroup = .false.
           exit
       end if
    end do
    if( newgroup ) then
        k = k + 1
        GAPGroupList(k)%groupname = GAPCVList(i)%groupname
        GAPCVList(i)%groupid = k

        call gp_basic_read(GAPGroupList(k)%gp, GAP_INP, trim(GAPCVList(i)%gpfile))
        call gp_basic_complete(GAPGroupList(k)%gp, SE_kernel_r_rr)

        GAPGroupList(k)%nindexes = GAPGroupList(k)%gp%n_dof

        allocate(GAPGroupList(k)%gapcvindx(GAPGroupList(k)%nindexes), stat = alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[GAP] Unable to allocate memory GAPGroupList%gapcvindx!')
        end if
        GAPGroupList(k)%gapcvindx(:) = 0

        allocate(GAPGroupList(k)%values(GAPGroupList(k)%nindexes), stat = alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[GAP] Unable to allocate memory GAPGroupList%values!')
        end if

        allocate(GAPGroupList(k)%forces(GAPGroupList(k)%nindexes), stat = alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[GAP] Unable to allocate memory GAPGroupList%forces!')
        end if

        GAPGroupList(k)%gapcvindx(GAPCVList(i)%indexid) = i
    else
        GAPCVList(i)%groupid = j

        GAPGroupList(j)%gapcvindx(GAPCVList(i)%indexid) = i
    end if
 end do

 ! check if GAP list is complete
 do k=1,NumOfGAPGroups
    do i=1,GAPGroupList(k)%nindexes
       if( GAPGroupList(k)%gapcvindx(i) .eq. 0 ) then
           write(PMF_OUT,200) trim(GAPGroupList(k)%groupname), i 
           call pmf_utils_exit(PMF_OUT, 1,'[GAP] Incomplete GAP!')
       end if
    end do
 end do

 ! get minimum and maximum values
 do k=1,NumOfGAPGroups
    allocate(GAPGroupList(k)%min_values(GAPGroupList(k)%nindexes), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[GAP] Unable to allocate memory GAPGroupList%min_values!')
    end if

    allocate(GAPGroupList(k)%max_values(GAPGroupList(k)%nindexes), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[GAP] Unable to allocate memory GAPGroupList%max_values!')
    end if

    do i=1,GAPGroupList(k)%nindexes
       GAPGroupList(k)%min_values(i) = GAPCVList(GAPCVList(i)%indexid)%min_value
       GAPGroupList(k)%max_values(i) = GAPCVList(GAPCVList(i)%indexid)%max_value
    end do

 end do

100 format(A,' and ',A,' have the same groupname(',A,') and indexid(',I2,')')
200 format('>>> ERROR: missing component in GAP groupname ',A,' in place ',I2,'!')

end subroutine gap_init_core

!===============================================================================

end module gap_init
