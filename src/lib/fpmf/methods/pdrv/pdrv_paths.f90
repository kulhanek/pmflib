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

module pdrv_paths

use pmf_sizes
use pmf_constants

contains

!===============================================================================
! Subroutine:  pdrv_paths_reset_pdrv
!===============================================================================

subroutine pdrv_paths_reset_pdrv(pdrv_item)

    use pdrv_dat

    implicit none
    type(CVTypePDRV)    :: pdrv_item
    ! --------------------------------------------------------------------------

    pdrv_item%mode                  = ' '  ! mode
    pdrv_item%initial_value_set     = .false.
    pdrv_item%set_attach_to_value   = .false.
    pdrv_item%initial_alpha         = 0.0
    pdrv_item%final_alpha           = 0.0

end subroutine pdrv_paths_reset_pdrv

!===============================================================================
! Subroutine:  pdrv_paths_read_pdrv
!===============================================================================

subroutine pdrv_paths_read_pdrv(prm_fin,pdrv_item)

    use prmfile
    use pdrv_dat
    use pmf_dat
    use pmf_unit
    use pmf_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypePDRV)                    :: pdrv_item
    ! -----------------------------------------------
    character(len=1)                    :: buffer
    integer                             :: alloc_stat
    ! --------------------------------------------------------------------------

    if(      prmfile_get_real8_by_key(prm_fin,'attach_to',pdrv_item%final_alpha) &
        .or. prmfile_get_string_by_key(prm_fin,'attach_to',buffer) ) then
        pdrv_item%mode = 'A'
        if( trim(buffer) .eq. '@' ) then
            pdrv_item%set_attach_to_value = .true.
            write(PMF_OUT,101)
        else
            write(PMF_OUT,100) pdrv_item%final_alpha
        end if
        allocate(pdrv_item%ipos(pdrv_item%path%ncvs), &
                pdrv_item%fpos(pdrv_item%path%ncvs), &
                stat = alloc_stat)

        if ( alloc_stat .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[PDRV] Unable to allocate memory for driven path!')
        end if
! --------------------------------------
    else if( prmfile_get_real8_by_key(prm_fin,'increment',pdrv_item%final_alpha) ) then
        if( prmfile_get_real8_by_key(prm_fin,'initial_value',pdrv_item%initial_alpha) ) then
            write(PMF_OUT,10) pdrv_item%initial_alpha
            pdrv_item%initial_value_set = .true.
        else
            write(PMF_OUT,15)
        end if
        pdrv_item%mode = 'I'
        write(PMF_OUT,110) pdrv_item%final_alpha
! --------------------------------------
    else if( prmfile_get_real8_by_key(prm_fin,'change_to',pdrv_item%final_alpha) ) then
        if( prmfile_get_real8_by_key(prm_fin,'initial_value',pdrv_item%initial_alpha) ) then
            write(PMF_OUT,10) pdrv_item%initial_alpha
            pdrv_item%initial_value_set = .true.
        else
            write(PMF_OUT,15)
        end if
        pdrv_item%mode = 'P'
        write(PMF_OUT,120) pdrv_item%final_alpha
! --------------------------------------
    else
        call pmf_utils_exit(PMF_OUT,1,'[PDRV] One of attach_to, increment, change_to keywords must be provided!')
    end if

    return

 10 format('   ** Initial value : ',F7.4)
 15 format('   ** Initial value : - value from input coordinates -')
100 format('   ** Attach to     : ',F7.4)
101 format('   ** Attach to     : - value from input coordinates -')
110 format('   ** Increment by  : ',F7.4)
120 format('   ** Change to     : ',F7.4)

end subroutine pdrv_paths_read_pdrv

!===============================================================================
! Subroutine:   pdrv_paths_pdrv_info
!===============================================================================

subroutine pdrv_paths_pdrv_info(pdrv_item)

    use pdrv_dat

    implicit none
    type(CVTypePDRV)    :: pdrv_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,10) pdrv_item%path%name
    write(PMF_OUT,20) pdrv_item%path%ncvs
    write(PMF_OUT,30) pdrv_item%path%nbeads

    write(PMF_OUT,40) pdrv_item%path%current_alpha
    write(PMF_OUT,50) pmf_paths_get_dist(pdrv_item%path,CVContext,pdrv_item%path%current_alpha)

    select case(pdrv_item%mode)
        case('A')
            write(PMF_OUT,60) 'A','attach CVs to the path'
            if( pdrv_item%set_attach_to_value ) then
                pdrv_item%final_alpha = pdrv_item%path%current_alpha
            end if
            pdrv_item%req_alpha = pdrv_item%final_alpha
            write(PMF_OUT,100) pdrv_item%final_alpha
            call pdrv_paths_pdrv_info_a(pdrv_item)
            pdrv_item%path%rpos = pdrv_item%path%cpos
        case('I')
            write(PMF_OUT,60) 'I','move CVs along the path by an increment'
            if( .not. pdrv_item%initial_value_set ) then
                pdrv_item%initial_alpha = pdrv_item%path%current_alpha
            end if
            call pmf_paths_get_intpoint(pdrv_item%path,pdrv_item%initial_alpha,pdrv_item%path%rpos)
            pdrv_item%final_alpha = pdrv_item%initial_alpha + pdrv_item%final_alpha
            write(PMF_OUT,70) pdrv_item%initial_alpha
            write(PMF_OUT,80) pdrv_item%final_alpha - pdrv_item%initial_alpha
            write(PMF_OUT,90) pdrv_item%final_alpha
            pdrv_item%req_alpha = pdrv_item%initial_alpha
            call pdrv_paths_pdrv_info_ip(pdrv_item)
        case('P')
            write(PMF_OUT,60) 'P','move CVs along the path to the specified point'
            if( .not. pdrv_item%initial_value_set ) then
                pdrv_item%initial_alpha = pdrv_item%path%current_alpha
            end if
            call pmf_paths_get_intpoint(pdrv_item%path,pdrv_item%initial_alpha,pdrv_item%path%rpos)
            write(PMF_OUT,70) pdrv_item%initial_alpha
            write(PMF_OUT,80) pdrv_item%final_alpha - pdrv_item%initial_alpha
            write(PMF_OUT,90) pdrv_item%final_alpha
            pdrv_item%req_alpha = pdrv_item%initial_alpha
            call pdrv_paths_pdrv_info_ip(pdrv_item)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[PDRV] Unsupported driving mode!')
    end select

 10 format('    ** Driven path name  : ',A)
 20 format('    ** Number of CVs     : ',I3)
 30 format('    ** Number of beads   : ',I3)
 40 format('    ** Current alpha     : ',F8.4)
 50 format('       >> CVs Deviation  : ',E14.6)
 60 format('    ** Path driving mode :   ',A,' - ',A)

 70 format('    ** Initial alpha     : ',F8.4)
 80 format('    ** Alpha increment   : ',F8.4)
 90 format('    ** Final alpha       : ',F8.4)
100 format('    ** Attach to alpha   : ',F8.4)

end subroutine pdrv_paths_pdrv_info

!===============================================================================
! Subroutine:   pdrv_paths_pdrv_info_header
!===============================================================================

subroutine pdrv_paths_pdrv_info_header(pdrv_item)

    use pdrv_dat

    implicit none
    type(CVTypePDRV)    :: pdrv_item
    ! --------------------------------------------
    integer             :: i
    ! --------------------------------------------------------------------------

    ! ------------------------
    ! print header
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,100,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        write(PMF_OUT,105,ADVANCE='NO') i
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,200,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        write(PMF_OUT,205,ADVANCE='NO')
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,300,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        write(PMF_OUT,250,ADVANCE='NO') trim(CVList(pdrv_item%path%cvindxs(i))%cv%name)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,400,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        write(PMF_OUT,250,ADVANCE='NO') trim(CVList(pdrv_item%path%cvindxs(i))%cv%ctype)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,500,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        write(PMF_OUT,250,ADVANCE='NO') trim(CVList(pdrv_item%path%cvindxs(i))%cv%get_ulabel())
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,200,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        write(PMF_OUT,205,ADVANCE='NO')
    end do
    write(PMF_OUT,*)

100 format('    #            ')
105 format(1X,4X,'CV#',I2.2,5X)

200 format('    # -----------')
205 format(' --------------')

250 format(1X,A14)

300 format('    # names      ')
400 format('    # types      ')
500 format('    # units      ')

end subroutine pdrv_paths_pdrv_info_header

!===============================================================================
! Subroutine:   pdrv_paths_pdrv_info_a
!===============================================================================

subroutine pdrv_paths_pdrv_info_a(pdrv_item)

    use pdrv_dat

    implicit none
    type(CVTypePDRV)    :: pdrv_item
    ! --------------------------------------------
    integer             :: i
    real(PMFDP)         :: value
    ! --------------------------------------------------------------------------

    call pdrv_paths_pdrv_info_header(pdrv_item)

    do i=1,pdrv_item%path%ncvs
        pdrv_item%ipos(i) = pdrv_item%path%cpos(i)
        pdrv_item%fpos(i) = pmf_paths_get_intcv(pdrv_item%path,i,pdrv_item%final_alpha)
    end do

    ! ------------------------
    write(PMF_OUT,100,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pdrv_item%ipos(i)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,150,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pdrv_item%fpos(i) - pdrv_item%ipos(i)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,200,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pdrv_item%fpos(i)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

100 format('    # Current CVs')
150 format('    # Increment  ')
200 format('    # Final CVs  ')

255 format(1X,E14.6)

end subroutine pdrv_paths_pdrv_info_a

!===============================================================================
! Subroutine:   pdrv_paths_pdrv_info_ip
!===============================================================================

subroutine pdrv_paths_pdrv_info_ip(pdrv_item)

    use pdrv_dat

    implicit none
    type(CVTypePDRV)    :: pdrv_item
    ! --------------------------------------------
    integer             :: i
    real(PMFDP)         :: value
    ! --------------------------------------------------------------------------

    call pdrv_paths_pdrv_info_header(pdrv_item)

    ! ------------------------
    write(PMF_OUT,100,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pdrv_item%path%cpos(i)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,110,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pmf_paths_get_intcv(pdrv_item%path,i,pdrv_item%initial_alpha) - pdrv_item%path%cpos(i)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,120,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pmf_paths_get_intcv(pdrv_item%path,i,pdrv_item%initial_alpha)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,130,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pmf_paths_get_intcv(pdrv_item%path,i,pdrv_item%final_alpha) &
              - pmf_paths_get_intcv(pdrv_item%path,i,pdrv_item%initial_alpha)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

    ! ------------------------
    write(PMF_OUT,140,ADVANCE='NO')
    do i=1,pdrv_item%path%ncvs
        value = pmf_paths_get_intcv(pdrv_item%path,i,pdrv_item%final_alpha)
        write(PMF_OUT,255,ADVANCE='NO') CVList(pdrv_item%path%cvindxs(i))%cv%get_rvalue(value)
    end do
    write(PMF_OUT,*)

100 format('    # Current CVs')
110 format('    # Deviation  ')
120 format('    # Initial CVs')
130 format('    # Increment  ')
140 format('    # Final CVs  ')

255 format(1X,E14.6)

end subroutine pdrv_paths_pdrv_info_ip

!===============================================================================
! Subroutine:   pdrv_paths_increment
!===============================================================================

subroutine pdrv_paths_increments

    use pdrv_dat

    implicit none
    integer             :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfPDRVItems
        call pdrv_paths_pdrv_increment(PDRVCVList(i))
    end do

end subroutine pdrv_paths_increments

!===============================================================================
! Subroutine:   pdrv_paths_pdrv_increment
!===============================================================================

subroutine pdrv_paths_pdrv_increment(pdrv_item)

    use pdrv_dat

    implicit none
    type(CVTypePDRV)    :: pdrv_item
    ! --------------------------------------------
    real(PMFDP)         :: alpha
    ! --------------------------------------------------------------------------

    select case(pdrv_item%mode)
        case('A')
            pdrv_item%req_alpha = pdrv_item%final_alpha
            alpha = real(fstep) / real(fnstlim)
            pdrv_item%path%rpos = pdrv_item%ipos + (pdrv_item%fpos - pdrv_item%ipos)*alpha
        case('I','P')
            alpha = pdrv_item%initial_alpha &
                  + (pdrv_item%final_alpha - pdrv_item%initial_alpha)* real(fstep) / real(fnstlim)
            pdrv_item%req_alpha = alpha
            call pmf_paths_get_intpoint(pdrv_item%path,alpha,pdrv_item%path%rpos)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[PDRV] Unsupported driving mode!')
    end select

end subroutine pdrv_paths_pdrv_increment

!===============================================================================

end module pdrv_paths
