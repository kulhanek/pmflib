!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

! MDIS - minimum distance between pairs of atoms

module cv_mdis

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeMDIS

    real(PMFDP) :: beta

    contains
        procedure :: load_cv        => load_mdis
        procedure :: calculate_cv   => calculate_mdis
end type CVTypeMDIS

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_mdis
!===============================================================================

subroutine load_mdis(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypeMDIS)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    logical                             :: ret
    integer                             :: count, atomid, alloc_failed, sbeg, send
    character(PRMFILE_MAX_LINE)         :: key
    character(PRMFILE_MAX_LINE)         :: mask1
    character(PRMFILE_MAX_LINE)         :: mask2
    character(PRMFILE_MAX_LINE)         :: locvalue
    ! --------------------------------------------------------------------------

    ! simple init and allocation --------------------
    cv_item%ctype         = 'MDIS'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! determine number of pairs ---------------------
    ret = prmfile_first_line(prm_fin)
    count = 0
    do while( prmfile_get_field_on_line(prm_fin,key) .and. ret)
        if( trim(key) .eq. 'pair' )  count = count + 1
        ret = prmfile_next_line(prm_fin)
    end do

    ! initialize groups -----------------------------
    cv_item%ngrps = 1

    ! allocate grps array
    allocate(cv_item%grps(cv_item%ngrps),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory!')
    end if

    cv_item%natoms  = 2*count
    cv_item%grps(1) = cv_item%natoms

    ! allocate arrays
    allocate(cv_item%rindexes(cv_item%natoms),cv_item%lindexes(cv_item%natoms),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory!')
    end if

    ! read distance pairs
    ret = prmfile_first_line(prm_fin)
    count = 1
    atomid = 1
    do while( prmfile_get_field_on_line(prm_fin,key) .and. ret)
        if( trim(key) .eq. 'pair' ) then
            if( .not. prmfile_get_field_on_line(prm_fin,locvalue) ) then
                write(locvalue,*) count
                call pmf_utils_exit(PMF_OUT,1, &
                        'Unable to get the first atom mask for the pair '//trim(locvalue)//'!')
            end if

            ! trim quotation
            sbeg = 1
            send = len(trim(locvalue))

            ! handle quotations
            if( send .ge. 2 ) then
                if( (locvalue(1:1) .eq. '''') .and. (locvalue(send:send) .eq. '''')  ) then
                    sbeg = sbeg + 1
                    send = send - 1
                else if( (locvalue(1:1) .eq. '"') .and. (locvalue(send:send) .eq. '"') ) then
                    sbeg = sbeg + 1
                    send = send - 1
                end if
            end if

            if( sbeg .le. send ) then
                mask1 = locvalue(sbeg:send)
            else
                mask1 = ''
            end if

            if( .not. prmfile_get_field_on_line(prm_fin,locvalue) ) then
                write(locvalue,*) count
                call pmf_utils_exit(PMF_OUT,1, &
                        'Unable to get the second atom mask for the pair '//trim(locvalue)//'!')
            end if

            ! trim quotation
            sbeg = 1
            send = len(trim(locvalue))

            ! handle quotations
            if( send .ge. 2 ) then
                if( (locvalue(1:1) .eq. '''') .and. (locvalue(send:send) .eq. '''')  ) then
                    sbeg = sbeg + 1
                    send = send - 1
                else if( (locvalue(1:1) .eq. '"') .and. (locvalue(send:send) .eq. '"') ) then
                    sbeg = sbeg + 1
                    send = send - 1
                end if
            end if

            if( sbeg .le. send ) then
                mask2 = locvalue(sbeg:send)
            else
                mask2 = ''
            end if

            write(PMF_OUT,20) trim(mask1),trim(mask2)
            call cv_common_set_atom_from_mask(mask1,cv_item%rindexes(atomid))

            atomid = atomid + 1
            call cv_common_set_atom_from_mask(mask2,cv_item%rindexes(atomid))

            atomid = atomid + 1
            count = count + 1
        end if
        ret = prmfile_next_line(prm_fin)
    end do

    ! read beta
    cv_item%beta = 1.0d0
    if( prmfile_get_real8_by_key(prm_fin,'beta',cv_item%beta) ) then
        write(PMF_OUT,10) cv_item%beta
    else
        write(PMF_OUT,15) cv_item%beta
    end if

    20 format('   ** atom pair          : ',A10,' <-> ',A10)
    10 format('   ** beta               : ',F5.2)
    15 format('   ** beta               : ',F5.2,' (default)')

end subroutine load_mdis

!===============================================================================
! Subroutine:  calculate_mdis
!===============================================================================

subroutine calculate_mdis(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils
    use pmf_pbc

    implicit none
    class(CVTypeMDIS)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: m,ai1,ai2
    real(PMFDP)    :: d1(3),d2(3),dd(3),dx2,dx,ex,beta,sc,lsc
    ! --------------------------------------------------------------------------

    ! calculate value

    m = 1
    beta = cv_item%beta
    ex = 0.0d0
    do while(m .lt. cv_item%natoms)
        ai1 = cv_item%lindexes(m)
        d1(:) = x(:,ai1)
        m = m + 1
        ai2 = cv_item%lindexes(m)
        d2(:) = x(:,ai2)
        m = m + 1
        dd(:) = d1(:) - d2(:)

        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dd)
        end if

        dx2 = dd(1)**2 + dd(2)**2 + dd(3)**2
        dx  = sqrt(dx2)
        if( dx .lt. 1.0e-7 ) then
            call pmf_utils_exit(PMF_OUT,1,'Distance is smaller than 1.0e-7 in calculate_mdis!')
        end if
        ex = ex + exp( - dx / beta)
    end do

    ctx%CVsValues(cv_item%idx) = - beta * log(ex)

    ! ------------------------------------------------

    ! calculate gradient
    m = 1
    sc = 1.0d0 / ex
    do while(m .lt. cv_item%natoms)
        ai1 = cv_item%lindexes(m)
        d1(:) = x(:,ai1)
        m = m + 1
        ai2 = cv_item%lindexes(m)
        d2(:) = x(:,ai2)
        m = m + 1

        dd(:) = d1(:) - d2(:)
        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dd)
        end if

        dx2 = dd(1)**2 + dd(2)**2 + dd(3)**2
        dx  = sqrt(dx2)
        lsc = sc * exp(- dx / beta) / dx
        ctx%CVsDrvs(:,ai1,cv_item%idx) = ctx%CVsDrvs(:,ai1,cv_item%idx) + lsc*dd(:)
        ctx%CVsDrvs(:,ai2,cv_item%idx) = ctx%CVsDrvs(:,ai2,cv_item%idx) - lsc*dd(:)
    end do

end subroutine calculate_mdis

!===============================================================================

end module cv_mdis

