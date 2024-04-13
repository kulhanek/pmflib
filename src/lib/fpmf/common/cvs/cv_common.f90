!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_common

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Function:  cv_common_find_cv
!===============================================================================

integer function cv_common_find_cv(cv_name)

    use pmf_dat
    use pmf_utils

    implicit none
    character(*)    :: cv_name
    ! -----------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    cv_common_find_cv = 0

    do i=1,NumOfCVs
        if( trim(cv_name) .eq. trim(CVList(i)%cv%name) ) then
            cv_common_find_cv = i
            return
        end if
    end do

    call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to find CV with name: '''//trim(cv_name)//'''!')

end function cv_common_find_cv

!===============================================================================
! Subroutine:  cv_common_read_name
!===============================================================================

subroutine cv_common_read_name(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    class(CVType)                       :: cv_item
    ! -----------------------------------------------
    integer                             :: i
    character(len=2)                    :: ws
    ! --------------------------------------------------------------------------

    ! load CV name
    if( .not. prmfile_get_string_by_key(prm_fin,'name',cv_item%name) ) then
        call pmf_utils_exit(PMF_OUT,1,'name is not specified!')
    end if

    write(PMF_OUT,100) trim(cv_item%name)

    ws(1:1)   = char(9)    ! tabulator
    ws(2:2)   = ' '        ! space

    if( scan(trim(cv_item%name),ws) .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'CV name contains white characters, which is not allowed!')
    end if

    ! check if name is unique
    do i=1,cv_item%idx-1
        if( trim(cv_item%name) .eq. trim(CVList(i)%cv%name) ) then
            call pmf_utils_exit(PMF_OUT,1,'CV name collision detected!')
        end if
    end do

    return

100 format('   Collective variable name : ''',a,'''')

end subroutine cv_common_read_name

!===============================================================================
! Subroutine:  cv_common_read_groups
!===============================================================================

subroutine cv_common_read_groups(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    class(CVType)                       :: cv_item
    ! -----------------------------------------------
    integer                             :: i
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    ! init groups
    call cv_common_init_groups(cv_item,prm_fin)

    ! read groups
    do i=1,cv_item%ngrps
        call cv_common_read_group(cv_item,prm_fin,i)
    end do

    return

end subroutine cv_common_read_groups

!===============================================================================
! Subroutine:  cv_common_init_groups
!===============================================================================

subroutine cv_common_init_groups(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    class(CVType)                       :: cv_item
    ! -----------------------------------------------
    integer                             :: i, alloc_failed
    character(len=PRMFILE_MAX_LINE)     :: mask
    integer,parameter                   :: group_index = 96   ! ascii code of 'a' - 1
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    ! allocate grps array
    allocate(cv_item%grps(cv_item%ngrps),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
                       'Unable to allocate memory for CV ('//trim(cv_item%name)//')!')
    end if

    ! determine number of atoms in all groups
    cv_item%natoms = 0
    do i=1,cv_item%ngrps
        if( .not. prmfile_get_string_by_key(prm_fin,'group_'//achar(group_index+i),mask) ) then
            call pmf_utils_exit(PMF_OUT,1, &
                                ''//'group_'//achar(group_index+i)//' mask is not defined!')
        end if
        cv_item%grps(i) = pmf_mask_natoms_in_mask(trim(mask))
        if( cv_item%grps(i) .eq. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                                ''//'group_'//achar(group_index+i)//' (' // trim(mask) // ') does not contain any atom!')
        end if
    end do

    ! update grp items
    cv_item%natoms = cv_item%grps(1)
    do i=2,cv_item%ngrps
        cv_item%natoms = cv_item%natoms + cv_item%grps(i)
        cv_item%grps(i) = cv_item%grps(i) + cv_item%grps(i-1)
    end do

    ! allocate arrays
    allocate(cv_item%rindexes(cv_item%natoms),cv_item%lindexes(cv_item%natoms),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
                       'Unable to allocate memory for CV ('//trim(cv_item%name)//')!')
    end if

    return

end subroutine cv_common_init_groups

!===============================================================================
! Subroutine:  cv_common_init_groups_I
!===============================================================================

subroutine cv_common_init_groups_I(cv_item)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    class(CVType)                       :: cv_item
    ! -----------------------------------------------
    integer                             :: alloc_failed
    integer,parameter                   :: group_index = 96   ! ascii code of 'a' - 1
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    ! allocate grps array
    allocate(cv_item%grps(cv_item%ngrps),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
                       'Unable to allocate memory for CV ('//trim(cv_item%name)//')!')
    end if

end subroutine cv_common_init_groups_I

!===============================================================================
! Subroutine:  cv_common_init_groups_II
!===============================================================================

subroutine cv_common_init_groups_II(cv_item,prm_fin,groupid,groupname)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    class(CVType)                       :: cv_item
    integer                             :: groupid
    character(len=PRMFILE_MAX_LINE)     :: groupname
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_LINE)     :: mask
    integer,parameter                   :: group_index = 96   ! ascii code of 'a' - 1
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    if( (groupid .le. 0) .or. (groupid .gt. cv_item%ngrps) ) then
        call pmf_utils_exit(PMF_OUT,1,'Group index out-of-range for CV ('//trim(cv_item%name)//')!')
    end if

    ! get mask
    if( .not. prmfile_get_string_by_key(prm_fin,trim(groupname),mask) ) then
        call pmf_utils_exit(PMF_OUT,1, &
                        '('//trim(groupname)//') group mask is not defined!')
    end if

    cv_item%grps(groupid) = pmf_mask_natoms_in_mask(trim(mask))
    if( cv_item%grps(groupid) .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
                            ''//'group_'//achar(group_index+groupid)//' (' // trim(mask) // ') does not contain any atom!')
    end if

    return

end subroutine cv_common_init_groups_II

!===============================================================================
! Subroutine:  cv_common_init_groups_II_bymask
!===============================================================================

subroutine cv_common_init_groups_II_bymask(cv_item,groupid,mask)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    class(CVType)                       :: cv_item
    integer                             :: groupid
    character(len=PRMFILE_MAX_LINE)     :: mask
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    if( (groupid .le. 0) .or. (groupid .gt. cv_item%ngrps) ) then
        call pmf_utils_exit(PMF_OUT,1,'Group index out-of-range for CV ('//trim(cv_item%name)//')!')
    end if

    cv_item%grps(groupid) = pmf_mask_natoms_in_mask(trim(mask))
    if( cv_item%grps(groupid) .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, 'Mask (' // trim(mask) // ') does not contain any atom!')
    end if

    return

end subroutine cv_common_init_groups_II_bymask

!===============================================================================
! Subroutine:  cv_common_init_groups_III
!===============================================================================

subroutine cv_common_init_groups_III(cv_item)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    class(CVType)                       :: cv_item
    ! -----------------------------------------------
    integer                             :: i, alloc_failed
    integer,parameter                   :: group_index = 96   ! ascii code of 'a' - 1
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    do i=1,cv_item%ngrps
        if( cv_item%grps(i) .eq. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                                ''//'group_'//achar(group_index+i)//' does not contain any atom!')
        end if
    end do

    ! update grp items
    cv_item%natoms = cv_item%grps(1)
    do i=2,cv_item%ngrps
        cv_item%natoms = cv_item%natoms + cv_item%grps(i)
        cv_item%grps(i) = cv_item%grps(i) + cv_item%grps(i-1)
    end do

    ! allocate arrays
    allocate(cv_item%rindexes(cv_item%natoms),cv_item%lindexes(cv_item%natoms),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
                       'Unable to allocate memory for CV ('//trim(cv_item%name)//')!')
    end if

    return

end subroutine cv_common_init_groups_III

!===============================================================================
! Subroutine:  cv_common_read_group
!===============================================================================

subroutine cv_common_read_group(cv_item,prm_fin,groupid)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    class(CVType)                       :: cv_item
    integer                             :: groupid
    ! -----------------------------------------------
    integer                             :: i, k, l, sel
    character(len=PRMFILE_MAX_LINE)     :: mask
    integer,parameter                   :: group_index = 96   ! ascii code of 'a' - 1
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    if( (groupid .le. 0) .or. (groupid .gt. cv_item%ngrps) ) then
        call pmf_utils_exit(PMF_OUT,1,'Group index out-of-range for CV ('//trim(cv_item%name)//')!')
    end if

    ! read group mask
    if( .not. prmfile_get_string_by_key(prm_fin,'group_'//achar(group_index+groupid),mask) ) then
        call pmf_utils_exit(PMF_OUT,1, &
                        ''//'group_'//achar(group_index+groupid)//' mask is not specified!')
    end if

    !set mask
    write(PMF_OUT,100) achar(96+groupid),trim(mask)
    call pmf_mask_set_mask(trim(mask))

    ! populate cv rindexes
    k = 1
    l = 0
    if( groupid .gt. 1 ) k = cv_item%grps(groupid-1) + 1
    do i=1,fnatoms
        call pmf_mask_is_atom_selected(i,sel)
        if( sel .eq. 1 ) then
            cv_item%rindexes(k) = i
            k = k + 1
            l = l + 1
        end if
    end do

    write(PMF_OUT,110) l
    ! print mask
    if( fprint_masks ) then
        call pmf_mask_print();
    end if

    return

100 format('   ** group_',A1,' mask       = ',A)
110 format('      >> Number of atoms : ',I6)

end subroutine cv_common_read_group

!===============================================================================
! Subroutine:  cv_common_read_group_by_name
!===============================================================================

subroutine cv_common_read_group_by_name(cv_item,prm_fin,groupid,groupname)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    class(CVType)                       :: cv_item
    integer                             :: groupid
    character(len=PRMFILE_MAX_LINE)     :: groupname
    ! -----------------------------------------------
    integer                             :: i, k, l, sel
    character(len=PRMFILE_MAX_LINE)     :: mask
    character(len=18)                   :: string
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    if( (groupid .le. 0) .or. (groupid .gt. cv_item%ngrps) ) then
        call pmf_utils_exit(PMF_OUT,1,'Group index out-of-range for CV ('//trim(cv_item%name)//')!')
    end if

    ! get mask
    if( .not. prmfile_get_string_by_key(prm_fin,trim(groupname),mask) ) then
        call pmf_utils_exit(PMF_OUT,1, &
                        '('//trim(groupname)//') group mask is not defined!')
    end if

    !set mask
    string = adjustl(trim(groupname)//' mask')
    write(PMF_OUT,100) string,trim(mask)
    call pmf_mask_set_mask(trim(mask))

    ! populate cv rindexes
    k = 1
    l = 0
    if( groupid .gt. 1 ) k = cv_item%grps(groupid-1) + 1
    do i=1,fnatoms
        call pmf_mask_is_atom_selected(i,sel)
        if( sel .eq. 1 ) then
            cv_item%rindexes(k) = i
            k = k + 1
            l = l + 1
        end if
    end do

    write(PMF_OUT,110) l
    ! print mask
    if( fprint_masks ) then
        call pmf_mask_print();
    end if

    return

100 format('   ** ',A18,' = ',A)
110 format('      >> Number of atoms : ',I6)

end subroutine cv_common_read_group_by_name

!===============================================================================
! Subroutine:  cv_common_read_group_by_mask
!===============================================================================

subroutine cv_common_read_group_by_mask(cv_item,groupid,title,mask)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_mask

    implicit none
    class(CVType)                       :: cv_item
    integer                             :: groupid
    character(len=PRMFILE_MAX_LINE)     :: title
    character(len=PRMFILE_MAX_LINE)     :: mask
    ! -----------------------------------------------
    integer                             :: i, k, l, sel
    character(len=18)                   :: string
    ! --------------------------------------------------------------------------

    if( cv_item%ngrps .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No groups defined for CV ('//trim(cv_item%name)//')!')
    end if

    if( (groupid .le. 0) .or. (groupid .gt. cv_item%ngrps) ) then
        call pmf_utils_exit(PMF_OUT,1,'Group index out-of-range for CV ('//trim(cv_item%name)//')!')
    end if

    !set mask
    string = adjustl(trim(title)//' mask')
    write(PMF_OUT,100) string,trim(mask)
    call pmf_mask_set_mask(trim(mask))

    ! populate cv rindexes
    k = 1
    l = 0
    if( groupid .gt. 1 ) k = cv_item%grps(groupid-1) + 1
    do i=1,fnatoms
        call pmf_mask_is_atom_selected(i,sel)
        if( sel .eq. 1 ) then
            cv_item%rindexes(k) = i
            k = k + 1
            l = l + 1
        end if
    end do

    write(PMF_OUT,110) l
    ! print mask
    if( fprint_masks ) then
        call pmf_mask_print();
    end if

    return

100 format('   ** ',A18,' = ',A)
110 format('      >> Number of atoms : ',I6)

end subroutine cv_common_read_group_by_mask

!===============================================================================
! Subroutine:  cv_common_set_atom_from_mask
!===============================================================================

subroutine cv_common_set_atom_from_mask(mask,atomidx)

    use prmfile
    use pmf_utils
    use pmf_dat
    use pmf_mask

    implicit none
    character(len=PRMFILE_MAX_LINE) :: mask
    integer                         :: atomidx
    ! -----------------------------------------------
    integer                         :: i,k,sel
    ! --------------------------------------------------------------------------

    ! set mask
    call pmf_mask_set_mask(trim(mask))

    ! populate cv rindexes
    k = 0
    do i=1,fnatoms
        call pmf_mask_is_atom_selected(i,sel)
        if( sel .eq. 1 ) then
            atomidx = i
            k = k + 1
        end if
    end do

    if( k .gt. 1 ) then
      call pmf_utils_exit(PMF_OUT,1,'More than one atom in the mask "'//trim(mask)//'"!')
    end if

    if( k .lt. 1 ) then
      call pmf_utils_exit(PMF_OUT,1,'No atom in the mask "'//trim(mask)//'"!')
    end if

    return

end subroutine cv_common_set_atom_from_mask

!===============================================================================
! Subroutine:  cv_common_check_grp_overlap
!===============================================================================

subroutine cv_common_check_grp_overlap(cv_item,groupa,groupb)

    use prmfile
    use pmf_utils
    use pmf_dat
    use pmf_mask

    implicit none
    class(CVType)       :: cv_item
    integer             :: groupa
    integer             :: groupb
    ! -----------------------------------------------
    integer             :: i,j,bi,ei,bj,ej
    integer,parameter   :: group_index = 96   ! ascii code of 'a' - 1
    ! --------------------------------------------------------------------------

    if( (groupa .le. 0) .or. (groupa .gt. cv_item%ngrps) ) then
        call pmf_utils_exit(PMF_OUT,1,'Group index (1) is out-of-range!')
    end if

    if( (groupb .le. 0) .or. (groupb .gt. cv_item%ngrps) ) then
        call pmf_utils_exit(PMF_OUT,1,'Group index (2) is out-of-range!')
    end if

    if( groupa .eq. 1 ) then
        bi = 1
        ei = cv_item%grps(1)
    else
        bi = cv_item%grps(groupa-1)+1
        ei = cv_item%grps(groupa)
    end if

    if( groupb .eq. 1 ) then
        bj = 1
        ej = cv_item%grps(1)
    else
        bj = cv_item%grps(groupb-1)+1
        ej = cv_item%grps(groupb)
    end if

    do i=bi,ei
        do j=bj,ej
            if( cv_item%rindexes(i) .eq. cv_item%rindexes(j) ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     'Atoms of group_'//achar(group_index+groupa)// &
                     ' and group_'//achar(group_index+groupb)//' overlaps!')
            end if
        end do
    end do

    return

end subroutine cv_common_check_grp_overlap

!===============================================================================
! Subroutine:  cv_get_group_com
!===============================================================================

subroutine cv_get_group_com(cv_item,grpid,x,com,totmass)

    use pmf_utils
    use pmf_dat

    implicit none
    class(CVType)   :: cv_item
    integer         :: grpid
    real(PMFDP)     :: x(:,:)
    real(PMFDP)     :: com(3)
    real(PMFDP)     :: totmass
    ! -----------------------------------------------
    integer             :: starti,stopi,i, ai
    real(PMFDP)         :: amass
    integer,parameter   :: group_index = 96   ! ascii code of 'a' - 1
    ! --------------------------------------------------------------------------

    com = 0.0d0
    totmass = 0.0d0

    starti = 1
    if( grpid .gt. 1 ) then
        starti = cv_item%grps(grpid - 1) + 1
    end if
    stopi = cv_item%grps(grpid)

    do i = starti, stopi
        ai = cv_item%lindexes(i)
        amass = mass(ai)
        totmass = totmass + amass
        com(:) = com(:) + amass*x(:,ai)
    end do
    if( totmass .le. 0 ) then
        write(PMF_OUT,*) 'CV name    : '// cv_item%name
        write(PMF_OUT,*) 'Atom group : group_'//achar(group_index+grpid)
        call pmf_utils_exit(PMF_OUT,1,'totmass is zero in cv_get_group_com!')
    end if
    com(:) = com(:) / totmass

end subroutine cv_get_group_com

!===============================================================================
! Function:  cv_get_group_rmass
!===============================================================================

real(PMFDP) function cv_get_group_rmass(cv_item,grpid)

    use pmf_utils
    use pmf_dat

    implicit none
    class(CVType)   :: cv_item
    integer         :: grpid
    ! -----------------------------------------------
    integer         :: starti,stopi,i,ar
    real(PMFDP)     :: amass,totmass
    ! --------------------------------------------------------------------------

    totmass = 0.0d0

    starti = 1
    if( grpid .gt. 1 ) then
        starti = cv_item%grps(grpid - 1) + 1
    end if
    stopi = cv_item%grps(grpid)

    do i = starti, stopi
        ar = cv_item%rindexes(i)
        amass = frmass(ar)
        totmass = totmass + amass
    end do
    cv_get_group_rmass = totmass

end function cv_get_group_rmass

!===============================================================================
! Function:  cv_get_group_natoms
!===============================================================================

integer  function cv_get_group_natoms(cv_item,grpid)

    use pmf_utils
    use pmf_dat

    implicit none
    class(CVType)   :: cv_item
    integer         :: grpid
    ! -----------------------------------------------
    integer         :: starti,stopi
    ! --------------------------------------------------------------------------

    starti = 1
    if( grpid .gt. 1 ) then
        starti = cv_item%grps(grpid - 1) + 1
    end if
    stopi = cv_item%grps(grpid)

    cv_get_group_natoms = stopi - starti + 1

end function cv_get_group_natoms

!===============================================================================
! Subroutine:  cv_get_angle
!===============================================================================

real(PMFDP) function cv_get_angle(v1,v2)

    use pmf_utils
    use pmf_dat

    implicit none
    real(PMFDP)     :: v1(3)
    real(PMFDP)     :: v2(3)
    ! -----------------------------------------------
    real(PMFDP)     :: v1len, v2len
    ! --------------------------------------------------------------------------

    v1len = sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
    if( v1len .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'v1len is zero in cv_get_angle!')
    end if
    v2len = sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
    if( v2len .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'v2len is zero in cv_get_angle!')
    end if

    cv_get_angle = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    cv_get_angle = cv_get_angle / (v1len*v2len)

    if ( cv_get_angle .gt.  1.0 ) then
        cv_get_angle =  1.0
    else if ( cv_get_angle .lt. -1.0 ) then
        cv_get_angle = -1.0
    end if

    cv_get_angle = acos(cv_get_angle)

end function cv_get_angle

!===============================================================================

end module cv_common
