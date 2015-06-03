!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_paths

use pmf_sizes
use pmf_constants
use pmf_utils
use pmf_dat

! list of paths ----------------------------------------------------------------

! core definition of paths -----------------------------------------------------
type PathType
    integer                     :: idx              ! index of PATH
    character(PMF_MAX_CV_NAME)  :: name             ! PATH name
    integer                     :: ncvs             ! number of CVs
    integer,pointer             :: cvindxs(:)       ! indexes of CVs
    type(CVPointer),pointer     :: cvs(:)           ! CVs - dimm (ncvs)
    integer                     :: nbeads           ! number of beads
    real(PMFDP),pointer         :: alphas(:)        ! alphas
    real(PMFDP),pointer         :: points(:,:)      ! points along path - dimm (nbeads,ncvs) !!!!
    logical,pointer             :: fixed(:)         ! point statuses - dimm (nbeads)
    real(PMFDP),pointer         :: minvalues(:)     ! minimum allowed CV values - dimm (ncvs)
    real(PMFDP),pointer         :: maxvalues(:)     ! maximum allowed CV values - dimm (ncvs)
    real(PMFDP),pointer         :: maxmoves(:)      ! maximum allowed CV value change - dimm (ncvs)
    ! cubic spline
    real(PMFDP),pointer         :: ypp(:,:)         ! points along path - dimm (nbeads,ncvs) !!!!
    real(PMFDP),pointer         :: spos(:)          ! helper - dimm (ncvs)
    ! runtime data
    logical                     :: driven_mode      ! someone else control this path
    real(PMFDP)                 :: current_alpha    ! current alpha value
    real(PMFDP),pointer         :: cpos(:)          ! current CV positions
    real(PMFDP),pointer         :: rpos(:)          ! required CV positions
end type PathType

type PathPointer
    class(PathType),pointer :: path                 ! path data
end type PathPointer

integer                         :: NumOfPaths       ! number of Paths
type(PathPointer),allocatable   :: PathList(:)      ! input definition of Paths

!-------------------------------------------------------------------------------

! how many points are used to calculate path length
integer,parameter       :: PMF_PATH_SEGMENT_DISCRETIZATION = 10

! what is the precission to find closest alpha 1/PMF_PATH_ALPHA_PRECISION
! in focusing segment
real(PMFDP),parameter   :: PMF_PATH_ALPHA_PRECISION = 20.0d0
! number of focusing
integer,parameter       :: PMF_PATH_ALPHA_NFOCUS = 4

!-------------------------------------------------------------------------------

 contains

!===============================================================================
! Subroutine:   pmf_paths_load_path
! read path from [PATH] section
!===============================================================================

subroutine pmf_paths_load_path(prm_fin,path_item)

    use prmfile
    use pmf_core
    use cv_common

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(PathType)                      :: path_item
    ! --------------------------------------------
    integer                             :: alloc_failed, user_nbeads, i, b
    type(PathType)                      :: incomplete_path_item
    ! -----------------------------------------------------------------------------

    ! flexible   cv1 cv2 cv3 ... cvn
    ! permanent  cv1 cv2 cv3 ... cvn

    write(PMF_OUT,'(/,a)') '=== [PATH] ====================================================================='

    ! load path name
    call pmf_paths_load_name(prm_fin,path_item)

    ! read number of points and CVs
    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'ncvs',path_item%ncvs) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of CVs (ncvs) must be specified!')
    end if
    write(PMF_OUT,10) path_item%ncvs

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbeads',path_item%nbeads) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of beads (nbeads) must be specified!')
    end if
    write(PMF_OUT,20) path_item%nbeads

    ! ========================
    ! determine number of provided beads
    user_nbeads = pmf_paths_get_nbeads(prm_fin)
    write(PMF_OUT,30) user_nbeads

    if( user_nbeads .gt. path_item%nbeads ) then
        call pmf_utils_exit(PMF_OUT,1,'More bead specifications than number provided via ''nbeads'' keyword!')
    end if

    if( user_nbeads .eq. path_item%nbeads ) then
        write(PMF_OUT,40)
    else
        write(PMF_OUT,50)
    end if

    ! allocate path data
    allocate(path_item%cvindxs(path_item%ncvs), &
             path_item%cvs(path_item%ncvs), &
             path_item%minvalues(path_item%ncvs), &
             path_item%maxvalues(path_item%ncvs), &
             path_item%maxmoves(path_item%ncvs), &
             path_item%points(path_item%nbeads,path_item%ncvs), &
             path_item%fixed(path_item%nbeads), &
             path_item%alphas(path_item%nbeads), &
             path_item%ypp(path_item%nbeads,path_item%ncvs), &
             path_item%spos(path_item%ncvs), &
             path_item%cpos(path_item%ncvs), &
             path_item%rpos(path_item%ncvs), &
             stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for the path!')
    end if

    path_item%fixed(:) = .false.
    path_item%driven_mode = .false.

    ! -----------------------------------------------
    ! write header
    call pmf_paths_print_header(path_item)

    ! -----------------------------------------------
    ! assign CVs
    call pmf_paths_load_cvs(prm_fin,path_item)

    ! -----------------------------------------------

    if( user_nbeads .eq. path_item%nbeads ) then
        call pmf_paths_load_beads(prm_fin,path_item)
    else
        ! allocate intermediate path data
        incomplete_path_item%nbeads = user_nbeads
        incomplete_path_item%ncvs = path_item%ncvs
        allocate(incomplete_path_item%cvindxs(path_item%ncvs), &
                 incomplete_path_item%cvs(path_item%ncvs), &
                 incomplete_path_item%minvalues(path_item%ncvs), &
                 incomplete_path_item%maxvalues(path_item%ncvs), &
                 incomplete_path_item%maxmoves(path_item%ncvs), &
                 incomplete_path_item%points(user_nbeads,path_item%ncvs), &
                 incomplete_path_item%fixed(user_nbeads), &
                 incomplete_path_item%alphas(user_nbeads), &
                 incomplete_path_item%ypp(user_nbeads,path_item%ncvs), &
                 incomplete_path_item%spos(path_item%ncvs), &
                 stat=alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for the intermediate path!')
        end if
        incomplete_path_item%cvindxs(:) = path_item%cvindxs(:)
        incomplete_path_item%cvs(:) = path_item%cvs(:)
        incomplete_path_item%minvalues(:) = path_item%minvalues(:)
        incomplete_path_item%maxvalues(:) = path_item%maxvalues(:)
        incomplete_path_item%maxmoves(:) = path_item%maxmoves(:)
        incomplete_path_item%fixed(:) = .false.

        ! load user data
        call pmf_paths_load_beads(prm_fin,incomplete_path_item)

        ! optimize path
        call pmf_paths_optimize_alphas(incomplete_path_item)

        do b=2,incomplete_path_item%nbeads-1
            if( incomplete_path_item%fixed(b) ) then
                call pmf_utils_exit(PMF_OUT,1,'Permanent beads can be only at path ends!')
            end if
        end do

        ! construct final path
        path_item%alphas(1) = 0.0
        path_item%fixed(1) = incomplete_path_item%fixed(1)
        do b=2,path_item%nbeads-1
            path_item%alphas(b) = real(b-1) / real(path_item%nbeads-1)
        end do
        path_item%alphas(path_item%nbeads) = 1.0
        path_item%fixed(path_item%nbeads) = incomplete_path_item%fixed(incomplete_path_item%nbeads)

        do b=1,path_item%nbeads
            do i=1,path_item%ncvs
                path_item%points(b,i) = pmf_paths_get_intcv(incomplete_path_item,i,path_item%alphas(b))
            end do
        end do

        ! relase resources
        deallocate(incomplete_path_item%cvindxs, &
                 incomplete_path_item%cvs, &
                 incomplete_path_item%minvalues, &
                 incomplete_path_item%maxvalues, &
                 incomplete_path_item%maxmoves, &
                 incomplete_path_item%points, &
                 incomplete_path_item%fixed, &
                 incomplete_path_item%alphas, &
                 incomplete_path_item%ypp, &
                 incomplete_path_item%spos )
    end if

    ! optimize path
    call pmf_paths_optimize_alphas(path_item)

    write(PMF_OUT,60)
    call pmf_paths_print_path(path_item)

 return

 10 format(  'Number of CVs            = ',I3)
 20 format(  'Number of beads          = ',I3)
 30 format(  'Number of provided beads = ',I3)
 40 format(/,'>>> Complete path provided ...')
 50 format(/,'>>> Incomplete path provided, it will be reconstructed ...')
 60 format(/,'>>> Final path with optimized alphas ...')

end subroutine pmf_paths_load_path

!===============================================================================
! Subroutine:   pmf_paths_load_name
!===============================================================================

subroutine pmf_paths_load_name(prm_fin,path_item)

    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(PathType)                      :: path_item
    ! --------------------------------------------
    integer                             :: i
    character(len=2)                    :: ws
    ! --------------------------------------------------------------------------

    ! load path name
    if( .not. prmfile_get_string_by_key(prm_fin,'name',path_item%name) ) then
        call pmf_utils_exit(PMF_OUT,1,'The name of path (name) is not specified!')
    end if

    write(PMF_OUT,10) trim(path_item%name)
    write(PMF_OUT,20) path_item%idx

    ws(1:1)   = char(9)    ! tabulator
    ws(2:2)   = ' '        ! space

    if( scan(trim(path_item%name),ws) .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'The name of path contains white characters, which is not allowed!')
    end if

    ! check if name is unique
    do i=1,path_item%idx-1
        if( trim(path_item%name) .eq. trim(PathList(i)%path%name) ) then
            call pmf_utils_exit(PMF_OUT,1,'Path name collision detected!')
        end if
    end do

    ! check collisions with CV names
    do i=1,NumOfCVs
        if( trim(path_item%name) .eq. trim(CVList(i)%cv%name) ) then
            call pmf_utils_exit(PMF_OUT,1,'Path name/CV name collision detected!')
        end if
    end do

    return

 10 format(  'Path name                = ',A)
 20 format(  'Path ID                  = ',I3)

end subroutine pmf_paths_load_name

!===============================================================================
! Subroutine:   pmf_paths_load_cvs
!===============================================================================

subroutine pmf_paths_load_cvs(prm_fin,path_item)

    use prmfile
    use pmf_core
    use cv_common

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(PathType)                      :: path_item
    ! --------------------------------------------
    character(PRMFILE_MAX_LINE)             :: text
    character(15)                           :: code
    integer                                 :: alloc_failed, i
    character(PMF_MAX_CV_NAME),allocatable  :: buff(:)
    logical                                 :: res
    logical                                 :: found
    ! -----------------------------------------------------------------------------

    allocate(buff(path_item%ncvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for the path!')
    end if

    ! -----------------------------------------------
    ! assign CVs
    res = prmfile_first_line(prm_fin)
    found = .false.
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( trim(code) .eq. 'names' ) then
            read(text,*,err=10,end=10) code,(buff(i),i=1,path_item%ncvs)
            found = .true.
            exit
        end if
    end do

    if( found .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'The keyword ''names'' is mandatory!')
    end if

    write(PMF_OUT,100,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,105,ADVANCE='NO') trim(buff(i))
    end do
    write(PMF_OUT,*)

    do i=1,path_item%ncvs
        path_item%cvindxs(i) = cv_common_find_cv(buff(i))
        path_item%cvs(i)%cv  => CVList(path_item%cvindxs(i))%cv
        if( CVList(path_item%cvindxs(i))%cv%pathidx .ne. 0 ) then
            write(code,'(I3)') CVList(path_item%cvindxs(i))%cv%pathidx
            call pmf_utils_exit(PMF_OUT,1,'Collective variable ''' // &
                                trim(CVList(path_item%cvindxs(i))%cv%name) // ''' is already member of path number ' // &
                                trim(code) // '!')
        end if
        CVList(path_item%cvindxs(i))%cv%pathidx = path_item%idx
    end do

    ! -----------------------------------------------
    ! check CVs types
    res = prmfile_first_line(prm_fin)
    found = .false.
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( trim(code) .eq. 'types' ) then
            read(text,*,err=20,end=20) code,(buff(i),i=1,path_item%ncvs)
            found = .true.
            exit
        end if
    end do

    if( found .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'The keyword ''types'' is mandatory!')
    end if

    write(PMF_OUT,200,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,205,ADVANCE='NO') trim(buff(i))
    end do
    write(PMF_OUT,*)

    do i=1,path_item%ncvs
        if( trim(buff(i)) .ne. trim(CVList(path_item%cvindxs(i))%cv%ctype) ) then
            write(code,'(I3)') i
            call pmf_utils_exit(PMF_OUT,1,'CV types do not match for CV number ''' // trim(code) // '''!')
        end if
    end do

    ! -----------------------------------------------
    ! print CV units
    write(PMF_OUT,250,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,255,ADVANCE='NO') trim(CVList(path_item%cvindxs(i))%cv%get_ulabel())
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    ! load min values
    res = prmfile_first_line(prm_fin)
    found = .false.
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( trim(code) .eq. 'min' ) then
            read(text,*,err=20,end=20) code,(path_item%minvalues(i),i=1,path_item%ncvs)
            found = .true.
            exit
        end if
    end do

    if( found .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'The keyword ''min'' is mandatory!')
    end if

    write(PMF_OUT,300,ADVANCE='NO')
    do i=1,path_item%ncvs
        call CVList(path_item%cvindxs(i))%cv%conv_to_ivalue(path_item%minvalues(i))
        write(PMF_OUT,305,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%minvalues(i))
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    ! load max values
    res = prmfile_first_line(prm_fin)
    found = .false.
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( trim(code) .eq. 'max' ) then
            read(text,*,err=20,end=20) code,(path_item%maxvalues(i),i=1,path_item%ncvs)
            found = .true.
            exit
        end if
    end do

    if( found .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'The keyword ''max'' is mandatory!')
    end if

    write(PMF_OUT,400,ADVANCE='NO')
    do i=1,path_item%ncvs
        call CVList(path_item%cvindxs(i))%cv%conv_to_ivalue(path_item%maxvalues(i))
        write(PMF_OUT,405,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%maxvalues(i))
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    ! load maxmov values
    path_item%maxmoves(:) = 0.0d0
    res = prmfile_first_line(prm_fin)
    found = .false.
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( trim(code) .eq. 'maxmov' ) then
            read(text,*,err=20,end=20) code,(path_item%maxmoves(i),i=1,path_item%ncvs)
            found = .true.
            exit
        end if
    end do

    ! maxmov is optional

    write(PMF_OUT,500,ADVANCE='NO')
    do i=1,path_item%ncvs
        call CVList(path_item%cvindxs(i))%cv%conv_to_ivalue(path_item%maxmoves(i))
        write(PMF_OUT,505,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%maxmoves(i))
    end do
    write(PMF_OUT,*)

    write(PMF_OUT,600,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,605,ADVANCE='NO')
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    ! clean up
    deallocate(buff)
    return

10 call pmf_utils_exit(PMF_OUT,1,'The keyword ''names'' is not provided correctly!')
20 call pmf_utils_exit(PMF_OUT,1,'The keyword ''types'' is not provided correctly!')
30 call pmf_utils_exit(PMF_OUT,1,'The keyword ''min'' is not provided correctly!')
40 call pmf_utils_exit(PMF_OUT,1,'The keyword ''max'' is not provided correctly!')
50 call pmf_utils_exit(PMF_OUT,1,'The keyword ''maxmov'' is not provided correctly!')

100 format('#      names        ')
105 format(1X,A14)

200 format('#      types        ')
205 format(1X,A14)

250 format('#      units        ')
255 format(1X,A14)

300 format('#      min          ')
305 format(1X,E14.6)

400 format('#      max          ')
405 format(1X,E14.6)

500 format('#      maxmov       ')
505 format(1X,E14.6)

600 format('# ---- ------ ------')
605 format(' --------------')

end subroutine pmf_paths_load_cvs

!===============================================================================
! Function:   pmf_paths_load_fake_cvs
!===============================================================================

subroutine pmf_paths_load_fake_cvs(prm_fin)

    use prmfile
    use pmf_core
    use pmf_alloc_cv

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------
    character(PRMFILE_MAX_LINE)             :: text
    character(15)                           :: code
    integer                                 :: ncvs, alloc_failed, i
    character(PMF_MAX_CV_NAME),allocatable  :: types(:)
    character(PMF_MAX_CV_NAME),allocatable  :: names(:)
    logical                                 :: res
    ! -----------------------------------------------------------------------------

    if( .not. prmfile_get_integer_by_key(prm_fin,'ncvs',ncvs) ) then
        return
    end if

    allocate(types(ncvs), names(ncvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for the path!')
    end if

    types(:) = ''
    names(:) = ''

    ! -----------------------------------------------
    ! load types and names
    res = prmfile_first_line(prm_fin)
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( trim(code) .eq. 'types' ) then
            read(text,*,err=20,end=20) code,(types(i),i=1,ncvs)
        end if
        if( trim(code) .eq. 'names' ) then
            read(text,*,err=10,end=10) code,(names(i),i=1,ncvs)
        end if
    end do

    ! generate CVs
    do i=1,ncvs
        NumOfFakeCVs = NumOfFakeCVs + 1
        if( NumOfFakeCVs .gt. NumOfCVs ) then
            call pmf_utils_exit(PMF_OUT,1,'Inconsistency in number of fake CVs and maximum number of CVs!')
        end if

        write(PMF_OUT,*)
        write(PMF_OUT,130) NumOfFakeCVs,trim(types(i))

        call pmf_alloc_cv_allocate(types(i),CVList(NumOfFakeCVs)%cv)

        ! init and load CV data
        call CVList(NumOfFakeCVs)%cv%reset_cv()
        CVList(NumOfFakeCVs)%cv%idx   = NumOfFakeCVs
        CVList(NumOfFakeCVs)%cv%name  = names(i)
        CVList(NumOfFakeCVs)%cv%ctype = types(i)

        write(PMF_OUT,140) trim(CVList(NumOfFakeCVs)%cv%name)
    end do

    ! -----------------------------------------------
    ! clean up
    deallocate(types,names)
    return

10 call pmf_utils_exit(PMF_OUT,1,'The keyword ''names'' is not provided correctly!')
20 call pmf_utils_exit(PMF_OUT,1,'The keyword ''types'' is not provided correctly!')

130 format('== Creating collective variable #',I4.4,' of type "',A,'"')
140 format('   ** Collective variable name : ''',a,'''')
end subroutine pmf_paths_load_fake_cvs

!===============================================================================
! Function:   pmf_paths_get_nbeads
!===============================================================================

real(PMFDP) function pmf_paths_get_nbeads(prm_fin)

    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------
    character(PRMFILE_MAX_LINE)             :: text
    character(15)                           :: code
    logical                                 :: res
    ! -----------------------------------------------------------------------------

    pmf_paths_get_nbeads = 0

    res = prmfile_first_line(prm_fin)
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( trim(code) .eq. 'permanent' ) pmf_paths_get_nbeads = pmf_paths_get_nbeads + 1
        if( trim(code) .eq. 'flexible' )  pmf_paths_get_nbeads = pmf_paths_get_nbeads + 1
    end do

end function pmf_paths_get_nbeads

!===============================================================================
! Subroutine:   pmf_paths_load_beads
!===============================================================================

subroutine pmf_paths_load_beads(prm_fin,path_item)

    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(PathType)                      :: path_item
    ! --------------------------------------------
    character(PRMFILE_MAX_LINE)         :: text
    character(15)                       :: code
    logical                             :: res
    integer                             :: b, i
    ! --------------------------------------------------------------------------

    res = prmfile_first_line(prm_fin)
    b = 1
    do while (prmfile_get_line(prm_fin,text))
        read(text,*) code
        if( (trim(code) .eq. 'permanent') .or. (trim(code) .eq. 'flexible') ) then
            path_item%fixed(b) = trim(code) .eq. 'permanent'
            if(  path_item%fixed(b) ) then
                write(PMF_OUT,25,ADVANCE='NO') b
            else
                write(PMF_OUT,20,ADVANCE='NO') b
            end if
            read(text,*,err=10,end=10) code, (path_item%points(b,i),i=1,path_item%ncvs)
            do i=1,path_item%ncvs
                call CVList(path_item%cvindxs(i))%cv%conv_to_ivalue(path_item%points(b,i))
                write(PMF_OUT,30,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%points(b,i))
            end do
            write(PMF_OUT,*)
            b = b + 1
        end if
    end do

    return

10 call pmf_utils_exit(PMF_OUT,1,'The keyword ''permanent'' or ''flexible'' is not provided correctly!')

20 format(I6,   ' flexible     ')
25 format(I6,   ' permanent    ')
30 format(1X,E14.6)

end subroutine pmf_paths_load_beads

!===============================================================================
! Subroutine:   pmf_paths_write_path
! write path including [PATH] section
!===============================================================================

subroutine pmf_paths_write_path(iounit,path_item)

    use prmfile
    use pmf_core
    use cv_common

    implicit none
    integer         :: iounit
    type(PathType)  :: path_item
    ! --------------------------------------------
    integer         :: i, b
    ! -----------------------------------------------------------------------------

    write(iounit,10)
    write(iounit,20) trim(path_item%name)
    write(iounit,30) path_item%ncvs
    write(iounit,40) path_item%nbeads

    write(iounit,100,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,50,ADVANCE='NO') trim(path_item%cvs(i)%cv%name)
    end do
    write(iounit,*)

    write(iounit,110,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,50,ADVANCE='NO') trim(path_item%cvs(i)%cv%ctype)
    end do
    write(iounit,*)

    write(iounit,120,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,60,ADVANCE='NO') path_item%minvalues(i)
    end do
    write(iounit,*)

    write(iounit,130,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,60,ADVANCE='NO') path_item%maxvalues(i)
    end do
    write(iounit,*)

    write(iounit,140,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,60,ADVANCE='NO') path_item%maxmoves(i)
    end do
    write(iounit,*)

    do b=1,path_item%nbeads
        if( path_item%fixed(b) ) then
            write(iounit,160,ADVANCE='NO')
        else
            write(iounit,150,ADVANCE='NO')
        end if
        do i=1,path_item%ncvs
            write(iounit,60,ADVANCE='NO') path_item%points(b,i)
        end do
        write(iounit,*)
    end do

    write(iounit,*)

 10 format('[PATH]')
 20 format('name     ',A)
 30 format('ncvs     ',I3)
 40 format('nbeads   ',I3)

 50 format(1X,A14)
 60 format(1X,E14.7)

100 format('names    ')
110 format('types    ')
120 format('min      ')
130 format('max      ')
140 format('maxmov   ')
150 format('flexible ')
160 format('permanent')

end subroutine pmf_paths_write_path

!===============================================================================
! Subroutine:   pmf_paths_write_path_derivatives
! write path derivatives including [PATH] section without header
!===============================================================================

subroutine pmf_paths_write_path_derivatives(iounit,path_item)

    use prmfile
    use pmf_core
    use cv_common

    implicit none
    integer         :: iounit
    type(PathType)  :: path_item
    ! --------------------------------------------
    integer         :: i, b
    ! -----------------------------------------------------------------------------

    write(iounit,10)
    write(iounit,20) trim(path_item%name)
    write(iounit,30) path_item%ncvs
    write(iounit,40) path_item%nbeads

    write(iounit,110,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,140,ADVANCE='NO') i
    end do
    do i=1,path_item%ncvs
        write(iounit,150,ADVANCE='NO') i
    end do
    write(iounit,*)

    write(iounit,120,ADVANCE='NO')
    do i=1,path_item%ncvs * 2
        write(iounit,130,ADVANCE='NO')
    end do
    write(iounit,*)

    write(iounit,100,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,50,ADVANCE='NO') trim(path_item%cvs(i)%cv%name)
    end do
    do i=1,path_item%ncvs
        write(iounit,50,ADVANCE='NO') trim(path_item%cvs(i)%cv%name)
    end do
    write(iounit,*)

    write(iounit,100,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(iounit,50,ADVANCE='NO') trim(path_item%cvs(i)%cv%ctype)
    end do
    do i=1,path_item%ncvs
        write(iounit,50,ADVANCE='NO') trim(path_item%cvs(i)%cv%ctype)
    end do
    write(iounit,*)

    write(iounit,120,ADVANCE='NO')
    do i=1,path_item%ncvs * 2
        write(iounit,130,ADVANCE='NO')
    end do
    write(iounit,*)

    write(iounit,125,ADVANCE='NO')
    do i=1,path_item%ncvs * 2
        write(iounit,80,ADVANCE='NO') i+1
    end do
    write(iounit,*)

    write(iounit,120,ADVANCE='NO')
    do i=1,path_item%ncvs * 2
        write(iounit,130,ADVANCE='NO')
    end do
    write(iounit,*)

    do b=1,path_item%nbeads
        write(iounit,70,ADVANCE='NO') path_item%alphas(b)
        call pmf_paths_get_intpoint_der(path_item,path_item%alphas(b),path_item%cpos)
        do i=1,path_item%ncvs
            write(iounit,60,ADVANCE='NO') path_item%points(b,i)
        end do
        do i=1,path_item%ncvs
            write(iounit,60,ADVANCE='NO') path_item%cpos(i)
        end do
        write(iounit,*)
    end do

    write(iounit,*)

 10 format('# [PATH]')
 20 format('# name   ',A)
 30 format('# ncvs   ',I3)
 40 format('# nbeads ',I3)

 50 format(1X,A14)
 60 format(1X,E14.7)
 70 format(F9.4)
 80 format(1X,I14)

100 format('#        ')
110 format('#  alpha ')
120 format('# -------')
125 format('#       1')
130 format(1X,'--------------')
140 format(1X,'     CV',I2.2,'     ')
150 format(1X,'  dCV',I2.2,'/dalpha')

end subroutine pmf_paths_write_path_derivatives

!===============================================================================
! Function:  pmf_paths_find_path
!===============================================================================

integer function pmf_paths_find_path(path_name)

    use pmf_dat
    use pmf_utils

    implicit none
    character(*)    :: path_name
    ! -----------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    pmf_paths_find_path = 0

    do i=1,NumOfPaths
        if( trim(path_name) .eq. trim(PathList(i)%path%name) ) then
            pmf_paths_find_path = i
            return
        end if
    end do

    call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to find PATH with name: '''//trim(path_name)//'''!')

end function pmf_paths_find_path

!===============================================================================
! Subroutine:   pmf_paths_print_path
!===============================================================================

subroutine pmf_paths_print_path(path_item)

    implicit none
    type(PathType)  :: path_item
    ! --------------------------------------------------------------------------

    call pmf_paths_print_header(path_item)
    call pmf_paths_print_controls(path_item)
    call pmf_paths_print_beads(path_item)

end subroutine pmf_paths_print_path

!===============================================================================
! Subroutine:   pmf_paths_print_header
!===============================================================================

subroutine pmf_paths_print_header(path_item)

    implicit none
    type(PathType)  :: path_item
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)

    write(PMF_OUT,10,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,15,ADVANCE='NO') i
    end do
    write(PMF_OUT,*)

    write(PMF_OUT,20,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,25,ADVANCE='NO')
    end do
    write(PMF_OUT,*)

10 format('#  ID   Type   alpha')
15 format(1X,4X,'CV#',I2.2,5X)

20 format('# ---- ------ ------')
25 format(' --------------')

end subroutine pmf_paths_print_header

!===============================================================================
! Subroutine:   pmf_paths_print_controls
!===============================================================================

subroutine pmf_paths_print_controls(path_item)

    use pmf_core
    use cv_common

    implicit none
    type(PathType)  :: path_item
    ! --------------------------------------------
    integer         :: i
    ! -----------------------------------------------------------------------------

    ! -----------------------------------------------
    write(PMF_OUT,100,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,105,ADVANCE='NO') trim(CVList(path_item%cvindxs(i))%cv%name)
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    write(PMF_OUT,200,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,205,ADVANCE='NO') trim(CVList(path_item%cvindxs(i))%cv%ctype)
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    write(PMF_OUT,250,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,255,ADVANCE='NO') trim(CVList(path_item%cvindxs(i))%cv%get_ulabel())
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    write(PMF_OUT,300,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,305,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%minvalues(i))
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    write(PMF_OUT,400,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,405,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%maxvalues(i))
    end do
    write(PMF_OUT,*)

    ! -----------------------------------------------
    write(PMF_OUT,500,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,505,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%maxmoves(i))
    end do
    write(PMF_OUT,*)

    write(PMF_OUT,600,ADVANCE='NO')
    do i=1,path_item%ncvs
        write(PMF_OUT,605,ADVANCE='NO')
    end do
    write(PMF_OUT,*)

    return

100 format('#      names        ')
105 format(1X,A14)

200 format('#      types        ')
205 format(1X,A14)

250 format('#      units        ')
255 format(1X,A14)

300 format('#      min          ')
305 format(1X,E14.6)

400 format('#      max          ')
405 format(1X,E14.6)

500 format('#      maxmov       ')
505 format(1X,E14.6)

600 format('# ---- ------ ------')
605 format(' --------------')

end subroutine pmf_paths_print_controls

!===============================================================================
! Subroutine:   pmf_paths_print_beads
!===============================================================================

subroutine pmf_paths_print_beads(path_item)

    implicit none
    type(PathType)  :: path_item
    ! --------------------------------------------
    integer         :: b, i
    ! --------------------------------------------------------------------------

    do b=1,path_item%nbeads
        if(  path_item%fixed(b) ) then
            write(PMF_OUT,25,ADVANCE='NO') b,path_item%alphas(b)
        else
            write(PMF_OUT,20,ADVANCE='NO') b,path_item%alphas(b)
        end if
        do i=1,path_item%ncvs
            write(PMF_OUT,30,ADVANCE='NO') CVList(path_item%cvindxs(i))%cv%get_rvalue(path_item%points(b,i))
        end do
        write(PMF_OUT,*)
    end do

    return

20 format(I6,   ' F     ',1X,F6.3)
25 format(I6,   ' P     ',1X,F6.3)
30 format(1X,E14.6)

end subroutine pmf_paths_print_beads

!===============================================================================
! Subroutine:   pmf_paths_optimize_alphas
!===============================================================================

subroutine pmf_paths_optimize_alphas(path_item)

    implicit none
    type(PathType)  :: path_item  
    ! --------------------------------------------
    integer                 :: i,b,alloc_failed
    real(PMFDP)             :: tot_length, slen, path_length, prev_length
    real(PMFDP),allocatable :: lalphas(:)
    ! --------------------------------------------------------------------------

    if( path_item%nbeads .lt. 2 ) then
        call pmf_utils_exit(PMF_OUT,1,'The path must have at least two beads!')
    end if

    allocate(lalphas(path_item%nbeads), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for the path optimization!')
    end if

    ! get initial path length from linear interpolation
    tot_length = 0
    do b=2,path_item%nbeads
        slen = 0
        do i=1,path_item%ncvs
            slen = slen + (path_item%points(b,i)-path_item%points(b-1,i))**2
        end do
        tot_length = tot_length + sqrt(slen)
    end do

    if( tot_length .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'The path has zero length!')
    end if

    ! get initial alphas from linear interpolation
    path_item%alphas(1) = 0.0d0
    path_length = 0
    do b=2,path_item%nbeads-1
        slen = 0
        do i=1,path_item%ncvs
            slen = slen + (path_item%points(b,i)-path_item%points(b-1,i))**2
        end do
        if( slen .eq. 0.0d0 ) then
            call pmf_utils_exit(PMF_OUT,1,'The path segment has zero length!')
        end if
        path_length =  path_length + sqrt(slen)
        path_item%alphas(b) = path_length/tot_length;
    end do
    path_item%alphas(path_item%nbeads) = 1.0

    do while( .true. )
        ! interpolate CVs
        call pmf_paths_get_setint(path_item)

        prev_length = tot_length

        ! determine new path length
        tot_length = 0
        do b=2,path_item%nbeads
            tot_length = tot_length + pmf_paths_get_seglen(path_item,path_item%alphas(b-1),path_item%alphas(b))
        end do

        if( abs(tot_length-prev_length) < 1e-7 ) then
            deallocate(lalphas)
            ! no significant change - quit optimization
            return
        end if

        ! determine new alphas
        lalphas(1) = 0.0
        path_length = 0
        do b=2,path_item%nbeads-1
            path_length =  path_length + pmf_paths_get_seglen(path_item,path_item%alphas(b-1),path_item%alphas(b))
            lalphas(b) = path_length / tot_length
        end do
        lalphas(path_item%nbeads) = 1.0
        path_item%alphas = lalphas

    end do

end subroutine pmf_paths_optimize_alphas

!===============================================================================
! Function:   pmf_paths_get_path_length
!===============================================================================

real(PMFDP) function pmf_paths_get_path_length(path_item)

    implicit none
    type(PathType)          :: path_item
    ! --------------------------------------------
    integer                 :: b
    ! --------------------------------------------------------------------------

    if( path_item%nbeads .lt. 2 ) then
        call pmf_utils_exit(PMF_OUT,1,'The path must have at least two beads!')
    end if

    ! determine path length
    pmf_paths_get_path_length = 0
    do b=2,path_item%nbeads
        pmf_paths_get_path_length = pmf_paths_get_path_length &
           + pmf_paths_get_seglen(path_item,path_item%alphas(b-1),path_item%alphas(b))
    end do

end function pmf_paths_get_path_length

!===============================================================================
! Subroutine:   pmf_paths_get_setint
!===============================================================================

subroutine pmf_paths_get_setint(path_item)

    use pmf_spline

    implicit none
    type(PathType)  :: path_item
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,path_item%ncvs
        call spline_cubic_set(path_item%nbeads, path_item%alphas, path_item%points(:,i), &
                              0, 0.0d0, 0, 0.0d0, path_item%ypp(:,i) )
    end do

end subroutine pmf_paths_get_setint

!===============================================================================
! Function:   pmf_paths_get_intpoint
!===============================================================================

subroutine pmf_paths_get_intpoint(path_item,alpha,point)

    use pmf_spline

    implicit none
    type(PathType)  :: path_item
    real(PMFDP)     :: alpha
    real(PMFDP)     :: point(:)
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,path_item%ncvs
        point(i) = pmf_paths_get_intcv(path_item,i,alpha)
    end do

end subroutine pmf_paths_get_intpoint

!===============================================================================
! Function:   pmf_paths_get_intcv
! cvidx - local CV index
!===============================================================================

real(PMFDP) function pmf_paths_get_intcv(path_item,cvidx,alpha)

    use pmf_spline

    implicit none
    type(PathType)  :: path_item
    integer         :: cvidx
    real(PMFDP)     :: alpha
    ! --------------------------------------------
    real(PMFDP)     :: yval, ypval, yppval
    ! --------------------------------------------------------------------------

    call spline_cubic_val(path_item%nbeads, path_item%alphas, path_item%points(:,cvidx), &
                          path_item%ypp(:,cvidx), &
                          alpha, yval, ypval, yppval)
    pmf_paths_get_intcv =  yval

end function pmf_paths_get_intcv

!===============================================================================
! Function:   pmf_paths_get_intpoint_der
!===============================================================================

subroutine pmf_paths_get_intpoint_der(path_item,alpha,der)

    use pmf_spline

    implicit none
    type(PathType)  :: path_item
    real(PMFDP)     :: alpha
    real(PMFDP)     :: der(:)
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,path_item%ncvs
        der(i) = pmf_paths_get_intcv_der(path_item,i,alpha)
    end do

end subroutine pmf_paths_get_intpoint_der

!===============================================================================
! Function:   pmf_paths_get_intcv
! cvidx - local CV index
!===============================================================================

real(PMFDP) function pmf_paths_get_intcv_der(path_item,cvidx,alpha)

    use pmf_spline

    implicit none
    type(PathType)  :: path_item
    integer         :: cvidx
    real(PMFDP)     :: alpha
    ! --------------------------------------------
    real(PMFDP)     :: yval, ypval, yppval
    ! --------------------------------------------------------------------------

    call spline_cubic_val(path_item%nbeads, path_item%alphas, path_item%points(:,cvidx), &
                          path_item%ypp(:,cvidx), &
                          alpha, yval, ypval, yppval)
    pmf_paths_get_intcv_der =  ypval

end function pmf_paths_get_intcv_der

!===============================================================================
! Function:   pmf_paths_get_seglen
!===============================================================================

real(PMFDP) function pmf_paths_get_seglen(path_item,alpha1,alpha2)

    implicit none
    type(PathType)  :: path_item
    real(PMFDP)     :: alpha1
    real(PMFDP)     :: alpha2
    ! --------------------------------------------
    integer         :: i
    real(PMFDP)     :: len, step, alpha, slen2, curr
    ! --------------------------------------------------------------------------

    len = 0
    do i=1,path_item%ncvs
        path_item%spos(i) = pmf_paths_get_intcv(path_item,i,alpha1)
    end do

    step = (alpha2-alpha1)/real(PMF_PATH_SEGMENT_DISCRETIZATION)
    alpha = alpha1 + step
    do while( alpha .lt. alpha2 )
        slen2 = 0
        do i=1,path_item%ncvs
            curr = pmf_paths_get_intcv(path_item,i,alpha)
            slen2 = slen2 +  (curr-path_item%spos(i))**2
            path_item%spos(i) = curr
        end do
        len = len + sqrt(slen2)
        alpha = alpha + step
    end do

    slen2 = 0
    do i=1,path_item%ncvs
        curr = pmf_paths_get_intcv(path_item,i,alpha2)
        slen2 = slen2 + (curr-path_item%spos(i))**2
    end do
    len = len + sqrt(slen2)

    pmf_paths_get_seglen = len

end function pmf_paths_get_seglen

!===============================================================================
! Function:   pmf_paths_get_alpha
! it numericaly determines alpha that is closest to current CV values
!===============================================================================

real(PMFDP) function pmf_paths_get_alpha(path_item,ctx)

    implicit none
    type(PathType)      :: path_item
    type(CVContextType) :: ctx
    ! --------------------------------------------
    integer             :: i
    real(PMFDP)         :: left_alpha, right_alpha
    real(PMFDP)         :: int_len
    real(PMFDP)         :: min_alpha
    ! --------------------------------------------------------------------------

    left_alpha = 0.0
    right_alpha = 1.0
    do i=0,PMF_PATH_ALPHA_NFOCUS
        int_len = (right_alpha-left_alpha)/PMF_PATH_ALPHA_PRECISION
        min_alpha = pmf_paths_get_alpha_guess(path_item,ctx,left_alpha,right_alpha)
        left_alpha = min_alpha - int_len
        right_alpha = min_alpha + int_len
    end do

    pmf_paths_get_alpha = min_alpha

end function pmf_paths_get_alpha

!===============================================================================
! Function:   pmf_paths_get_alpha_guess
! it numericaly determines alpha that is closest to current CV values
! in the specified interval
!===============================================================================

real(PMFDP) function pmf_paths_get_alpha_guess(path_item,ctx,left_alpha,right_alpha)

    implicit none
    type(PathType)      :: path_item
    type(CVContextType) :: ctx
    real(PMFDP)         :: left_alpha, right_alpha
    ! --------------------------------------------
    integer             :: i
    real(PMFDP)         :: step,dist,min_dist,alpha,min_alpha
    ! --------------------------------------------------------------------------

    if( left_alpha .lt. 0.0 ) left_alpha = 0.0
    if( right_alpha .gt. 1.0 ) right_alpha = 1.0

    alpha = left_alpha
    min_alpha = alpha
    step = (right_alpha - left_alpha)/PMF_PATH_ALPHA_PRECISION
    i = 0
    do while( alpha .le. right_alpha )
        dist = pmf_paths_get_dist(path_item,ctx,alpha)
        if( i .eq. 0 ) then
            min_dist = dist
            min_alpha = alpha
        end if
        if( dist .lt. min_dist) then
            min_dist = dist
            min_alpha = alpha
        end if
        i = i + 1
        alpha = alpha + step
    end do

    pmf_paths_get_alpha_guess = min_alpha

end function pmf_paths_get_alpha_guess

!===============================================================================
! Function:   pmf_paths_get_dist
! get distance from path point to CVs
!===============================================================================

real(PMFDP) function pmf_paths_get_dist(path_item,ctx,alpha)

    implicit none
    type(PathType)      :: path_item
    type(CVContextType) :: ctx
    real(PMFDP)         :: alpha
    ! --------------------------------------------
    integer             :: i
    real(PMFDP)         :: dist,pvalue
    ! --------------------------------------------------------------------------

    dist = 0.0d0

    do i=1,path_item%ncvs
        pvalue = pmf_paths_get_intcv(path_item,i,alpha)
        dist = dist + (pvalue - ctx%CVsValues(path_item%cvindxs(i)))**2
    end do
    dist = sqrt(dist)

    pmf_paths_get_dist = dist

end function pmf_paths_get_dist

!===============================================================================
! Subroutine:   pmf_paths_get_path_current_alpha
!===============================================================================

subroutine pmf_paths_get_path_current_alpha(path_item,ctx)

    use pmf_spline

    implicit none
    type(PathType)      :: path_item
    type(CVContextType) :: ctx
    ! --------------------------------------------
    integer             :: i
    ! --------------------------------------------------------------------------

    do i=1,path_item%ncvs
        path_item%cpos(i) = ctx%CVsValues(path_item%cvindxs(i))
    end do

    call pmf_paths_get_setint(path_item)

    path_item%current_alpha = pmf_paths_get_alpha(path_item,ctx)

end subroutine pmf_paths_get_path_current_alpha

!===============================================================================
! Function:   pmf_paths_get_rpos
! cvidx - global CV index
!===============================================================================

real(PMFDP) function pmf_paths_get_rpos(cvidx)

    use pmf_spline

    implicit none
    integer                 :: cvidx
    ! --------------------------------------------
    class(PathType),pointer :: path_item
    integer                 :: loc_cvidx
    integer                 :: i
    ! --------------------------------------------------------------------------

    pmf_paths_get_rpos = 0.0d0

    if( CVList(cvidx)%cv%pathidx .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'The CV is not associated with the path!')
    end if

    path_item => PathList(CVList(cvidx)%cv%pathidx)%path

    loc_cvidx = 0
    do i=1,path_item%ncvs
        if( path_item%cvindxs(i) .eq. cvidx ) then
            loc_cvidx = i
            exit
        end if
    end do

    if( loc_cvidx .eq. 0 ) then
        write(PMF_OUT,*) 'cvidx = ', cvidx, ' pathidx = ', CVList(cvidx)%cv%pathidx, &
                         ' pathidx = ', path_item%idx
        call pmf_utils_exit(PMF_OUT,1,'Unable to translate CV index to path local index!')
    end if

    pmf_paths_get_rpos = path_item%rpos(loc_cvidx)

end function pmf_paths_get_rpos

!===============================================================================

end module pmf_paths
