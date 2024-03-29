!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2008 Silvia Cereda, sc578@cam.ac.uk &
!                       Petr Kulhanek, kulhanek@enzim.hu
!
!    This library is free software; you can repostribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is postributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module cv_opos

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeOPOS

    integer :: z_direction
    integer :: direction

    contains
        procedure :: load_cv        => load_opos
        procedure :: calculate_cv   => calculate_opos
end type CVTypeOPOS

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_opos
!===============================================================================

subroutine load_opos(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypeOPOS)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: m
    character(1)                        :: cdir
    logical                             :: found
    character(len=PRMFILE_MAX_LINE)     :: mask
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'OPOS'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! init groups and control array -----------------
    cv_item%ngrps = 2
    call cv_common_init_groups(cv_item,prm_fin)

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    ! load z_direction ------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'z_direction',mask) ) then
        call pmf_utils_exit(PMF_OUT,1,'z_direction atom is not specified!')
    end if
    write(PMF_OUT,70) trim(mask)
    call cv_common_set_atom_from_mask(mask,cv_item%z_direction)

    ! find atom in cv_item%rindexes
    found = .false.
    do m=1,cv_item%grps(1)
        if( cv_item%rindexes(m) .eq. cv_item%z_direction ) then
            cv_item%z_direction = m
            found = .true.
            exit
        end if
    end do

    if( .not. found ) then
        call pmf_utils_exit(PMF_OUT,1,'z_direction atom has to be member of group_a!')
    end if

    ! read group b ----------------------------------
    write(PMF_OUT,90)
    call cv_common_read_group(cv_item,prm_fin,2)

    ! load orientation ------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'direction',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'direction is not specified!')
    else
        write(PMF_OUT,60) cdir
    end if

    select case(cdir)
        case('x')
            cv_item%direction = 1
        case('y')
            cv_item%direction = 2
        case('z')
            cv_item%direction = 3
        case default
            call pmf_utils_exit(PMF_OUT,1,&
                                'direction has to be x, y, or z!')
    end select

    return

50 format('   == System A ===================================')
70 format('   ** z-axis ctrl atom   : ',A)
90 format('   == Point B ====================================')
60 format('   ** z-axis direction   : ',A1)

end subroutine load_opos

!===============================================================================
! Subroutine:  calculate_opos
!===============================================================================

subroutine calculate_opos(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeOPOS)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: m,mi,mj,ai,di,info
    real(PMFDP)    :: amass,zvec_s,ac,msign(3)
    real(PMFDP)    :: totmass1,com(3)
    real(PMFDP)    :: totmass2,pos(3)
    real(PMFDP)    :: px(3), mit(3,3)
    real(PMFDP)    :: eigenvalues(3)
    real(PMFDP)    :: work(26*3)
    real(PMFDP)    :: vx(3)
    real(PMFDP)    :: zvec(3)
    real(PMFDP)    :: v(3,3),api(3,3),cij(3),xij(3,3,3),bint(3,3)
    ! --------------------------------------------------------------------------

    ! calculate centre of mass for first group
    totmass1 = 0.0d0
    com(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        com(:) = com(:) + x(:,ai)*amass
        totmass1 = totmass1 + amass
    end do

    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_opos!')
    end if

    com(:) = com(:) / totmass1

    ! calculate centre of mass for second group
    totmass2 = 0.0d0
    pos(:) = 0.0
    do  m = cv_item%grps(1)+1, cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        pos(:) = pos(:) + x(:,ai)*amass
        totmass2 = totmass2 + amass
    end do

    if( totmass2 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass2 is zero in calculate_opos!')
    end if

    pos(:) = pos(:) / totmass2

    ! calculate moment of inertia tensor
    mit(:,:) = 0.0d0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        px(:) = x(:,ai) - com(:)
        mit(1,1) = mit(1,1) + amass*(px(2)**2 + px(3)**2)
        mit(2,2) = mit(2,2) + amass*(px(1)**2 + px(3)**2)
        mit(3,3) = mit(3,3) + amass*(px(1)**2 + px(2)**2)
        mit(1,2) = mit(1,2) - amass*px(1)*px(2)
        mit(1,3) = mit(1,3) - amass*px(1)*px(3)
        mit(2,3) = mit(2,3) - amass*px(2)*px(3)
    end do
    mit(2,1) = mit(1,2)
    mit(3,1) = mit(1,3)
    mit(3,2) = mit(2,3)

    ! find principal moments
    ! calculate eignevalues and eigenvectors of mit
    eigenvalues(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 3, mit, 3, eigenvalues, work, 26*3, info)

    if( info .ne. 0 ) then
       call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize MIT matrix in calculate_opos!')
    end if

    ! determine z-axis orientation - z-control vector
    zvec(:) = x(:,cv_item%lindexes(cv_item%z_direction)) - com(:)
    zvec_s = sqrt(zvec(1)**2 + zvec(2)**2 + zvec(3)**2)

    if( zvec_s .lt. 0.1d0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Control vector in z-direction is too short (< 0.1A)!')
    end if

    zvec(:) = zvec(:)/zvec_s

    ! angle between z-axis and principal z-axis
    ac = zvec(1)*mit(1,3) + zvec(2)*mit(2,3) + zvec(3)*mit(2,3)

    ! TODO - make system right-handed
    msign(1) = 1.0d0
    msign(2) = 1.0d0
    msign(3) = sign(1.0d0,ac)

    ! move origin to com and rotate coordinate system
    px(:) = pos(:) - com(:)
    vx(:) = 0.0d0
    do m=1,3
        vx(1) = vx(1) + px(m)*mit(m,1)
        vx(2) = vx(2) + px(m)*mit(m,2)
        vx(3) = vx(3) + px(m)*mit(m,3)
    end do

    di = cv_item%direction
    ctx%CVsValues(cv_item%idx) = vx(di)*msign(di)

    ! calculate forces -------------------------------

    ! first part, e.g. a*dx'
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - msign(di)*mit(:,di)*amass/totmass1
    end do

    do  m = cv_item%grps(1) + 1 , cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + msign(di)*mit(:,di)*amass/totmass2
    end do

    ! eigenvector derivatives -----------------------

    ! construct pseudoinverse matrix, api
    v(:,:) = mit(:,:)
    api(:,:) = 0.0d0
    do mi=1,3
       if( mi .ne. di ) api(mi,mi) = 1.0d0/(eigenvalues(mi) - eigenvalues(di))
    end do
    call dgemm('N','N',3,3,3,1.0d0,v,3,api,3,0.0d0,bint,3)
    call dgemm('N','T',3,3,3,1.0d0,bint,3,v,3,0.0d0,api,3)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,3
       do mj=1,3
           ! construct cij
           cij(:) = 0.0d0
           cij(mi) = cij(mi) + mit(mj,di)

           ! find eigenvector derivatives
           ! xi contains derivatives of eigenvector by A_ij element
           call dgemv('N',3,3,-1.0d0,api,3,cij,1,0.0d0,xij(:,mi,mj),1)

           ! multiply by px
           xij(:,mi,mj) = xij(:,mi,mj)*px(:)*msign(di)
       end do
    end do

    ! and finaly gradients --------------------------
    do m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + amass*( &
             2.0d0*(x(1,ai)-com(1))*(xij(1,2,2)+xij(1,3,3)+xij(2,2,2)+xij(2,3,3)+xij(3,2,2)+xij(3,3,3)) &
           -       (x(2,ai)-com(2))*(xij(1,1,2)+xij(1,2,1)+xij(2,1,2)+xij(2,2,1)+xij(3,1,2)+xij(3,2,1)) &
           -       (x(3,ai)-com(3))*(xij(1,1,3)+xij(1,3,1)+xij(2,1,3)+xij(2,3,1)+xij(3,1,3)+xij(3,3,1)) )

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + amass*( &
           -       (x(1,ai)-com(1))*(xij(1,1,2)+xij(1,2,1)+xij(2,1,2)+xij(2,2,1)+xij(3,1,2)+xij(3,2,1)) &
           + 2.0d0*(x(2,ai)-com(2))*(xij(1,1,1)+xij(1,3,3)+xij(2,1,1)+xij(2,3,3)+xij(3,1,1)+xij(3,3,3)) &
           -       (x(3,ai)-com(3))*(xij(1,2,3)+xij(1,3,2)+xij(2,2,3)+xij(2,3,2)+xij(3,2,3)+xij(3,3,2)) )

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + amass*( &
           -       (x(1,ai)-com(1))*(xij(1,1,3)+xij(1,3,1)+xij(2,1,3)+xij(2,3,1)+xij(3,1,3)+xij(3,3,1)) &
           -       (x(2,ai)-com(2))*(xij(1,2,3)+xij(1,3,2)+xij(2,2,3)+xij(2,3,2)+xij(3,2,3)+xij(3,3,2)) &
           + 2.0d0*(x(3,ai)-com(3))*(xij(1,1,1)+xij(1,2,2)+xij(2,1,1)+xij(2,2,2)+xij(3,1,1)+xij(3,2,2)) )

    end do

end subroutine calculate_opos

!===============================================================================

end module cv_opos

