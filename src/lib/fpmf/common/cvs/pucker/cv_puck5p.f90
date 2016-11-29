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

! 5-membered ring puckering - phi - <0,2*PI)

module cv_puck5p

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePUCK5P
    contains
        procedure :: load_cv        => load_puck5p
        procedure :: calculate_cv   => calculate_puck5p
        ! PBC related methods
        procedure   :: is_periodic_cv       => is_periodic_cv_puck5p
        procedure   :: get_min_cv_value     => get_min_cv_value_puck5p
        procedure   :: get_max_cv_value     => get_max_cv_value_puck5p
end type CVTypePUCK5P

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_puck5p
!===============================================================================

subroutine load_puck5p(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePUCK5P)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: i,j
    ! --------------------------------------------------------------------------

    ! simple init and allocation --------------------
    cv_item%ctype         = 'PUCK5P'
    cv_item%unit          = AngleUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 5
    call cv_common_read_groups(cv_item,prm_fin)

    ! every group can contain only one atom
    do i=1,cv_item%ngrps
        if( cv_item%grps(i) .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,&
                                'Each atom group can contain only one atom!')
        end if
    end do

    ! check overlaps
    do i=1,cv_item%ngrps
        do j=i+1,cv_item%ngrps
            if( cv_item%rindexes(i) .eq. cv_item%rindexes(j) ) then
                call pmf_utils_exit(PMF_OUT,1,&
                                    'Atom groups cannot share atoms!')
            end if
        end do
    end do

end subroutine load_puck5p

!===============================================================================
! Function:  is_periodic_cv_puck5p
!===============================================================================

logical function is_periodic_cv_puck5p(cv_item)

    implicit none
    class(CVTypePUCK5P) :: cv_item
    ! --------------------------------------------------------------------------

    is_periodic_cv_puck5p = .true.

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function is_periodic_cv_puck5p

!===============================================================================
! Function:  get_min_cv_value_puck5p
!===============================================================================

real(PMFDP) function get_min_cv_value_puck5p(cv_item)

    implicit none
    class(CVTypePUCK5P) :: cv_item
    ! --------------------------------------------------------------------------

    get_min_cv_value_puck5p = 0.0d0

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function get_min_cv_value_puck5p

!===============================================================================
! Function:  get_max_cv_value_puck5p
!===============================================================================

real(PMFDP) function get_max_cv_value_puck5p(cv_item)

    implicit none
    class(CVTypePUCK5P) :: cv_item
    ! --------------------------------------------------------------------------

    get_max_cv_value_puck5p = 2.0d0*PMF_PI

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function get_max_cv_value_puck5p

!===============================================================================
! Subroutine:  calculate_puck5p
!===============================================================================

subroutine calculate_puck5p(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypePUCK5P) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ai,j,aj,k
    real(PMFDP)         :: com(3)
    real(PMFDP)         :: r1(3),r2(3),dx(3)
    real(PMFDP)         :: n(3),n_len
    real(PMFDP)         :: zj,qx,qy,qxt,qyt
    real(PMFDP)         :: p1(3),p25(3),p26(3),p27(3),p28(3),p29(3),p30(3)
    real(PMFDP)         :: p32(3),p33(3),p34(3)
    real(PMFDP)         :: p31,wx,vx,ddx,ddy,ddz,fac
    ! --------------------------------------------------------------------------

    ! calculate center of the ring
    com(:) = 0.0d0
    do  i = 1, 5
        ai = cv_item%lindexes(i)
        com(:) = com(:) + x(:,ai)
    end do
    com(:) = com(:) / 5.0d0

    ! calculate r1 and r2
    r1(:) = 0.0d0
    r2(:) = 0.0d0
    do  i = 1, 5
        ai = cv_item%lindexes(i)
        dx(:) = x(:,ai) - com(:)
        r1(:) = r1(:) + dx(:)*sin(2.0d0/5.0d0*PMF_PI*(i-1))
        r2(:) = r2(:) + dx(:)*cos(2.0d0/5.0d0*PMF_PI*(i-1))
    end do

    ! calculate plane axis - normalized cross-product
    n = vcross(r1,r2)
    n_len = sqrt( n(1)**2 + n(2)**2 + n(3)**2 )
    n(:) = n(:) / n_len

    ! calculate Q
    qxt = 0.0d0
    qyt = 0.0d0
    do  j = 1, 5
        aj = cv_item%lindexes(j)
        dx(:) = x(:,aj) - com(:)
        zj = vdot(dx,n)
        qxt = qxt + zj*cos(4.0d0/5.0d0*PMF_PI*(j-1))
        qyt = qyt - zj*sin(4.0d0/5.0d0*PMF_PI*(j-1))
    end do

    qx = qxt * sqrt(2.0d0/5.0d0)
    qy = qyt * sqrt(2.0d0/5.0d0)

    ! phi
    ctx%CVsValues(cv_item%idx) = atan2(qy,qx)

    ! transform to 0,2PI interval
    if( ctx%CVsValues(cv_item%idx) .lt. 0 ) then
        ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + 2*PMF_PI
    end if

    ! ------------------------------------------------
    ! derivatives

    p25 = vcross(r2,r1)
    p27 = vcross(p25,r1)
    p31 = n_len

    do  j = 1, 5  ! over z_j
        aj = cv_item%lindexes(j)

        dx(:) = x(:,aj) - com(:)
        zj = vdot(dx,n)

        p1 = dx

        wx =   cos(4.0d0/5.0d0*PMF_PI*(j-1))
        vx = - sin(4.0d0/5.0d0*PMF_PI*(j-1))

        fac = (vx*qxt - wx*qyt)/(qxt**2+qyt**2)

        do  i = 1, 5  ! over atoms
            ai = cv_item%lindexes(i)

            wx = 0.0
            vx = 0.0
            do k = 1, 5
                wx = wx + (kdelta(k,i)-1.0d0/5.0d0)*sin(2.0d0/5.0d0*PMF_PI*(k-1))
                vx = vx + (kdelta(k,i)-1.0d0/5.0d0)*cos(2.0d0/5.0d0*PMF_PI*(k-1))
            end do

            ! -------------------------
            p26 = 0;
            p26(1) = vx
            p29 = 0;
            p29(1) = wx
            p28 = vcross(p26,r1)
            p30 = vcross(r2,p29)
            p32 = vcross(p25,p29)
            p33 = p30+p28
            p34 = vcross(p33,r1)
            ddx = (vdot(p1,p25)*(- vdot(r2,p32) - vdot(r2,p34) - vdot(p26,p27)))/(2.0d0*p31**3)
            ddx = ddx -(vdot(p1,p30) + vdot(p1,p28) + (kdelta(j,i)-1.0d0/5.0d0)*p25(1))/(p31)

            ! -------------------------
            p26 = 0;
            p26(2) = vx
            p29 = 0;
            p29(2) = wx
            p28 = vcross(p26,r1)
            p30 = vcross(r2,p29)
            p32 = vcross(p25,p29)
            p33 = p30+p28
            p34 = vcross(p33,r1)
            ddy = (vdot(p1,p25)*(- vdot(r2,p32) - vdot(r2,p34) - vdot(p26,p27)))/(2.0d0*p31**3)
            ddy = ddy -(vdot(p1,p30) + vdot(p1,p28) + (kdelta(j,i)-1.0d0/5.0d0)*p25(2))/(p31)

            ! -------------------------
            p26 = 0;
            p26(3) = vx
            p29 = 0;
            p29(3) = wx
            p28 = vcross(p26,r1)
            p30 = vcross(r2,p29)
            p32 = vcross(p25,p29)
            p33 = p30+p28
            p34 = vcross(p33,r1)
            ddz = (vdot(p1,p25)*(- vdot(r2,p32) - vdot(r2,p34) - vdot(p26,p27)))/(2.0d0*p31**3)
            ddz = ddz -(vdot(p1,p30) + vdot(p1,p28) + (kdelta(j,i)-1.0d0/5.0d0)*p25(3))/(p31)

            ! final
            ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                    + fac*ddx
            ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                    + fac*ddy
            ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                    + fac*ddz

        end do
    end do

 return

end subroutine calculate_puck5p

!===============================================================================

end module cv_puck5p

