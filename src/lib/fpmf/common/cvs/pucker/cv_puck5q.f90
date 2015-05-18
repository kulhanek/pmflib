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

! 5-membered ring puckering - total puckering amplitude Q

module cv_puck5q

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePUCK5Q
    contains
        procedure :: load_cv        => load_puck5q
        procedure :: calculate_cv   => calculate_puck5q
end type CVTypePUCK5Q

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_puck5q
!===============================================================================

subroutine load_puck5q(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common
    use pmf_utils

    implicit none
    class(CVTypePUCK5Q)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: i,j
    ! --------------------------------------------------------------------------

    ! simple init and allocation --------------------
    cv_item%ctype         = 'PUCK5Q'
    cv_item%unit          = LengthUnit
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

end subroutine load_puck5q

!===============================================================================
! Subroutine:  calculate_puck5q
!===============================================================================

subroutine calculate_puck5q(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypePUCK5Q) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ai,j,aj,k
    real(PMFDP)         :: com(3)
    real(PMFDP)         :: r1(3),r2(3),dx(3)
    real(PMFDP)         :: n(3),n_len
    real(PMFDP)         :: zj,q2,sc,value
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
    q2 = 0.0d0
    do  j = 1, 5
        aj = cv_item%lindexes(j)
        dx(:) = x(:,aj) - com(:)
        zj = vdot(dx,n)
        q2 = q2 + zj**2
    end do

    value = sqrt(q2)
    ctx%CVsValues(cv_item%idx) = value

    ! ------------------------------------------------
    ! derivatives
    if( value .gt. 1.0d-7 ) then
        sc = 1.0d0/value
    else
        sc = 0.0d0
    end if

    p25 = vcross(r2,r1)
    p27 = vcross(p25,r1)
    p31 = n_len

    do  j = 1, 5  ! over z_j
        aj = cv_item%lindexes(j)

        dx(:) = x(:,aj) - com(:)
        zj = vdot(dx,n)

        ! %1:[r1x-comx,r1y-comy,r1z-comz]
        p1 = dx

        fac = sc*zj

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

end subroutine calculate_puck5q

!===============================================================================

end module cv_puck5q

