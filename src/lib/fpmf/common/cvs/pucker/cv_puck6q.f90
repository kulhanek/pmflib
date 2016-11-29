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

! 6-membered ring puckering - total puckering amplitude Q

module cv_puck6q

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePUCK6Q
    contains
        procedure :: load_cv        => load_puck6q
        procedure :: calculate_cv   => calculate_puck6q
end type CVTypePUCK6Q

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_puck6q
!===============================================================================

subroutine load_puck6q(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePUCK6Q)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: i,j
    ! --------------------------------------------------------------------------

    ! simple init and allocation --------------------
    cv_item%ctype         = 'PUCK6Q'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 6
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

end subroutine load_puck6q

!===============================================================================
! Subroutine:  calculate_puck6q
!===============================================================================

subroutine calculate_puck6q(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypePUCK6Q) :: cv_item
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
    do  i = 1, 6
        ai = cv_item%lindexes(i)
        com(:) = com(:) + x(:,ai)
    end do
    com(:) = com(:) / 6.0d0

    ! calculate r1 and r2
    r1(:) = 0.0d0
    r2(:) = 0.0d0
    do  i = 1, 6
        ai = cv_item%lindexes(i)
        dx(:) = x(:,ai) - com(:)
        r1(:) = r1(:) + dx(:)*sin(2.0d0/6.0d0*PMF_PI*(i-1))
        r2(:) = r2(:) + dx(:)*cos(2.0d0/6.0d0*PMF_PI*(i-1))
    end do

    ! calculate plane axis - normalized cross-product
    n = vcross(r1,r2)
    n_len = sqrt( n(1)**2 + n(2)**2 + n(3)**2 )
    n(:) = n(:) / n_len

    ! calculate Q
    q2 = 0.0d0
    do  j = 1, 6
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


    ! block([%1,%2,%3,%4,%5,%6,%7,%8,%9,%10,%11,%12,%13,%14,%15,%16,%17,%18,%19,%20,%21,%22,%23,%24,%25,%26,%27,%28,%29,%30,%31]
    ! %2:-(r6x+r5x+r4x+r3x+r2x+r1x)/6
    ! %3:%2+r1x
    ! %4:%2+r2x
    ! %5:%2+r3x
    ! %6:%2+r4x
    ! %7:%2+r5x
    ! %8:%2+r6x
    ! %9:-(r6y+r5y+r4y+r3y+r2y+r1y)/6
    ! %10:%9+r1y
    ! %11:%9+r2y
    ! %12:%9+r3y
    ! %13:%9+r4y
    ! %14:%9+r5y
    ! %15:%9+r6y
    ! %16:-(r6z+r5z+r4z+r3z+r2z+r1z)/6
    ! %17:%16+r1z
    ! %18:%16+r2z
    ! %19:%16+r3z
    ! %20:%16+r4z
    ! %21:%16+r5z
    ! %22:%16+r6z
    ! %23:[%8*v6+%7*v5+%6*v4+%5*v3+%4*v2+%3*v1,%15*v6+%14*v5+%13*v4+%12*v3+%11*v2+%10*v1,%22*v6+%21*v5+%20*v4+%19*v3+%18*v2+%17*v1]
    ! %23 = r2
    ! %24:[%8*w6+%7*w5+%6*w4+%5*w3+%4*w2+%3*w1,%15*w6+%14*w5+%13*w4+%12*w3+%11*w2+%10*w1,%22*w6+%21*w5+%20*w4+%19*w3+%18*w2+%17*w1]
    ! %24 = r1
    ! %25:%23~%24
    p25 = vcross(r2,r1)
    ! %27:%25~%24
    p27 = vcross(p25,r1)

    ! %31:sqrt(-%23 . %27)
    ! p31 = sqrt( - vdot(r2,p27) )
    p31 = n_len

    do  j = 1, 6  ! over z_j
        aj = cv_item%lindexes(j)

        dx(:) = x(:,aj) - com(:)
        zj = vdot(dx,n)

        ! %1:[r1x-comx,r1y-comy,r1z-comz]
        p1 = dx

        fac = sc*zj

        do  i = 1, 6  ! over atoms
            ai = cv_item%lindexes(i)

            wx = 0.0
            vx = 0.0
            do k = 1, 6
                wx = wx + (kdelta(k,i)-1.0d0/6.0d0)*sin(2.0d0/6.0d0*PMF_PI*(k-1))
                vx = vx + (kdelta(k,i)-1.0d0/6.0d0)*cos(2.0d0/6.0d0*PMF_PI*(k-1))
            end do

            ! -------------------------
            ! %26:[-v6/6-v5/6-v4/6-v3/6-v2/6+(5*v1)/6,0,0]
            p26 = 0;
            p26(1) = vx

            ! %29:[-w6/6-w5/6-w4/6-w3/6-w2/6+(5*w1)/6,0,0]
            p29 = 0;
            p29(1) = wx

            ! %28:%26~%24
            p28 = vcross(p26,r1)
            ! %30:%23~%29
            p30 = vcross(r2,p29)

            ! ((%1 . %25)*(-%23 . %25~%29 -%23 . (%30+%28)~%24 - %26 . %27))/(2*%31^3)-(%1 . %30 + %1 . %28 +[1,0,0] . %25)/%31)
            p32 = vcross(p25,p29)
            p33 = p30+p28
            p34 = vcross(p33,r1)
            ddx = (vdot(p1,p25)*(- vdot(r2,p32) - vdot(r2,p34) - vdot(p26,p27)))/(2.0d0*p31**3)
            ddx = ddx -(vdot(p1,p30) + vdot(p1,p28) + (kdelta(j,i)-1.0d0/6.0d0)*p25(1))/(p31)

            ! -------------------------
            ! %26:[-v6/6-v5/6-v4/6-v3/6-v2/6+(5*v1)/6,0,0]
            p26 = 0;
            p26(2) = vx

            ! %29:[-w6/6-w5/6-w4/6-w3/6-w2/6+(5*w1)/6,0,0]
            p29 = 0;
            p29(2) = wx

            ! %28:%26~%24
            p28 = vcross(p26,r1)
            ! %30:%23~%29
            p30 = vcross(r2,p29)

            ! ((%1 . %25)*(-%23 . %25~%29 -%23 . (%30+%28)~%24 - %26 . %27))/(2*%31^3)-(%1 . %30 + %1 . %28 +[1,0,0] . %25)/%31)
            p32 = vcross(p25,p29)
            p33 = p30+p28
            p34 = vcross(p33,r1)
            ddy = (vdot(p1,p25)*(- vdot(r2,p32) - vdot(r2,p34) - vdot(p26,p27)))/(2.0d0*p31**3)
            ddy = ddy -(vdot(p1,p30) + vdot(p1,p28) + (kdelta(j,i)-1.0d0/6.0d0)*p25(2))/(p31)

            ! -------------------------
            ! %26:[-v6/6-v5/6-v4/6-v3/6-v2/6+(5*v1)/6,0,0]
            p26 = 0;
            p26(3) = vx

            ! %29:[-w6/6-w5/6-w4/6-w3/6-w2/6+(5*w1)/6,0,0]
            p29 = 0;
            p29(3) = wx

            ! %28:%26~%24
            p28 = vcross(p26,r1)
            ! %30:%23~%29
            p30 = vcross(r2,p29)

            ! ((%1 . %25)*(-%23 . %25~%29 -%23 . (%30+%28)~%24 - %26 . %27))/(2*%31^3)-(%1 . %30 + %1 . %28 +[1,0,0] . %25)/%31)
            p32 = vcross(p25,p29)
            p33 = p30+p28
            p34 = vcross(p33,r1)
            ddz = (vdot(p1,p25)*(- vdot(r2,p32) - vdot(r2,p34) - vdot(p26,p27)))/(2.0d0*p31**3)
            ddz = ddz -(vdot(p1,p30) + vdot(p1,p28) + (kdelta(j,i)-1.0d0/6.0d0)*p25(3))/(p31)

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

end subroutine calculate_puck6q

!===============================================================================

end module cv_puck6q

