!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module pmf_pbc

 implicit none
 contains

!===============================================================================
! Subroutine:  pmf_pbc_set_box
!===============================================================================

subroutine pmf_pbc_set_box(a,b,c,alpha,beta,gamma)

    use pmf_dat
    use pmf_utils

    real(PMFDP)    :: a,b,c
    real(PMFDP)    :: alpha,beta,gamma
    ! -----------------------------------------------
    integer        :: i
    real(PMFDP)    :: dist,u12(3),u23(3),u31(3)
    real(PMFDP)    :: fbox_abc(3),fbox_ang(3)
    ! -----------------------------------------------------------------------------

    fbox_abc(1) = a*LengthConv
    fbox_abc(2) = b*LengthConv
    fbox_abc(3) = c*LengthConv

    fbox_ang(1) = alpha*AngleConv
    fbox_ang(2) = beta*AngleConv
    fbox_ang(3) = gamma*AngleConv

    ! calculate cell parameters ---------------------
    fucell(:,:) = 0.0d0
    frecip(:,:) = 0.0d0
    fbox_volume = 0.0d0
    fbox_sphere = 0.0d0

    fbox_type = BOX_ISOLATED_SYSTEM

    if( fsystype .eq. SYS_UNK .or. fsystype .eq. SYS_NT ) return

    ! determine box type ----------------------------
    if( (fbox_abc(1) .eq. 0.0d0) .or. &
     (fbox_abc(2) .eq. 0.0d0) .or. &
     (fbox_abc(3) .eq. 0.0d0) .or. &
     (fbox_ang(1) .eq. 0.0d0) .or. &
     (fbox_ang(2) .eq. 0.0d0) .or. &
     (fbox_ang(3) .eq. 0.0d0) ) then
        fbox_type = BOX_ISOLATED_SYSTEM
    else if( (abs(fbox_ang(1) - PMF_HPI) .lt. 0.01d0) .and. &
          (abs(fbox_ang(2) - PMF_HPI) .lt. 0.01d0) .and. &
          (abs(fbox_ang(3) - PMF_HPI) .lt. 0.01d0) ) then
        fbox_type = BOX_ORTHOGONAL
        fucell(1,1) = fbox_abc(1)
        fucell(2,2) = fbox_abc(2)
        fucell(3,3) = fbox_abc(3)
        frecip(1,1) = 1.0d0 / fucell(1,1)
        frecip(2,2) = 1.0d0 / fucell(2,2)
        frecip(3,3) = 1.0d0 / fucell(3,3)

        fbox_volume = fucell(1,1)*fucell(2,2)*fucell(3,3)
        fbox_sphere = 0.5d0*min(fucell(1,1),fucell(2,2),fucell(3,3))
    else
        fbox_type = BOX_GENERAL

        fucell(1,1) = fbox_abc(1)

        fucell(1,2) = fbox_abc(2)*cos(fbox_ang(3))
        fucell(2,2) = fbox_abc(2)*sin(fbox_ang(3))

        fucell(1,3) = fbox_abc(3)*cos(fbox_ang(2))
        fucell(2,3) = (fbox_abc(2)*fbox_abc(3)*cos(fbox_ang(1)) - fucell(1,3)*fucell(1,2))/fucell(2,2)
        fucell(3,3) = sqrt(fbox_abc(3)**2 - fucell(1,3)**2 - fucell(2,3)**2)

        u23(1) = fucell(2,2)*fucell(3,3) - fucell(3,2)*fucell(2,3)
        u23(2) = fucell(3,2)*fucell(1,3) - fucell(1,2)*fucell(3,3)
        u23(3) = fucell(1,2)*fucell(2,3) - fucell(2,2)*fucell(1,3)

        u31(1) = fucell(2,3)*fucell(3,1) - fucell(3,3)*fucell(2,1)
        u31(2) = fucell(3,3)*fucell(1,1) - fucell(1,3)*fucell(3,1)
        u31(3) = fucell(1,3)*fucell(2,1) - fucell(2,3)*fucell(1,1)

        u12(1) = fucell(2,1)*fucell(3,2) - fucell(3,1)*fucell(2,2)
        u12(2) = fucell(3,1)*fucell(1,2) - fucell(1,1)*fucell(3,2)
        u12(3) = fucell(1,1)*fucell(2,2) - fucell(2,1)*fucell(1,2)

        fbox_volume = fucell(1,1)*u23(1) + fucell(2,1)*u23(2)+ fucell(3,2)*u23(3)

        frecip(1,:) = u23(:) / fbox_volume
        frecip(2,:) = u31(:) / fbox_volume
        frecip(3,:) = u12(:) / fbox_volume

        fbox_sphere = fbox_abc(1)*fbox_abc(2)*fbox_abc(3)
        do i = 1, 3
            dist = frecip(i,1)*fucell(1,i) + frecip(i,2)*fucell(2,i) + frecip(i,3)*fucell(3,i)
            dist = dist / sqrt( frecip(i,1)**2 + frecip(i,2)**2 + frecip(i,3)**2)
            if( dist .lt. fbox_sphere ) fbox_sphere = dist
        end do
        fbox_sphere = 0.5d0 * fbox_sphere
    end if

    p0VEne = fpressure * fbox_volume * PMF_PVCONV

    return

end subroutine pmf_pbc_set_box

!===============================================================================
! Subroutine:  pmf_pbc_set_box_from_lvectors
!===============================================================================

subroutine pmf_pbc_set_box_from_lvectors(lattice)

    use pmf_dat
    use pmf_utils

    real(PMFDP)    :: lattice(3,3)
    ! -----------------------------------------------
    integer        :: i
    real(PMFDP)    :: u12(3),u23(3),u31(3)
    real(PMFDP)    :: a, b, c, dist
    real(PMFDP)    :: calpha,cbeta,cgamma
    ! -----------------------------------------------------------------------------

    ! direct mapping
    fucell(:,:) = lattice(:,:)*LengthConv

    ! calculate cell parameters ---------------------
    frecip(:,:) = 0.0d0
    fbox_volume = 0.0d0
    fbox_sphere = 0.0d0

    fbox_type = BOX_ISOLATED_SYSTEM

    if( fsystype .eq. SYS_UNK .or. fsystype .eq. SYS_NT ) return

    a = sqrt(fucell(1,1)*fucell(1,1) + fucell(2,1)*fucell(2,1) + fucell(3,1)*fucell(3,1))
    b = sqrt(fucell(1,2)*fucell(1,2) + fucell(2,2)*fucell(2,2) + fucell(3,2)*fucell(3,2))
    c = sqrt(fucell(1,3)*fucell(1,3) + fucell(2,3)*fucell(2,3) + fucell(3,3)*fucell(3,3))

    calpha = fucell(1,1)*fucell(1,2) + fucell(2,1)*fucell(2,2) + fucell(3,1)*fucell(3,2)
    calpha = calpha / (a*b)

    cbeta = fucell(1,2)*fucell(1,3) + fucell(2,2)*fucell(2,3) + fucell(3,2)*fucell(3,3)
    calpha = calpha / (b*c)

    cgamma = fucell(1,1)*fucell(1,3) + fucell(2,1)*fucell(2,3) + fucell(3,1)*fucell(3,3)
    calpha = calpha / (a*c)

    ! determine box type ----------------------------
    if( (abs(a) .lt. 0.00001d0) .or. &
     (abs(b) .lt. 0.00001d0) .or. &
     (abs(c) .lt. 0.00001d0)  ) then
        fbox_type = BOX_ISOLATED_SYSTEM
    else if( (abs(calpha) .lt. 0.00001d0) .and. &
          (abs(cbeta)  .lt. 0.00001d0) .and. &
          (abs(cgamma) .lt. 0.00001d0) ) then
        fbox_type = BOX_ORTHOGONAL

        frecip(1,1) = 1.0d0 / fucell(1,1)
        frecip(2,2) = 1.0d0 / fucell(2,2)
        frecip(3,3) = 1.0d0 / fucell(3,3)

        fbox_volume = fucell(1,1)*fucell(2,2)*fucell(3,3)
        fbox_sphere = 0.5d0*min(fucell(1,1),fucell(2,2),fucell(3,3))
    else
        fbox_type = BOX_GENERAL

        u23(1) = fucell(2,2)*fucell(3,3) - fucell(3,2)*fucell(2,3)
        u23(2) = fucell(3,2)*fucell(1,3) - fucell(1,2)*fucell(3,3)
        u23(3) = fucell(1,2)*fucell(2,3) - fucell(2,2)*fucell(1,3)

        u31(1) = fucell(2,3)*fucell(3,1) - fucell(3,3)*fucell(2,1)
        u31(2) = fucell(3,3)*fucell(1,1) - fucell(1,3)*fucell(3,1)
        u31(3) = fucell(1,3)*fucell(2,1) - fucell(2,3)*fucell(1,1)

        u12(1) = fucell(2,1)*fucell(3,2) - fucell(3,1)*fucell(2,2)
        u12(2) = fucell(3,1)*fucell(1,2) - fucell(1,1)*fucell(3,2)
        u12(3) = fucell(1,1)*fucell(2,2) - fucell(2,1)*fucell(1,2)

        fbox_volume = fucell(1,1)*u23(1) + fucell(2,1)*u23(2)+ fucell(3,2)*u23(3)

        frecip(1,:) = u23(:) / fbox_volume
        frecip(2,:) = u31(:) / fbox_volume
        frecip(3,:) = u12(:) / fbox_volume

        fbox_sphere = a*b*c
        do i = 1, 3
            dist = frecip(i,1)*fucell(1,i) + frecip(i,2)*fucell(2,i) + frecip(i,3)*fucell(3,i)
            dist = dist / sqrt( frecip(i,1)**2 + frecip(i,2)**2 + frecip(i,3)**2)
            if( dist .lt. fbox_sphere ) fbox_sphere = dist
        end do
        fbox_sphere = 0.5d0 * fbox_sphere
    end if

    p0VEne = fpressure * fbox_volume * PMF_PVCONV

end subroutine pmf_pbc_set_box_from_lvectors

!===============================================================================
! Subroutine:  pmf_pbc_get_cbox
! get box flag and center of box (cbox_x,cbox_y,cbox_z)
!===============================================================================

subroutine pmf_pbc_get_cbox(has_box,cbox)

    use pmf_dat
    use pmf_utils

    integer        :: has_box
    real(PMFDP)    :: cbox(3)
    ! --------------------------------------------------------------------------

    cbox = 0.0

    if( fbox_type .eq. BOX_ISOLATED_SYSTEM ) then
        has_box = 0
        return
    end if

    has_box    = 1
    cbox(:) = 0.5*(fucell(:,1) + fucell(:,2) + fucell(:,3))

end subroutine pmf_pbc_get_cbox

!===============================================================================
! Subroutine:  pmf_pbc_image_point
!===============================================================================

subroutine pmf_pbc_image_point(pos,familiar)

    use pmf_dat
    use pmf_utils
    use pmf_constants

    implicit none
    real(PMFDP)     :: pos(3)
    logical         :: familiar
    ! --------------------------------------------
    real(PMFDP)     :: fx,fy,fz
    real(PMFDP)     :: ffx,ffy,ffz
    real(PMFDP)     :: rfx,rfy,rfz
    real(PMFDP)     :: ofx,ofy,ofz
    real(PMFDP)     :: px,py,pz
    real(PMFDP)     :: cbox(3)
    integer         :: has_box
    integer         :: mx,my,mz
    integer         :: lx,ly,lz
    real(PMFDP)     :: min_dis, ds
    ! --------------------------------------------------------------------------

    select case(fbox_type)
        case(BOX_ISOLATED_SYSTEM)
            ! nothing to do
            return

        case(BOX_ORTHOGONAL)
            pos(1) = pos(1) - fucell(1,1)*floor(pos(1)*frecip(1,1))
            pos(2) = pos(2) - fucell(2,2)*floor(pos(2)*frecip(2,2))
            pos(3) = pos(3) - fucell(3,3)*floor(pos(3)*frecip(3,3))

        case(BOX_GENERAL)
            ! calculate fractional coordinates
            fx = pos(1)*frecip(1,1) + pos(2)*frecip(1,2) + pos(3)*frecip(1,3);
            fy = pos(1)*frecip(2,1) + pos(2)*frecip(2,2) + pos(3)*frecip(2,3);
            fz = pos(1)*frecip(3,1) + pos(2)*frecip(3,2) + pos(3)*frecip(3,3);

            ! in which are we?
            ffx = floor(fx);
            ffy = floor(fy);
            ffz = floor(fz);

            if( familiar ) then
                ! do familiar imaging
                rfx = fx - ffx;
                rfy = fy - ffy;
                rfz = fz - ffz;

                ! coordinates of cell centre
                call pmf_pbc_get_cbox(has_box,cbox)

                ! this is just upper guess of possible distance
                min_dis = 1e10;
                mx = 0;
                my = 0;
                mz = 0;

                do lx = -1, 1
                    do ly = -1, 1
                        do lz = -1, 1
                            ofx = rfx + lx
                            ofy = rfy + ly
                            ofz = rfz + lz

                            px = ofx*fucell(1,1) + ofy*fucell(1,2) + ofz*fucell(1,3)
                            py = ofx*fucell(2,1) + ofy*fucell(2,2) + ofz*fucell(2,3)
                            pz = ofx*fucell(3,1) + ofy*fucell(3,2) + ofz*fucell(3,3)

                            ds = (px-cbox(1))*(px-cbox(1)) &
                               + (py-cbox(2))*(py-cbox(2)) &
                               + (pz-cbox(3))*(pz-cbox(3))

                            if( ds .lt. min_dis ) then
                                mx = lx
                                my = ly
                                mz = lz
                                min_dis = ds
                            end if
                        end do
                    end do
                end do
                ffx = ffx - mx
                ffy = ffy - my
                ffz = ffz - mz
            end if

            ! image point
            pos(1) = pos(1) - ffx*fucell(1,1) - ffy*fucell(1,2) - ffz*fucell(1,3)
            pos(2) = pos(2) - ffx*fucell(2,1) - ffy*fucell(2,2) - ffz*fucell(2,3)
            pos(3) = pos(3) - ffx*fucell(3,1) - ffy*fucell(3,2) - ffz*fucell(3,3)

    end select

end subroutine pmf_pbc_image_point

!===============================================================================
! Subroutine:  pmf_pbc_image_crd
!===============================================================================

subroutine pmf_pbc_image_crd(orig,imaged)

    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: orig(3)
    real(PMFDP)    :: imaged(3)
    ! --------------------------------------------------------------------------

    imaged = orig
    call pmf_pbc_image_point(imaged,.true.)

end subroutine pmf_pbc_image_crd

!===============================================================================
! Subroutine:  pmf_pbc_image_vector
!===============================================================================

subroutine pmf_pbc_image_vector(vec)

    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)     :: vec(3)
    ! -----------------------------------------------
    real(PMFDP)     :: cbox(3)
    integer         :: has_box
    ! --------------------------------------------------------------------------

    if( fbox_type .eq. BOX_ISOLATED_SYSTEM ) then
        return
    end if

    call pmf_pbc_get_cbox(has_box,cbox)

    ! image vector
    vec = vec + cbox
    call pmf_pbc_image_point(vec,.true.)
    vec = vec - cbox

end subroutine pmf_pbc_image_vector

!===============================================================================
! Subroutine:  pmf_pbc_image_vector
!===============================================================================

subroutine pmf_pbc_image_vector3(vx,vy,vz)

    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: vx,vy,vz
    ! -----------------------------------------------
    real(PMFDP)    :: vec(3)
    ! --------------------------------------------------------------------------

    vec(1) = vx
    vec(2) = vy
    vec(3) = vz

    call pmf_pbc_image_vector(vec)

    vx = vec(1)
    vy = vec(2)
    vz = vec(3)

end subroutine pmf_pbc_image_vector3

!===============================================================================
! Subroutine:   pmf_pbc_write_xyz
!===============================================================================

subroutine pmf_pbc_write_xyz(coords, name, image)

    use pmf_constants
    use pmf_dat
    use smf_xyzfile
    use smf_xyzfile_type
    use smf_periodic_table_dat
    use smf_periodic_table

    implicit none
    real(PMFDP)         :: coords(:,:)  ! input coordinates
    character(*)        :: name         ! output file name
    logical             :: image        ! shall we image coordinates?
    ! --------------------------------------------
    type(XYZFILE_TYPE)  :: xyz
    integer             :: i
    ! --------------------------------------------------------------------------

    ! init xyz file
    call init_xyz(xyz)
    call allocate_xyz(xyz,NumOfLAtoms)

    ! copy coordinates and symbols
    xyz%cvs = coords
    do i=1,NumOfLAtoms
        xyz%symbols(i) = pt_symbols(SearchZByMass(Mass(i),0.5d0))
    end do

    ! image coordinates if necessary
    if( image .eqv. .true. ) then
        do i=1,NumOfLAtoms
            call pmf_pbc_image_point(xyz%cvs(:,i),.true.)
        end do
    end if

    ! write xyz file and close everything
    call open_xyz(PMF_XYZ,name,xyz,'UNKNOWN')
    call write_xyz(PMF_XYZ,xyz)
    call close_xyz(PMF_XYZ,xyz)
    call free_xyz(xyz)

end subroutine pmf_pbc_write_xyz

!===============================================================================

end module pmf_pbc
