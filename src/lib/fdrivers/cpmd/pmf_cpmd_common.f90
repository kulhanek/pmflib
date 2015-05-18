!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
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

module pmf_cpmd_common

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  remap_masses_in
!===============================================================================

subroutine remap_masses_in(NSP,NSX,NA,PMA0,fmass)

    use pmf_dat

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: PMA0(NSX)            ! masses
    real(PMFDP)    :: fmass(:)             ! local masses
    ! -----------------------------------------------
    integer        :: iat,is,ia
    ! --------------------------------------------------------------------------

    iat=0
    do is=1,NSP
        do ia=1,NA(IS)
            iat=iat+1
            fmass(iat) = PMA0(is)
        end do
    end do

end subroutine remap_masses_in

!===============================================================================
! Subroutine:  remap_qmmm_masses_in
!===============================================================================

subroutine remap_qmmm_masses_in(NTAT,NSAT,NDAT,NSX,CPSP,PMA0,fmass)

    use pmf_dat

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NSX                  ! max number of atom types
    integer        :: CPSP(NTAT)           ! atom types
    real(PMFDP)    :: PMA0(NSX)            ! masses
    real(PMFDP)    :: fmass(:)             ! local masses
    ! -----------------------------------------------
    integer        :: iat,lat
    ! --------------------------------------------------------------------------

    lat = 1
    do iat=1,NTAT
        if( (iat .gt. NSAT) .and. (iat .le. NSAT+NDAT) ) cycle ! skip dummy atoms
        fmass(lat) = PMA0(CPSP(iat))
        lat = lat + 1
    end do

end subroutine remap_qmmm_masses_in

!===============================================================================
! Subroutine:  remap_coords_in
!===============================================================================

subroutine remap_coords_in(NSP,NAX,NSX,NA,TAUP,lx)

    use pmf_dat
    use pmf_constants

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates
    real(PMFDP)    :: lx(:,:)              ! local coordinates
    ! -----------------------------------------------
    integer        :: iat,is,ia
    ! --------------------------------------------------------------------------

    iat=0
    do is=1,NSP
        do ia=1,NA(IS)
            iat=iat+1
            lx(:,iat) = TAUP(:,ia,is)
        end do
    end do

end subroutine remap_coords_in

!===============================================================================
! Subroutine:  remap_qmmm_coords_in
!===============================================================================

subroutine remap_qmmm_coords_in(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAUP,lx)

    use pmf_dat
    use pmf_constants

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)           ! atom indexes
    integer        :: CPSP(NTAT)           ! atom types
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates
    real(PMFDP)    :: lx(:,:)              ! local coordinates
    ! -----------------------------------------------
    integer        :: iat,lat
    ! --------------------------------------------------------------------------

    lat = 1
    do iat=1,NTAT
        if( (iat .gt. NSAT) .and. (iat .le. NSAT+NDAT) ) cycle ! skip dummy atoms
        lx(:,lat) = TAUP(:,CPAT(iat),CPSP(iat))
        lat = lat + 1
    end do

end subroutine remap_qmmm_coords_in

!===============================================================================
! Subroutine:  remap_coords_back
!===============================================================================

subroutine remap_coords_back(NSP,NAX,NSX,NA,TAUP,lxp)

    use pmf_dat
    use pmf_constants

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates
    real(PMFDP)    :: lxp(:,:)              ! local coordinates
    ! -----------------------------------------------
    integer        :: iat,is,ia
    ! --------------------------------------------------------------------------

    iat=0
    do is=1,NSP
        do ia=1,NA(IS)
            iat=iat+1
            TAUP(:,ia,is) = lxp(:,iat)
        end do
    end do

end subroutine remap_coords_back

!===============================================================================
! Subroutine:  remap_qmmm_coords_back
!===============================================================================

subroutine remap_qmmm_coords_back(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAUP,lxp)

    use pmf_dat
    use pmf_constants

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)            ! atom indexes
    integer        :: CPSP(NTAT)            ! atom types
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates
    real(PMFDP)    :: lxp(:,:)             ! local coordinates
    ! -----------------------------------------------
    integer        :: iat,lat
    ! --------------------------------------------------------------------------

    lat = 1
    do iat=1,NTAT
        if( (iat .gt. NSAT) .and. (iat .le. NSAT+NDAT) ) cycle ! skip dummy atoms
        TAUP(:,CPAT(iat),CPSP(iat)) = lxp(:,lat)
        lat = lat + 1
    end do

end subroutine remap_qmmm_coords_back

!===============================================================================
! Subroutine:  remap_force_in_add
!===============================================================================

subroutine remap_force_in_add(NSP,NAX,NSX,NA,FION,lf)

    use pmf_dat

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: FION(3,NAX,NSX)      ! forces
    real(PMFDP)    :: lf(:,:)              ! local forces
    ! -----------------------------------------------
    integer        :: iat,is,ia
    ! --------------------------------------------------------------------------

    iat=0
    do is=1,NSP
        do ia=1,NA(IS)
            iat=iat+1
            lf(:,iat) = lf(:,iat) + FION(:,ia,is)
        end do
    end do

end subroutine remap_force_in_add

!===============================================================================
! Subroutine:  remap_qmmm_force_in_add
!===============================================================================

subroutine remap_qmmm_force_in_add(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,FION,lf)

    use pmf_dat

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)            ! atom indexes
    integer        :: CPSP(NTAT)            ! atom types
    real(PMFDP)    :: FION(3,NAX,NSX)      ! forces
    real(PMFDP)    :: lf(:,:)              ! local forces
    ! -----------------------------------------------
    integer        :: iat,lat
    ! --------------------------------------------------------------------------

    lat = 1
    do iat=1,NTAT
        if( (iat .gt. NSAT) .and. (iat .le. NSAT+NDAT) ) cycle ! skip dummy atoms
        lf(:,lat) = lf(:,lat) + FION(:,CPAT(iat),CPSP(iat))
        lat = lat + 1
    end do

end subroutine remap_qmmm_force_in_add

!===============================================================================
! Subroutine:  remap_force_back
!===============================================================================

subroutine remap_force_back(NSP,NAX,NSX,NA,FION,lf)

    use pmf_dat

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: FION(3,NAX,NSX)      ! forces
    real(PMFDP)    :: lf(:,:)              ! local forces
    ! -----------------------------------------------
    integer        :: iat,is,ia
    ! --------------------------------------------------------------------------

    iat=0
    do is=1,NSP
        do ia=1,NA(IS)
            iat=iat+1
            FION(:,ia,is) = lf(:,iat)
        end do
    end do

end subroutine remap_force_back

!===============================================================================
! Subroutine:  remap_force_back
!===============================================================================

subroutine remap_qmmm_force_back(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,FION,lf)

    use pmf_dat

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)            ! atom indexes
    integer        :: CPSP(NTAT)            ! atom types
    real(PMFDP)    :: FION(3,NAX,NSX)      ! forces
    real(PMFDP)    :: lf(:,:)              ! local forces
    ! -----------------------------------------------
    integer        :: iat,lat
    ! --------------------------------------------------------------------------

    lat = 1
    do iat=1,NTAT
        if( (iat .gt. NSAT) .and. (iat .le. NSAT+NDAT) ) cycle ! skip dummy atoms
        FION(:,CPAT(iat),CPSP(iat)) = lf(:,lat)
        lat = lat + 1
    end do

end subroutine remap_qmmm_force_back

!===============================================================================

end module pmf_cpmd_common
