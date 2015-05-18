! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License along
!     with this program; if not, write to the Free Software Foundation, Inc.,
!     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ==============================================================================

module pmfdyn_restraints

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine init_restraints

    use pmfdyn_restraints_dat

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,10)
    write(PMF_OUT,20) nrstr_seq

    10 format(/,'Initializing restraints ...')
    20 format(  '   Number of sequence restraints           = ',i10)

end subroutine init_restraints

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine potene_restraints(x,d,ene)

    use smf_profiling
    use pmfdyn_restraints_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)     :: x(:,:)
    real(PMFDP)     :: d(:,:)
    real(PMFDP)     :: ene
    ! --------------------------------------------------------------------------

    call start_timer(RE_TIMER)

    call potene_seqrestraints(x,d,ene)

    call stop_timer(RE_TIMER)

    return

end subroutine potene_restraints

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine potene_seqrestraints(x,d,ene)

    use pmfdyn_restraints_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)        :: x(:,:)
    real(PMFDP)        :: d(:,:)
    real(PMFDP)        :: ene
    ! -------------------------------------
    integer            :: ir, i, n_ctr
    real(PMFDP)        :: fk, dr(3), r2, edum, totmass
    ! -----------------------------------------------------------------------------

    ! sequence restraints (independent of Q-state)
    do ir = 1, nrstr_seq

        ! get real fk
        if( rstseq(ir)%change == 0 ) then
            fk = rstseq(ir)%fk
        else
            fk = rstseq(ir)%fk + (rstseq(ir)%fk_final - rstseq(ir)%fk)*istep/nsteps
        end if

        select case (rstseq(ir)%to_centre)
            case (SEQRES_FIX_TO_XYZ)
                ! restrain each atom to its topology co-ordinate
                do i = rstseq(ir)%i, rstseq(ir)%j
                    if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
                        dr(:)   = x(:,i) - xtop(:,i)
                        r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
                        edum    = 0.5*fk*r2
                        ene  = ene + edum
                        d(:,i) = d(:,i) + fk*dr(:)
                    end if
                end do

            case (SEQRES_FIX_TO_COG)
                ! restrain to geometrical centre
                ! reset dr & atom counter
                dr(:) = 0.
                n_ctr = 0
                ! calculate deviation from center
                do i = rstseq(ir)%i, rstseq(ir)%j
                    if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
                        n_ctr = n_ctr + 1
                        dr(:) = dr(:) + x(:,i) - xtop(:,i)
                    end if
                end do
                if(n_ctr > 0) then
                    ! only if atoms were found:
                    ! form average
                    dr(:) = dr(:) / n_ctr
                    r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
                    edum    = 0.5*fk*r2
                    ene  = ene + edum
                    ! apply same force to all atoms
                    do i = rstseq(ir)%i, rstseq(ir)%j
                        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
                                d(:,i) = d(:,i) + fk*dr(:)
                        end if
                    end do
                end if

            case (SEQRES_FIX_TO_COM)
                ! restrain to mass centre
                ! reset dr & variable to put masses
                dr(:) = 0.
                totmass = 0.

                ! calculate deviation from mass center
                do i = rstseq(ir)%i, rstseq(ir)%j
                    if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
                        totmass = totmass + mass(i)                              ! Add masses
                        dr(:) = dr(:) + (x(:,i) - xtop(:,i))*mass(i) ! Massweight distances
                    end if
                end do

                if(totmass > 0.0d0) then
                    ! only if atoms were found: (i.e has a total mass)
                    ! from average
                    dr(:) = dr(:)/totmass                                  ! divide by total mass
                    r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
                    edum    = 0.5*fk*r2
                    ene  = ene + edum
                    ! apply same force to all atoms
                    do i = rstseq(ir)%i, rstseq(ir)%j
                        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
                                d(:,i) = d(:,i) + fk*dr(:)*mass(i)/totmass
                        end if
                    end do
                end if
        end select

    end do

    return

end subroutine potene_seqrestraints

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmfdyn_restraints
