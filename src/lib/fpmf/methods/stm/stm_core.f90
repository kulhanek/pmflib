!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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
!    Lesser General Public License for more detajls.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module stm_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  stm_core_main
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine stm_core_main

    use stm_output
    use stm_client
    use stm_dat

    ! --------------------------------------------------------------------------

    call stm_client_exchange_data
    call stm_core_force
    call stm_output_write

end subroutine stm_core_main

!===============================================================================
! Subroutine:  stm_core_force
! this is leap-frog and velocity-verlet STM version
!===============================================================================

subroutine stm_core_force

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use stm_dat
    use stm_cvs

    implicit none
    integer         :: i,j,k,l
    real(PMFDP)     :: rvalue
    ! --------------------------------------------------------------------------

    curstep = curstep + 1

    ! process increments ----------------------------
    do i=1,NumOfSTMCVs
        call stm_restraints_increment(STMCVList(i))
    end do

    ! calculate energy and gradients ----------------
    TotalSTMEnergy = 0.0
    do i=1,NumOfSTMCVs
        rvalue = CVContext%CVsValues(STMCVList(i)%cvindx)
        STMCVList(i)%deviation = STMCVList(i)%cv%get_deviation(rvalue,STMCVList(i)%target_value)
        STMCVList(i)%energy = 0.5d0*STMCVList(i)%force_constant*STMCVList(i)%deviation**2
        TotalSTMEnergy = TotalSTMEnergy + STMCVList(i)%energy
        ! correct forces -----------------------------
        Frc(:,:) = Frc(:,:) - STMCVList(i)%force_constant*STMCVList(i)%deviation*CVContext%CVsDrvs(:,:,STMCVList(i)%cvindx)
        ! get PMF
        PMF(i) = PMF(i) - STMCVList(i)%force_constant*STMCVList(i)%deviation
    end do

    ! accumulate MTZ
    select case(ftensor)
        case (0)
            do i=1,NumOfSTMCVs
                MTZ(i,i) = MTZ(i,i) + 1.0
            end do
        case (1)
            do i=1,NumOfSTMCVs
                do j=1,NumOfSTMCVs
                    do k=1,NumOfLAtoms
                        do l=1,3
                            MTZ(j,i) = MTZ(j,i) + &
                                       CVContext%CVsDrvs(l,k,STMCVList(i)%cvindx) &
                                      *CVContext%CVsDrvs(l,k,STMCVList(j)%cvindx)
                        end do
                    end do
                end do
            end do
        case (2)
            do i=1,NumOfSTMCVs
                do j=1,NumOfSTMCVs
                    do k=1,NumOfLAtoms
                        do l=1,3
                            MTZ(j,i) = MTZ(j,i) + &
                                       MassInv(k)*CVContext%CVsDrvs(l,k,STMCVList(i)%cvindx) &
                                                 *CVContext%CVsDrvs(l,k,STMCVList(j)%cvindx)
                        end do
                    end do
                end do
            end do
    end select

    return

end subroutine stm_core_force

!===============================================================================

end module stm_core
