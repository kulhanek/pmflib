!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtd_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_core_main
!===============================================================================

subroutine mtd_core_main

 use mtd_output
 use mtd_client
 ! ------------------------------------------------------------------------------------------

 call mtd_core_force
 call mtd_output_write
 call mtd_gp_output
 call mtd_client_exchange_data

end subroutine mtd_core_main

!===============================================================================
! Subroutine:  mtd_core_force
!===============================================================================

subroutine mtd_core_force

 use pmf_utils
 use mtd_dat

 implicit none
 ! -----------------------------------------------------------------------------

 select case(fmode)
    case(0)
        return
    case(1)
        call mtd_core_direct()
    case(2)
        call mtd_core_grid()
    case(3)
        call mtd_core_gpmtd()
    case default
        call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported fmode!')
 end select

end subroutine mtd_core_force

!===============================================================================
! Function:  mtd_core_is_value_in_range
!===============================================================================

logical function mtd_core_is_value_in_range(values)

    use mtd_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: values(:)
    ! -----------------------------------------------
    integer                :: i
    real(PMFDP)            :: val 
    ! -----------------------------------------------

    mtd_core_is_value_in_range = .true.

    do i=1,NumOfMTDCVs
       if( .not. MTDCVList(i)%cv%is_periodic_cv() ) then
           if( (values(i) < MTDCVList(i)%min_value) .or. (values(i) > MTDCVList(i)%max_value) ) then
               mtd_core_is_value_in_range = .false.
               return
           end if
       else
           ! move value into default range
           val = values(i) - MTDCVList(i)%cv%get_period_cv_value() / 2.0d0 - MTDCVList(i)%cv%get_min_cv_value()
           val = val - nint(values(i)/MTDCVList(i)%cv%get_period_cv_value())*MTDCVList(i)%cv%get_period_cv_value()
           val = val + MTDCVList(i)%cv%get_period_cv_value() / 2.0d0 + MTDCVList(i)%cv%get_min_cv_value()
           ! min_value < max_value
           if( MTDCVList(i)%min_value < MTDCVList(i)%max_value ) then
               if( (val < MTDCVList(i)%min_value) .or. (val > MTDCVList(i)%max_value) ) then
                   mtd_core_is_value_in_range = .false.
                   return
               end if
           else
           ! min_value > max_value
               if( (val < MTDCVList(i)%min_value) .and. (val > MTDCVList(i)%max_value) ) then
                   mtd_core_is_value_in_range = .false.
                   return
               end if
           end if
       end if
    end do

    return

end function mtd_core_is_value_in_range

!===============================================================================
! Subroutine:  mtd_core_direct
! direct metadynamics algorithm
!===============================================================================

subroutine mtd_core_direct

 use pmf_dat
 use mtd_dat
 use mtd_history

 implicit none
 integer                    :: i,j,k,ci
 type(MTDHistType),pointer  :: current_history
 real(PMFDP)                :: fexparg,mf,fh,diff
 real(PMFDP)                :: values(NumOfMTDCVs)
 real(PMFDP)                :: hfactor
 logical                    :: skip
 ! -----------------------------------------------------------------------------

!!do we need to update history list, eg. to add new hill?
!if( fmetastep .gt. 0 ) then
!   if( mod(fstep,fmetastep) .eq. 0 ) then
!       ! add new hill
!       meta_step = meta_step + 1
!       call mtd_history_add_new_hill
!   end if
!end if

 TotalMTDEnergy = 0.0d0

 ! ==========================================================
 ! now calculate meta forces and total history potential value
!do i = 1,NumOfMTDCVs
!   mf = 0.0d0
!   ! energy will be recalculated several times :-) maybe some further optimization of code?
!   TotalMTDEnergy = 0.0d0
!   ! go through history
!   current_history => hill_history
!   do while( associated(current_history) )
!       ! for every item in buffer
!       do j = 1,current_history%nrst_of_values
!           fexparg = 0.0d0
!           do k = 1,NumOfMTDCVs
!               diff = MTDCVList(k)%cv%get_deviation(CVContext%CVsValues(MTDCVList(k)%cvindx),current_history%values(j,k))
!               fexparg = fexparg + diff**2 / (2.0d0 * current_history%widths(j,k)**2)
!           end do
!           fh = current_history%heights(j)*exp(-fexparg)
!           diff = MTDCVList(i)%cv%get_deviation(CVContext%CVsValues(MTDCVList(i)%cvindx),current_history%values(j,i))
!           mf = mf + fh * diff / current_history%widths(j,i)**2
!           TotalMTDEnergy = TotalMTDEnergy + fh
!       end do
!       current_history => current_history%next_history_buffer
!   end do
!   MTDCVList(i)%meta_force = mf
!end do

 ! here comes the wished optimization
 ! clear  MTDCVList(:)%meta_force array
 do i = 1,NumOfMTDCVs
    MTDCVList(i)%meta_force = 0.0d0
 end do

 ! go through history
 current_history => hill_history
 do while( associated(current_history) )
    ! for every item in buffer
    do j = 1,current_history%nrst_of_values
       call mtd_core_hill_scaling_factor(current_history%steps(j),hfactor)
       if( hfactor .eq. 0.0d0 ) cycle
       skip = .false.
       fexparg = 0.0d0
       do i = 1,NumOfMTDCVs
           diff = MTDCVList(i)%cv%get_deviation(CVContext%CVsValues(MTDCVList(i)%cvindx),current_history%values(j,i))
           ! skip this hill if it is further than specified
           if( diff .gt. MTDCVList(i)%max_dist ) then
               skip = .true.
               exit
           end if
           fexparg = fexparg + diff**2 / (2.0d0 * current_history%widths(j,i)**2)
       end do
       ! skip this hill and the rest of the calculation related to it 
       if( skip ) exit
       fh = hfactor*current_history%heights(j)*exp(-fexparg)
       TotalMTDEnergy = TotalMTDEnergy + fh
       do i = 1,NumOfMTDCVs
           diff = MTDCVList(i)%cv%get_deviation(CVContext%CVsValues(MTDCVList(i)%cvindx),current_history%values(j,i))
           MTDCVList(i)%meta_force = MTDCVList(i)%meta_force + fh * diff / current_history%widths(j,i)**2
       end do
    end do
    current_history => current_history%next_history_buffer
 end do

 ! lam81: well-tempered MTD is implemented
 ! do we need to update history list, eg. to add new hill?
 if( fmetastep .gt. 0 ) then
    if( (fmetatemp .eq. 0) .or. (fmetavary .eq. 1) ) then
        if( mod(fstep,fmetastep) .eq. 0 ) then
            ! add new hill
            meta_step = meta_step + 1
            call mtd_history_add_new_hill(TotalMTDEnergy)
        end if
    else
        if( fstep .eq. meta_next_fstep ) then
            ! add new hill
            meta_step = meta_step + 1
            call mtd_history_add_new_hill(TotalMTDEnergy)
        end if
    end if
 end if

 ! check if values are in the range

 do i = 1,NumOfMTDCVs
    ci = MTDCVList(i)%cvindx
    values(i) = CVContext%CVsValues(ci)
 end do

 if ( mtd_core_is_value_in_range(values) ) then

   ! now update forces =======================================
   do i = 1,NumOfMTDCVs
      ci = MTDCVList(i)%cvindx
      mf = MTDCVList(i)%meta_force
      do j = 1,NumOfLAtoms
          Frc(:,j) = Frc(:,j) + mf*CVContext%CVsDrvs(:,j,ci)
      end do
      CVContext%CVsFrc(ci) = CVContext%CVsFrc(ci) + mf  ! lam81
   end do

 else

   TotalMTDEnergy = 0.0d0

 end if

!  write(*,*) 'Analytical='
!  write(*,*) fdx(:,:)
!  call numerical_derivatives
!  write(*,*) 'Numerical='
!  write(*,*) fdx(:,:)

return

end subroutine mtd_core_direct

!===============================================================================
! Subroutine:  mtd_core_grid
! grid metadynamics algorithm
!===============================================================================

subroutine mtd_core_grid

 use pmf_dat
 use mtd_dat
 use mtd_history
 use mtd_grid

 implicit none
 integer                    :: i,j,ci
 real(PMFDP)                :: values(NumOfMTDCVs)
 real(PMFDP)                :: mf
 ! -----------------------------------------------------------------------------

 TotalMTDEnergy = 0.0d0

 ! clear  MTDCVList(:)%meta_force array
 do i = 1,NumOfMTDCVs
    MTDCVList(i)%meta_force = 0.0d0
 end do

 do i = 1,NumOfMTDCVs
    ci = MTDCVList(i)%cvindx
    values(i) = CVContext%CVsValues(ci)
 end do

 call mtd_grid_get_data(values, TotalMTDEnergy, MTDCVList(:)%meta_force)

 ! lam81: well-tempered MTD is implemented
 ! do we need to update history list, eg. to add new hill?
 if( fmetastep .gt. 0 ) then
    if( (fmetatemp .eq. 0) .or. (fmetavary .eq. 1) ) then
        if( mod(fstep,fmetastep) .eq. 0 ) then
            ! add new hill
            meta_step = meta_step + 1
            call mtd_history_add_new_hill(TotalMTDEnergy)
        end if
    else
        if( fstep .eq. meta_next_fstep ) then
            ! add new hill
            meta_step = meta_step + 1
            call mtd_history_add_new_hill(TotalMTDEnergy)
        end if
    end if
 end if

 ! now update forces =======================================
 do i = 1,NumOfMTDCVs
    ci = MTDCVList(i)%cvindx
    mf = MTDCVList(i)%meta_force
    do j = 1,NumOfLAtoms
        Frc(:,j) = Frc(:,j) + mf*CVContext%CVsDrvs(:,j,ci)
    end do
    CVContext%CVsFrc(ci) = CVContext%CVsFrc(ci) + mf  ! lam81
 end do

end subroutine mtd_core_grid

!===============================================================================
! Subroutine:  mtd_core_gpmtd
! gp metadynamics algorithm
!===============================================================================

subroutine mtd_core_gpmtd

 use pmf_dat
 use mtd_dat
 use mtd_history
 use clustering_mod
 use gp_dat_mod
 use gp_basic_mod

 implicit none
 integer                    :: i,j,k,ci
 type(MTDHistType),pointer  :: current_history
 real(PMFDP)                :: fexparg,mf,fh,diff
 real(PMFDP)                :: values(NumOfMTDCVs)
 ! -----------------------------------------------------------------------------

 TotalMTDEnergy = 0.0d0

 ! clear  MTDCVList(:)%meta_force array
 do i = 1,NumOfMTDCVs
    MTDCVList(i)%meta_force = 0.0d0
 end do

 do i = 1,NumOfMTDCVs
    ci = MTDCVList(i)%cvindx
    values(i) = CVContext%CVsValues(ci)
 end do

 ! go through available history
 current_history => hill_history
 do while( associated(current_history) )
    ! for every item in buffer
    do j = 1,current_history%nrst_of_values
       fexparg = 0.0d0
       do i = 1,NumOfMTDCVs
           diff = MTDCVList(i)%cv%get_deviation(CVContext%CVsValues(MTDCVList(i)%cvindx),current_history%values(j,i))
           fexparg = fexparg + diff**2 / (2.0d0 * current_history%widths(j,i)**2)
       end do
       fh = current_history%heights(j)*exp(-fexparg)
       TotalMTDEnergy = TotalMTDEnergy + fh
       do i = 1,NumOfMTDCVs
           diff = MTDCVList(i)%cv%get_deviation(CVContext%CVsValues(MTDCVList(i)%cvindx),current_history%values(j,i))
           MTDCVList(i)%meta_force = MTDCVList(i)%meta_force + fh * diff / current_history%widths(j,i)**2
       end do
    end do
    current_history => current_history%next_history_buffer
 end do

 call mtd_core_get_data_gp(values, TotalMTDEnergy, MTDCVList(:)%meta_force)

 ! lam81: well-tempered MTD is implemented
 ! do we need to update history list, eg. to add new hill?
 if( fmetastep .gt. 0 ) then
    if( (fmetatemp .eq. 0) .or. (fmetavary .eq. 1) ) then
        if( mod(fstep,fmetastep) .eq. 0 ) then
            ! add new hill
            meta_step = meta_step + 1
            call mtd_history_add_new_hill(TotalMTDEnergy)
        end if
    else
        if( fstep .eq. meta_next_fstep ) then
            ! add new hill
            meta_step = meta_step + 1
            call mtd_history_add_new_hill(TotalMTDEnergy)
        end if
    end if
 end if

 ! now update forces =======================================
 do i = 1,NumOfMTDCVs
    ci = MTDCVList(i)%cvindx
    mf = MTDCVList(i)%meta_force
    do j = 1,NumOfLAtoms
        Frc(:,j) = Frc(:,j) + mf*CVContext%CVsDrvs(:,j,ci)
    end do
    CVContext%CVsFrc(ci) = CVContext%CVsFrc(ci) + mf  ! lam81
 end do

 call mtd_core_gp_sparsification(values)

 call mtd_core_gp_data()

end subroutine mtd_core_gpmtd

!===============================================================================
! Subroutine:  mtd_core_get_data_gp
!===============================================================================

subroutine mtd_core_get_data_gp(values, mtdene, mtdforce)

 use pmf_dat
 use mtd_dat
 use gp_dat_mod
 use gp_basic_mod

 implicit none
 real(PMFDP), intent(in)       :: values(:)
 real(PMFDP), intent(inout)    :: mtdene
 real(PMFDP), intent(inout)    :: mtdforce(:)
 ! -----------------------------------------------

 ! return until we do not have any estimate
 if( .not. gpmtd%gp%initialised ) return

 ! get energy from GP prediction
 mtdene = mtdene + f_predict(gpmtd%gp,values,SE_kernel_r_rr)

 ! get forces from GP prediction
 mtdforce = mtdforce - f_predict_grad(gpmtd%gp,values,SE_kernel_r_rr)

end subroutine mtd_core_get_data_gp

!===============================================================================
! Subroutine:  mtd_core_gp_sparsification
!===============================================================================

subroutine mtd_core_gp_sparsification(values)

 use pmf_dat
 use pmf_utils
 use mtd_dat

 implicit none
 real(PMFDP)    :: values(:)

 select case ( fgpsparsification )
    case(1)
      call mtd_core_gp_sparsification_kmeans(values)
    case(2)
      call mtd_core_gp_sparsification_grid(values)
 end select

end subroutine mtd_core_gp_sparsification

!===============================================================================
! Subroutine:  mtd_core_gp_sparsification_kmeans
!===============================================================================

subroutine mtd_core_gp_sparsification_kmeans(values)

    use mtd_dat
    use pmf_dat
    use clustering_mod

    implicit none
    real(PMFDP)    :: values(:)

    integer :: i
    integer, allocatable :: ppsparseindex(:)
    integer, allocatable :: sparseindex(:)
    ! --------------------------------------------------------------------------

    ! first: add value to the gpmtd%clusterarray list
    if ( mod(fstep,fgpclusterfreq) .eq. 0 ) then
       fgpncluster = fgpncluster + 1
       gpmtd%clusterarray(:,fgpncluster) = values(:)
       if(fgpncluster .eq. fgpsparse) gpmtd%sparseindices(1:fgpncluster) = (/(i, i=1, fgpsparse)/)
    end if

    ! second: do kmeans calculation to find sparse points
    if ( (fgpncluster .le. fgpsparse) .and. (mod(fstep,fgpclusterfreq) .eq. 0) ) then
       fgpsparsechange = .true.
       if (associated(gpmtd%sparsepoints)) nullify(gpmtd%sparsepoints)
       allocate(gpmtd%sparsepoints(NumOfMTDCVs,fgpncluster))
       gpmtd%sparsepoints = gpmtd%clusterarray(:,1:fgpncluster)
    else if ( mod(fstep,fgpsparsefreq) .eq. 0 ) then
       fgpsparsechange = .true.
       allocate(sparseindex(fgpsparse))
       call k_meanspp_clustering_pick(x=gpmtd%clusterarray(:,1:fgpncluster), periodicity=gpmtd%periodicities, &
                                      cluster_indices=sparseindex,preselected_indices=gpmtd%sparseindices)
       call k_means_clustering_pick(x=gpmtd%clusterarray(:,1:fgpncluster),periodicity=gpmtd%periodicities, &
                                    cluster_indices=gpmtd%sparseindices,starting_indices=sparseindex)
       if (associated(gpmtd%sparsepoints)) nullify(gpmtd%sparsepoints)
       allocate(gpmtd%sparsepoints(NumOfMTDCVs,fgpsparse))
       do i=1, fgpsparse
          gpmtd%sparsepoints(:,i) = gpmtd%clusterarray(:,gpmtd%sparseindices(i))
       end do
       deallocate(sparseindex)
    end if

end subroutine mtd_core_gp_sparsification_kmeans

!===============================================================================
! Subroutine:  mtd_core_gp_sparsification_grid
!===============================================================================

subroutine mtd_core_gp_sparsification_grid(values)

    use mtd_dat
    use pmf_dat

    implicit none
    real(PMFDP)    :: values(:)

    logical :: gridchange
    integer :: i, j, k
    real(PMFDP)    :: val, diffmin, diffmax
    integer        :: multibins
    real(PMFDP)    :: width
    ! --------------------------------------------------------------------------

    ! first: check if grid is changed
    gridchange = .false.

    ! if it is the first sample than its values will specify the first min and max positions as well
    if ( fstep .eq. 1 ) then
       gridchange = .true.
       do i=1,NumOfMTDCVs
          if ( gpmtd%fixgrid(i) .eq. 1 ) then
               gpmtd%min_positions(i) = MTDCVList(i)%min_value
               gpmtd%max_positions(i) = MTDCVList(i)%max_value
               gpmtd%range_positions(i) = gpmtd%max_positions(i) - gpmtd%min_positions(i)
               if ( gpmtd%range_positions(i) < 0 ) then
                    gpmtd%range_positions(i) = gpmtd%range_positions(i) + gpmtd%periodicities(i)
               end if
          else
               gpmtd%min_positions(i) = values(i)
               gpmtd%max_positions(i) = values(i)
               gpmtd%range_positions(i) = 0.0d0
          end if
       end do
    ! otherwise update min and max positions according to the new values
    else
       do i=1,NumOfMTDCVs
          if ( gpmtd%fixgrid(i) .eq. 1 ) cycle
          ! cv is not periodic
          if ( gpmtd%periodicities(i) .eq. 0 ) then
               if ( values(i) < gpmtd%min_positions(i) ) then
                    gpmtd%min_positions(i) = values(i)
                    gridchange = .true.
                    gpmtd%range_positions(i) = gpmtd%max_positions(i) - gpmtd%min_positions(i)
               else if ( values(i) > gpmtd%max_positions(i) ) then
                    gpmtd%max_positions(i) = values(i)
                    gridchange = .true.
                    gpmtd%range_positions(i) = gpmtd%max_positions(i) - gpmtd%min_positions(i)
               end if
          ! cv is periodic
          else if ( gpmtd%range_positions(i) .ne. gpmtd%periodicities(i) ) then
               ! move value into default range: useful when range is e.g. [0:2*pi]]
               val = values(i) - MTDCVList(i)%cv%get_period_cv_value() / 2.0d0 - MTDCVList(i)%cv%get_min_cv_value()
               val = val - nint(values(i)/MTDCVList(i)%cv%get_period_cv_value())*MTDCVList(i)%cv%get_period_cv_value()
               val = val + MTDCVList(i)%cv%get_period_cv_value() / 2.0d0 + MTDCVList(i)%cv%get_min_cv_value()
               diffmin = abs(val - gpmtd%min_positions(i))
               diffmax = abs(val - gpmtd%max_positions(i))
               ! criterion of new min/max
               ! min < max
               if ( gpmtd%min_positions(i) < gpmtd%max_positions(i) ) then
                    if( (val < gpmtd%min_positions(i)) .or. (val > gpmtd%max_positions(i)) ) then
                        if( diffmin < diffmax ) then
                            gpmtd%min_positions(i) = val
                            gridchange = .true.
                        else
                            gpmtd%max_positions(i) = val
                            gridchange = .true.
                        end if
                    end if
               ! min > max
               else if ( gpmtd%min_positions(i) > gpmtd%max_positions(i) ) then
                    if( (val < gpmtd%min_positions(i)) .and. (val > gpmtd%max_positions(i)) ) then
                        if( diffmin < diffmax ) then
                            gpmtd%min_positions(i) = val
                            gridchange = .true.
                        else
                            gpmtd%max_positions(i) = val
                            gridchange = .true.
                        end if
                    end if
               ! min = max
               else
                    ! min/max and val are close to each other
                    if( diffmin .le. (gpmtd%periodicities(i) / 2.0d0) ) then
                        if( val < gpmtd%min_positions(i) ) then
                            gpmtd%min_positions(i) = val
                            gridchange = .true.
                        else
                            gpmtd%max_positions(i) = val
                            gridchange = .true.
                        end if
                    ! min/max and val are far from each other 
                    else
                        ! min/max is close to right edge, val appears on left
                        if( val < gpmtd%min_positions(i) ) then
                            gpmtd%max_positions(i) = val
                            gridchange = .true.
                        ! min/max is close to left edge, val appears on right
                        else
                            gpmtd%min_positions(i) = val
                            gridchange = .true.
                        end if
                    end if
               end if
               gpmtd%range_positions(i) = gpmtd%max_positions(i) - gpmtd%min_positions(i)
               if ( gpmtd%range_positions(i) .lt. 0 ) then
                    gpmtd%range_positions(i) = gpmtd%range_positions(i) + gpmtd%periodicities(i)
               end if
          end if
       end do
    end if

    ! calculate new sparse points
    if ( gridchange ) then
       do i=1,NumOfMTDCVs

          width = gpmtd%range_positions(i) / real( gpmtd%numgrid(i) )

          multibins = 1

          do j=1,i-1
              multibins = multibins * gpmtd%numgrid(j)
          end do

          do j=1,fgpsparse
             k = mod(((j-1)/multibins),gpmtd%numgrid(i))
             gpmtd%sparsepoints(i,j) = gpmtd%min_positions(i) + real(k) * width + width / 2.0d0
             if ( gpmtd%periodicities(i) .ne. 0 ) then
                  gpmtd%sparsepoints(i,j) = gpmtd%sparsepoints(i,j) &
                                       - nint(gpmtd%sparsepoints(i,j)/gpmtd%periodicities(i))*gpmtd%periodicities(i)
             end if
          end do
       end do
       fgpsparsechange = .true.
    end if

end subroutine mtd_core_gp_sparsification_grid

!===============================================================================
! Subroutine:  mtd_core_gp_data
!===============================================================================

subroutine mtd_core_gp_data

    use mtd_dat
    use pmf_dat
    use gp_dat_mod
    use gp_basic_mod

    implicit none
    integer        :: n, i, j
    real(PMFDP), allocatable   :: mtd_energy(:), mtd_var(:)
    type(MTDHistType),pointer  :: current_history
    real(PMFDP)                :: fexparg,fh,diff
    ! --------------------------------------------------------------------------

    if ( mod(fstep,fgpcovfreq) .ne. 0 ) return

    if( hill_history%nrst_of_values .eq. 0 ) return

    allocate(mtd_energy(size(gpmtd%sparsepoints,2)))
    mtd_energy = 0.0d0

    allocate(mtd_var(size(gpmtd%sparsepoints,2)))
    mtd_var = 0.1d0

    do n=1,size(gpmtd%sparsepoints,2)

       ! go through available history
       current_history => hill_history
       do while( associated(current_history) )
          ! for every item in buffer
          do j = 1,current_history%nrst_of_values
             fexparg = 0.0d0
             do i = 1,NumOfMTDCVs
                diff = MTDCVList(i)%cv%get_deviation(gpmtd%sparsepoints(i,n),current_history%values(j,i))
                fexparg = fexparg + diff**2 / (2.0d0 * current_history%widths(j,i)**2)
             end do
             fh = current_history%heights(j)*exp(-fexparg)
          mtd_energy(n) = mtd_energy(n) + fh
          end do
          current_history => current_history%next_history_buffer
       end do

       ! add gp energy
       mtd_energy(n) = mtd_energy(n) + f_predict(gpmtd%gp,gpmtd%sparsepoints(:,n),SE_kernel_r_rr)
    end do

    ! do GP covariance matrix calculation
    call finalise(gpmtd%gp)

    call gp_basic_teach(self=gpmtd%gp, len_scale=gpmtd%thetas, periodicity=gpmtd%periodicities, &
                        f_var=gpmtd%fvar, kernel=SE_kernel_r_rr, &
                        f_r=gpmtd%sparsepoints, &
                        f_v=mtd_energy, &
                        f_n=mtd_var, &
                        jitter=fgpjitter)

    ! clear hill_history
    hill_history%nrst_of_values = 0

end subroutine mtd_core_gp_data

!===============================================================================
! Subroutine:  mtd_core_hill_scaling_factor
!===============================================================================

subroutine mtd_core_hill_scaling_factor(step,height_factor)

 use pmf_dat
 use mtd_dat

 implicit none
 integer, intent(in)                :: step
 real(PMFDP), intent(out)           :: height_factor
 integer                            :: dstep
 ! -----------------------------------------------------------------------------

 height_factor = 1.0d0

 dstep = fstep-step

 if ( dstep .ge. ftransition ) return

 select case(fscaling)
    ! none
    case(0)
        return
    ! step signal
    case(1)
        if( dstep .lt. fdelaying ) then
            height_factor = 0.0d0
        end if
    ! ramp signal
    case(2)
        height_factor = (dstep-fdelaying)/(ftransition-fdelaying)
 end select

return

end subroutine mtd_core_hill_scaling_factor

! ! ------------------------------------------------------------------------------
! 
! subroutine numerical_derivatives
! 
! use mtd_dat
! use pmf_dat
! use mtd_coordinates
! 
!  implicit none
!  integer                        :: i,j,k,jj
!  real(PMFDP)                    :: energy, loc_x(3,fcatom),mv1,mv2,fnumdiff
! 
!  fdx(:,:) = 0.0d0
! 
!  fnumdiff = 0.0001
! 
!  do i=1,NumOfMTDCVs
! 
!    do j=1,fcatom
!       do k=1,3
!          ! x =======================================================
! 
!          ! right point ----------------------------------------  
!          loc_x = lx
!          loc_x(k,j) = loc_x(k,j) + fnumdiff
! 
!          call calculate_meta_energy(loc_x,mv1)
! 
!          ! left point ----------------------------------------
!          loc_x = lx
!          loc_x(k,j) = loc_x(k,j) - fnumdiff
! 
!          call calculate_meta_energy(loc_x,mv2)
! 
!          fdx(k,j) = fdx(k,j) + (mv1 - mv2)/(2.0*fnumdiff)
! 
!       end do
!    end do 
! 
! end do
! 
! end subroutine numerical_derivatives
! 
! ! ------------------------------------------------------------------------------
! 
! subroutine calculate_meta_energy(x,energy)
! use mtd_dat
! use pmf_dat
! use mtd_coordinates
!  implicit none
!  real(PMFDP)                    :: x(:,:)
!  real(PMFDP)                    :: energy
! 
! ! Local variables -------------------------------------------------------------
! 
!  integer                            :: i,j,k,jj
!  type(history_buffer_rec),pointer   :: current_history
!  real(PMFDP)                        :: fexparg, mf, fh, wf
! 
!  ! ==========================================================
! 
!  do i = 1,NumOfMTDCVs
!     call calculate_coordinate(MTDCVList(i),x)
!  end do
! 
!  ! energy will be recalculated several times :-) maybe some further optimization of code?
!  energy = 0.0d0
!  ! go through history
!  current_history => hill_history
!  do while( associated(current_history) )
!     !for every item in buffer
!     do j = 1,current_history%nrst_of_values
!         fexparg = 0.0d0
!         do k = 1,NumOfMTDCVs
!             fexparg = fexparg + (MTDCVList(k)%value - current_history%values(j,k))**2 / &
!                                 (2.0d0 * current_history%widths(j,k)**2)
!         end do
!         fh = current_history%heights(j)*exp(-fexparg)
!         energy = energy + fh
!     end do
!     current_history => current_history%next_history_buffer
!  end do
! 
! end subroutine calculate_meta_energy

!===============================================================================

end module mtd_core
