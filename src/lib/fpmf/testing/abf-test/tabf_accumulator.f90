!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module tabf_accumulator

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  tabf_accumulator_init
!===============================================================================

subroutine tabf_accumulator_init()

    use tabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    allocate(accumulator%sizes(NumOfABFCVs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[TABF] Unable to allocate memory for abf accumulator!')
    endif

    tot_nbins = 1
    do i=1,NumOfABFCVs
        accumulator%sizes(i)%min_value  = ABFCVList(i)%min_value
        accumulator%sizes(i)%max_value  = ABFCVList(i)%max_value
        accumulator%sizes(i)%nbins      = ABFCVList(i)%nbins
        accumulator%sizes(i)%width      = abs(accumulator%sizes(i)%max_value - accumulator%sizes(i)%min_value)
        accumulator%sizes(i)%bin_width  = accumulator%sizes(i)%width / accumulator%sizes(i)%nbins
        tot_nbins = tot_nbins * accumulator%sizes(i)%nbins
    end do

    accumulator%tot_nbins = tot_nbins
    accumulator%tot_cvs = NumOfABFCVs

    ! ABF force arrays
    allocate(  &
            accumulator%nsamples(accumulator%tot_nbins), &
            accumulator%micf(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%m2icf(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%micf_pot(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%m2icf_pot(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%micf_kin(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%m2icf_kin(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%metot(accumulator%tot_nbins), &
            accumulator%m2etot(accumulator%tot_nbins), &
            accumulator%mepot(accumulator%tot_nbins), &
            accumulator%m2epot(accumulator%tot_nbins), &
            accumulator%mekin(accumulator%tot_nbins), &
            accumulator%m2ekin(accumulator%tot_nbins), &
            accumulator%merst(accumulator%tot_nbins), &
            accumulator%m2erst(accumulator%tot_nbins), &
            accumulator%cds_pp(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%cds_pk(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%cds_pr(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%cds_kp(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%cds_kk(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%cds_kr(NumOfABFCVs,accumulator%tot_nbins), &
            stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[TABF] Unable to allocate memory for abf accumulator (abfforce)!')
    endif

    call tabf_accumulator_clear()

    return

end subroutine tabf_accumulator_init

!===============================================================================
! Subroutine:  tabf_accumulator_clear
!===============================================================================

subroutine tabf_accumulator_clear()

    use tabf_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    accumulator%nsamples(:)     = 0

    accumulator%micf(:,:)       = 0.0d0
    accumulator%m2icf(:,:)      = 0.0d0
    accumulator%micf_pot(:,:)   = 0.0d0
    accumulator%m2icf_pot(:,:)  = 0.0d0
    accumulator%micf_kin(:,:)   = 0.0d0
    accumulator%m2icf_kin(:,:)  = 0.0d0
    accumulator%metot(:)        = 0.0d0
    accumulator%m2etot(:)       = 0.0d0
    accumulator%mepot(:)        = 0.0d0
    accumulator%m2epot(:)       = 0.0d0
    accumulator%mekin(:)        = 0.0d0
    accumulator%m2ekin(:)       = 0.0d0
    accumulator%merst(:)        = 0.0d0
    accumulator%m2erst(:)       = 0.0d0
    accumulator%cds_pp(:,:)     = 0.0d0
    accumulator%cds_pk(:,:)     = 0.0d0
    accumulator%cds_pr(:,:)     = 0.0d0
    accumulator%cds_kp(:,:)     = 0.0d0
    accumulator%cds_kk(:,:)     = 0.0d0
    accumulator%cds_kr(:,:)     = 0.0d0

end subroutine tabf_accumulator_clear

!===============================================================================
! Function:  accumulator_index
! Arguments:
!               idxcoord ... number of ksi coordinate
!               accuvalue    ... value that is used to compute the bin index
! compute index for one accumulator coordinate
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function tabf_accumulator_index(idxcoord,accuvalue)

    use tabf_dat
    use pmf_dat

    implicit none
    integer        :: idxcoord
    real(PMFDP)    :: accuvalue
    ! --------------------------------------------------------------------------

    ! we need number from zero - therefore we use floor(x)
    tabf_accumulator_index = floor((accuvalue - accumulator%sizes(idxcoord)%min_value) / &
                               accumulator%sizes(idxcoord)%bin_width)

    if( tabf_accumulator_index .lt. 0 .or. tabf_accumulator_index .ge.  accumulator%sizes(idxcoord)%nbins) then
        tabf_accumulator_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

    return

end function tabf_accumulator_index

!===============================================================================
! Function:  accumulator_globalindex
! Description:  Compute globalindex for accumulator, based on accuvalues of all coordinates
! Arguments:    none
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function tabf_accumulator_globalindex(lvalues)

    use tabf_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: lvalues(:)
    ! -----------------------------------------------
    integer                :: idx_local,i
    ! --------------------------------------------------------------------------

    tabf_accumulator_globalindex = 0

    do i=1,NumOfABFCVs
        idx_local = tabf_accumulator_index(i,lvalues(i))

        if (idx_local .eq. -1) then
            tabf_accumulator_globalindex = -1
            return
        end if

        tabf_accumulator_globalindex = tabf_accumulator_globalindex*accumulator%sizes(i)%nbins + idx_local
    end do

    tabf_accumulator_globalindex = tabf_accumulator_globalindex + 1

    return

end function tabf_accumulator_globalindex

!===============================================================================
! Subroutine:  tabf_accumulator_write
!===============================================================================

subroutine tabf_accumulator_write(iounit)

    use tabf_dat
    use pmf_dat
    use pmf_utils
    use pmf_unit

    implicit none
    integer                     :: iounit
    ! -----------------------------------------------
    integer                     :: i,j
    character(len=PMF_MAX_KEY)  :: key
    !---------------------------------------------------------------------------

    ! write header --------------------------
    write(iounit,10) 'ABF ', 'V6', NumOfABFCVs

    key = 'CVS'
    write(iounit,5) adjustl(key)
    do i=1, NumOfABFCVs
        write(iounit,20) i,trim(ABFCVList(i)%cv%ctype), &
                          accumulator%sizes(i)%min_value,accumulator%sizes(i)%max_value, &
                          accumulator%sizes(i)%nbins
        write(iounit,25) i,trim(ABFCVList(i)%cv%name)
        write(iounit,26) i,pmf_unit_get_rvalue(ABFCVList(i)%cv%unit,1.0d0),trim(pmf_unit_label(ABFCVList(i)%cv%unit))
    end do

    key = 'TEMPERATURE'
    write(iounit,5) adjustl(key)
    write(iounit,6) ftemp

    key = 'ENERGY'
    write(iounit,5) adjustl(key)
    write(iounit,27) pmf_unit_get_rvalue(EnergyUnit,1.0d0),trim(pmf_unit_label(EnergyUnit))

    key = 'NSAMPLES'
    write(iounit,5) adjustl(key)
    write(iounit,30) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)

    key = 'MICF'
    write(iounit,5) adjustl(key)
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%micf(i,j),j=1,accumulator%tot_nbins)
    end do
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%m2icf(i,j),j=1,accumulator%tot_nbins)
    end do

    key = 'MICF-POT'
    write(iounit,5) adjustl(key)
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%micf_pot(i,j),j=1,accumulator%tot_nbins)
    end do
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%m2icf_pot(i,j),j=1,accumulator%tot_nbins)
    end do

    key = 'MICF-KIN'
    write(iounit,5) adjustl(key)
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%micf_kin(i,j),j=1,accumulator%tot_nbins)
    end do
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%m2icf_kin(i,j),j=1,accumulator%tot_nbins)
    end do

    if( fenthalpy .or. fentropy ) then
        key = 'METOT'
        write(iounit,5) adjustl(key)
        write(iounit,40) (accumulator%metot(i),i=1,accumulator%tot_nbins)
        write(iounit,40) (accumulator%m2etot(i),i=1,accumulator%tot_nbins)
        key = 'MEPOT'
        write(iounit,5) adjustl(key)
        write(iounit,40) (accumulator%metot(i),i=1,accumulator%tot_nbins)
        write(iounit,40) (accumulator%m2etot(i),i=1,accumulator%tot_nbins)
        key = 'MEKIN'
        write(iounit,5) adjustl(key)
        write(iounit,40) (accumulator%metot(i),i=1,accumulator%tot_nbins)
        write(iounit,40) (accumulator%m2etot(i),i=1,accumulator%tot_nbins)
    end if

    if( fentropy ) then
        key = 'CDS-PP'
        write(iounit,5) adjustl(key)
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%cds_pp(i,j),j=1,accumulator%tot_nbins)
        end do
        key = 'CDS-PK'
        write(iounit,5) adjustl(key)
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%cds_pk(i,j),j=1,accumulator%tot_nbins)
        end do
        key = 'CDS-PR'
        write(iounit,5) adjustl(key)
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%cds_pr(i,j),j=1,accumulator%tot_nbins)
        end do
        key = 'CDS-KP'
        write(iounit,5) adjustl(key)
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%cds_kp(i,j),j=1,accumulator%tot_nbins)
        end do
        key = 'CDS-KK'
        write(iounit,5) adjustl(key)
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%cds_kk(i,j),j=1,accumulator%tot_nbins)
        end do
        key = 'CDS-KR'
        write(iounit,5) adjustl(key)
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%cds_kr(i,j),j=1,accumulator%tot_nbins)
        end do
    end if

    return

 5  format(A20)
 6  format(F10.4)
10  format(A4,1X,A2,1X,I2)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,A36)
27  format(3X,E18.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

end subroutine tabf_accumulator_write

!===============================================================================
! Subroutine:  tabf_accumulator_add_data_online
!===============================================================================

subroutine tabf_accumulator_add_data_online(cvs,gfx_pot,gfx_kin,epot,ekin,erst)

    use tabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx_pot(:)
    real(PMFDP)    :: gfx_kin(:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: ekin
    real(PMFDP)    :: erst
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, etot, icf
    real(PMFDP)    :: detot1, detot2
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: dekin1, dekin2
    real(PMFDP)    :: derst1, derst2
    real(PMFDP)    :: dicf1, dicf2
    real(PMFDP)    :: dicf_pot1, dicf_pot2
    real(PMFDP)    :: dicf_kin1, dicf_kin2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = tabf_accumulator_globalindex(cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    accumulator%nsamples(gi0) = accumulator%nsamples(gi0) + 1
    invn = 1.0d0 / real(accumulator%nsamples(gi0),PMFDP)

    etot = epot + ekin + erst

    ! total energy
    detot1 = etot - accumulator%metot(gi0)
    accumulator%metot(gi0)  = accumulator%metot(gi0)  + detot1 * invn
    detot2 = etot - accumulator%metot(gi0)
    accumulator%m2etot(gi0) = accumulator%m2etot(gi0) + detot1 * detot2

    ! potential energy
    depot1 = epot - accumulator%mepot(gi0)
    accumulator%mepot(gi0)  = accumulator%mepot(gi0)  + depot1 * invn
    depot2 = epot - accumulator%mepot(gi0)
    accumulator%m2epot(gi0) = accumulator%m2epot(gi0) + depot1 * depot2

    ! potential energy
    dekin1 = ekin - accumulator%mekin(gi0)
    accumulator%mekin(gi0)  = accumulator%mekin(gi0)  + dekin1 * invn
    dekin2 = ekin - accumulator%mekin(gi0)
    accumulator%m2ekin(gi0) = accumulator%m2ekin(gi0) + dekin1 * dekin2

    ! restraint energy
    derst1 = erst - accumulator%merst(gi0)
    accumulator%merst(gi0)  = accumulator%merst(gi0)  + derst1 * invn
    derst2 = erst - accumulator%merst(gi0)
    accumulator%m2erst(gi0) = accumulator%m2erst(gi0) + derst1 * derst2

    do i=1,NumOfABFCVs
        icf = gfx_pot(i) + gfx_kin(i)

        dicf1 = - icf - accumulator%micf(i,gi0)
        accumulator%micf(i,gi0)  = accumulator%micf(i,gi0)  + dicf1 * invn
        dicf2 = - icf -  accumulator%micf(i,gi0)
        accumulator%m2icf(i,gi0) = accumulator%m2icf(i,gi0) + dicf1 * dicf2

        dicf_pot1 = - gfx_pot(i) - accumulator%micf_pot(i,gi0)
        accumulator%micf_pot(i,gi0)  = accumulator%micf_pot(i,gi0)  + dicf_pot1 * invn
        dicf_pot2 = - gfx_pot(i) -  accumulator%micf_pot(i,gi0)
        accumulator%m2icf_pot(i,gi0) = accumulator%m2icf_pot(i,gi0) + dicf_pot1 * dicf_pot2

        dicf_kin1 = - gfx_kin(i) - accumulator%micf_kin(i,gi0)
        accumulator%micf_kin(i,gi0)  = accumulator%micf_kin(i,gi0)  + dicf_kin1 * invn
        dicf_kin2 = - gfx_kin(i) -  accumulator%micf_kin(i,gi0)
        accumulator%m2icf_kin(i,gi0) = accumulator%m2icf_kin(i,gi0) + dicf_kin1 * dicf_kin2

        accumulator%cds_pp(i,gi0)  = accumulator%cds_pp(i,gi0) + dicf_pot1 * depot2
        accumulator%cds_pk(i,gi0)  = accumulator%cds_pk(i,gi0) + dicf_pot1 * dekin2
        accumulator%cds_pr(i,gi0)  = accumulator%cds_pr(i,gi0) + dicf_pot1 * derst2
        accumulator%cds_kp(i,gi0)  = accumulator%cds_kp(i,gi0) + dicf_kin1 * depot2
        accumulator%cds_kk(i,gi0)  = accumulator%cds_kk(i,gi0) + dicf_kin1 * dekin2
        accumulator%cds_kr(i,gi0)  = accumulator%cds_kr(i,gi0) + dicf_kin1 * derst2
    end do

end subroutine tabf_accumulator_add_data_online

!===============================================================================
! Subroutine:  tabf_accumulator_get_data_lramp
!===============================================================================

subroutine tabf_accumulator_get_data_lramp(values,gfx)

    use tabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0, n
    real(PMFDP)    :: sc_ramp
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = tabf_accumulator_globalindex(values)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get number of samples
    n = accumulator%nsamples(gi0)
    if( n .gt. 0 ) then
        sc_ramp = 1.0d0
        if( n .le. fhramp_max ) then
            sc_ramp = 0.0d0
            if( n .gt. fhramp_min ) then
                sc_ramp = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
            end if
        end if
        gfx(:) = sc_ramp * accumulator%micf(:,gi0)
    end if

end subroutine tabf_accumulator_get_data_lramp

!===============================================================================

end module tabf_accumulator

