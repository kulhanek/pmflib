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

module pmfdyn_thermostat

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine init_thermostat_subsystem(x,v)

    use pmf_utils
    use pmfdyn_thermostat_dat
    use pmfdyn_system_dat
    use cst_dat

    implicit none
    real(PMFDP)     :: x(:,:)
    real(PMFDP)     :: v(:,:)
    !----------------------------------------------------
    real(PMFDP)     :: ncom, nbm, tmpE
    !---------------------------------------------------------------------------

    write(PMF_OUT,10)

    select case (ThermostatType)
        case(THERMOSTAT_NONE)
            write(PMF_OUT,50) 'none'
        case(THERMOSTAT_BERENDSEN)
            write(PMF_OUT,50) 'berendsen'
            write(PMF_OUT,60) Tau_T
            write(PMF_OUT,80) Iseed
        case(THERMOSTAT_CHEATHAM)
            write(PMF_OUT,50) 'cheatham'
            write(PMF_OUT,60) Tau_T
            write(PMF_OUT,80) Iseed
        case(THERMOSTAT_LANGEVIN)
            write(PMF_OUT,50) 'langevin'
            write(PMF_OUT,70) Gamma_T
            write(PMF_OUT,80) Iseed
        case default
            call pmf_utils_exit(PMF_OUT,1,'>> ERROR: Thermostst is not implemented yet!')
    end select

    ! langevin stuff --------------------------------------------------------------
    select case (ThermostatType)
        case(THERMOSTAT_LANGEVIN)
            lgam_c = 1.d0 + 0.5d0 * Gamma_T * Stepsize
            lgam_p = 1.d0 / lgam_c
            lgam_m = 1.d0 - 0.5d0 * Gamma_T * Stepsize
    end select

    if( Temp1 < 0 ) then
        write(PMF_OUT,20) Temp0
    else
        write(PMF_OUT,25)
        write(PMF_OUT,30) Temp0
        write(PMF_OUT,40) Temp1
    end if
    write(PMF_OUT,72) prmfile_onoff(MonitorHotAtoms)
    if( MonitorHotAtoms ) then
        write(PMF_OUT,76) HotAtomTreshold
        write(PMF_OUT,74) prmfile_onoff(RemoveHotAtoms)
    end if

    select case (ThermostatType)
        case(THERMOSTAT_NONE,THERMOSTAT_CHEATHAM,THERMOSTAT_BERENDSEN)
            ncom = 6.0d0
        case(THERMOSTAT_LANGEVIN)
            ncom = 0.0d0
    end select

    nbm = real(NumOfCONs)

    NOF = 3.0d0*natoms - ncom - nbm

    ! calculate number of freedom
    write(PMF_OUT,110) 3.0d0*real(natoms)
    write(PMF_OUT,120) ncom
    write(PMF_OUT,140) nbm
    write(PMF_OUT,150) NOF

    ! set up random generator
    call RLUXGO(3,Iseed,0,0)

    if( (restart .eqv. .false.) .or. (Tmaxw .ge. 0) ) then
        call generate_MB_velocities(x,v)
    end if

    call adjust_thermostat_controls(v)

    call get_temperature(v,v,tmpE)

    write(PMF_OUT,90) TempA

    return

 10 format(/,'Initializing thermostat ...')
 50 format  ('   Thermostat                              = ',a10)
 60 format  ('   Bath coupling                           = ',f12.1,' [fs]')
 70 format  ('   Friction coefficient                    = ',f15.4,' [fs^-1]')
 80 format  ('   Random seed                             = ',i10)
 72 format  ('   Monitor hot atoms                       = ',a10)
 74 format  ('   Remove hot atoms                        = ',a10)
 76 format  ('   Hot atom treshold                       = ',f12.1,' [K]')

110 format(/,'   Degree of freedom (3*N_atoms)           = ',f12.1)
120 format  ('   Number of COM constraints               = ',f12.1)
140 format  ('   Number of PMFLib constraints            = ',f12.1)
150 format  ('   Degree of freedom                       = ',f12.1)

 20 format  ('   Temperature                             = ',f12.1,' [K]')
 25 format  ('   Temperature will be changed in following interval:')
 30 format  ('   Starting temperature                    = ',f12.1,' [K]')
 40 format  ('   Target temperature                      = ',f12.1,' [K]')

 90 format(/,'   Current temperature T                   = ',f12.1,' [K]')

end subroutine init_thermostat_subsystem

! ------------------------------------------------------------------------------

subroutine generate_MB_velocities(x,v)

    use pmfdyn_thermostat_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)     :: x(:,:)
    real(PMFDP)     :: v(:,:)
    !----------------------------------------------------
    integer         :: i,j
    real(PMFDP)     :: sd,scale_fac
    real(PMFDP)     :: tmpE
    ! --------------------------------------------------------------------------

    write(PMF_OUT,10) Tmaxw

    ! Generate Maxwell-Boltzman velocities ---------------
    do i=1,natoms
    sd = sqrt(PMF_Rgas * Tmaxw * winv(i))
        do j=1,3
            v(j,i)= rand_gaussian(sd)
        end do
    end do

    ! get current temp -----------------------------------
    call get_temperature(v,v,tmpE)
    write(PMF_OUT,20) TempA

    ! remove COM translational and rotational moments -----
    write(PMF_OUT,40)
    call stop_com(x,v)

    ! get current temp ------------------------------------
    call get_temperature(v,v,tmpE)
    write(PMF_OUT,20) TempA

    ! fix temperature -------------------------------------
    scale_fac = sqrt ( Tmaxw/TempA )
    write(PMF_OUT,60)
    v(:,:) = v(:,:) * scale_fac

    return

10 format(/,'   Generating new velocities -> T          = ',f12.1)
 20 format ('        Real temperature T                 = ',f12.1)
 40 format ('     Removing translational and rotational moment ...')
 50 format ('     Initial SHAKE ...')
 60 format ('     Rescaling to correct temperature ...')

end subroutine generate_MB_velocities

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine adjust_thermostat_controls(v)

    use pmfdyn_thermostat_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)                    :: v(:,:)
    !----------------------------------------------------
    integer                        :: i
    real(PMFDP)                    :: ekin0,scale_fac
    ! -----------------------------------------------------------------------------

    ! calculate target temperature
    TempT = Temp0
    if( Temp1 .ge. 0 ) then
        TempT = Temp0 + (Temp1-Temp0)*real(Istep)/real(Nsteps)
    end if

    ! recalculate treshold for hot atoms
    ! 1.5 = (3/2), eg. three degree of freedoms per atom
    HotAtomMaxEkin = HotAtomTreshold*1.5d0*PMF_Rgas*TempT

    ! find and remove hot atoms --------------------------------------------------
    if( MonitorHotAtoms ) then
        do i=1,natoms
            ekin0 = 0.5*mass(i)*(v(1,i)**2 + &
                                 v(2,i)**2 + &
                                 v(3,i)**2)
            ! is it a hot atom?
            if( ekin0 .gt. HotAtomMaxEkin ) then
                ! report hot atom
                write (PMF_OUT,10) istep, i, ekin0/(1.5d0*PMF_Rgas)
                if( RemoveHotAtoms ) then
                    ! remove this hot atom - rescale its velocities to target temperature
                    scale_fac = sqrt(ekin0/(1.5d0*PMF_Rgas*TempT))
                    v(:,i) = v(:,i) * scale_fac
                end if
            end if
        end do
    end if

    return

10 format ('>>> WARNING: Step = ',i10,' Hot atom i =',i10,' Temp(i)=',f10.2)

end subroutine adjust_thermostat_controls

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine get_temperature(oldv,v,E)

    use pmfdyn_thermostat_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)     :: oldv(:,:)
    real(PMFDP)     :: v(:,:)
    real(PMFDP)     :: E
    !----------------------------------------------------
    integer         :: i
    real(PMFDP)     :: ekin0,temp_fac
    ! -----------------------------------------------------------------------------

    ! calculate kinetic energies
    E = 0.0

    temp_fac = 1.0d0
    if( ThermostatType .eq. THERMOSTAT_LANGEVIN ) temp_fac = lgam_c

    do i=1,natoms
        ekin0 = temp_fac*0.5d0*0.25d0*mass(i)*((oldv(1,i)+v(1,i))**2 + &
                                           (oldv(2,i)+v(2,i))**2 + &
                                           (oldv(3,i)+v(3,i))**2)
        E = E + ekin0
    end do

    ! calculate temperatures
    TempA = 2.0*E / (PMF_Rgas*real(NOF))

    return

end subroutine get_temperature

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine berendsen_thermostat(v)

    use pmfdyn_thermostat_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)     :: v(:,:)
    !----------------------------------------------------
    real(PMFDP)     :: scale_fac
    ! --------------------------------------------------------------------------

    ! calculate scaling factor
    scale_fac = 0.0d0
    if ( TempA .ne. 0 ) scale_fac = TempT/TempA - 1.0
    scale_fac = sqrt ( 1 + stepsize * scale_fac / Tau_T  )

    ! scale velocities
    v(:,:) = v(:,:) * scale_fac
    return

end subroutine berendsen_thermostat

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine cheatham_thermostat(oldv,v)

    use pmfdyn_thermostat_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)     :: oldv(:,:)
    real(PMFDP)     :: v(:,:)
    !----------------------------------------------------
    integer         :: i
    real(PMFDP)     :: scale_fac
    real(PMFDP)     :: tekin, ekin0, ekin1, ekin2
    ! --------------------------------------------------------------------------

    ! kinetic energy representation of target temperature
    tekin = 0.5 * NOF * PMF_Rgas * TempT

    ! recalculate kinetic energies - we need more, optimization is not important now
    ekin0 = 0.0
    ekin1 = 0.0
    ekin2 = 0.0

    do i=1,natoms
        ekin0 = ekin0 + 0.5*mass(i)*(oldv(1,i)**2 + &
                                    oldv(2,i)**2 + &
                                    oldv(3,i)**2)
        ekin1 = ekin1 + 0.5*0.25*mass(i)*((oldv(1,i)+v(1,i))**2 + &
                                        (oldv(2,i)+v(2,i))**2 + &
                                        (oldv(3,i)+v(3,i))**2)
        ekin2 = ekin2 + 0.5*mass(i)*(v(1,i)**2 + &
                                    v(2,i)**2 + &
                                    v(3,i)**2)
    end do

    ! calculate scaling factor
    scale_fac = 0.0d0
    if ( (ekin0 + ekin2) .ne. 0 ) scale_fac = (tekin - ekin1)/(ekin0 + ekin2)
    scale_fac = sqrt ( 1 + 2.0d0 * scale_fac * stepsize / tau_T )

    ! scale velocities
    v(:,:) = v(:,:) * scale_fac
    return

end subroutine cheatham_thermostat

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine langevin_thermostat(v,d)

    use pmfdyn_thermostat_dat
    use pmfdyn_system_dat

    implicit none
    real(PMFDP)     :: v(:,:)
    real(PMFDP)     :: d(:,:)
    !----------------------------------------------------
    integer         :: i
    real(PMFDP)     :: randF(3)
    real(PMFDP)     :: rfac
    real(PMFDP)     :: rsi
    ! -----------------------------------------------------------------------------

    ! do langevin integration step
    rfac = sqrt(2.0d0 * Gamma_T * PMF_Rgas * TempT / Stepsize) / PMF_DT2VDT

    do i=1,natoms
        rsi = rfac * sqrt(mass(i))
        randF(1) = rand_gaussian(rsi)
        randF(2) = rand_gaussian(rsi)
        randF(3) = rand_gaussian(rsi)
        v(:,i) = ( v(:,i)*lgam_m + (randF(:) - d(:,i))*winv(i)*dt ) * lgam_p
    end do

    return

end subroutine langevin_thermostat

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine stop_com(x,v)

    use pmfdyn_system_dat

    implicit none
    real(PMFDP)                    :: x(:,:)
    real(PMFDP)                    :: v(:,:)
    !----------------------------------------------------
    integer                        :: j,m,n
    real(PMFDP)                    :: x1, x2, x3
    real(PMFDP)                    :: xx, xy, xz, yy, yz, zz
    real(PMFDP)                    :: xcom(3)
    real(PMFDP)                    :: vcom(3)
    real(PMFDP)                    :: acom(3)
    real(PMFDP)                    :: ocom(3)
    real(PMFDP)                    :: tcom(3,3)
    real(PMFDP)                    :: mmass, totmass
    real(PMFDP)                    :: crit
    integer                        :: indx(3),work(9)
    integer                        :: info
    ! -----------------------------------------------------------------------------

    crit = 1.0d-06

    ! calculate COM, translational velocity of COM, and angular momentum

    xcom(:) = 0.d0
    vcom(:) = 0.d0
    acom(:) = 0.d0
    ocom(:) = 0.d0
    totmass = 0.0

    do j = 1,natoms
        mmass = mass(j)
        xcom(:)  = xcom(:) + x(:,j)*mmass
        vcom(:)  = vcom(:) + v(:,j)*mmass
        acom(1) = acom(1) + (x(2,j)*v(3,j)-x(3,j)*v(2,j)) * mmass
        acom(2) = acom(2) + (x(3,j)*v(1,j)-x(1,j)*v(3,j)) * mmass
        acom(3) = acom(3) + (x(1,j)*v(2,j)-x(2,j)*v(1,j)) * mmass
        totmass = totmass + mmass
    end do

    xcom(:)  = xcom(:) / totmass
    vcom(:)  = vcom(:) / totmass

    acom(1) = acom(1) - (xcom(2)*vcom(3) - xcom(3)*vcom(2)) * totmass
    acom(2) = acom(2) - (xcom(3)*vcom(1) - xcom(1)*vcom(3)) * totmass
    acom(3) = acom(3) - (xcom(1)*vcom(2) - xcom(2)*vcom(1)) * totmass

    xx = 0.d0
    xy = 0.d0
    xz = 0.d0
    yy = 0.d0
    yz = 0.d0
    zz = 0.d0

    do j = 1,natoms
        x1 = x(1,j) - xcom(1)
        x2 = x(2,j) - xcom(2)
        x3 = x(3,j) - xcom(3)
        mmass = mass(j)
        xx = xx + x1*x1*mmass
        xy = xy + x1*x2*mmass
        xz = xz + x1*x3*mmass
        yy = yy + x2*x2*mmass
        yz = yz + x2*x3*mmass
        zz = zz + x3*x3*mmass
    end do

    tcom(1,1) = yy+zz
    tcom(2,1) = -xy
    tcom(3,1) = -xz
    tcom(1,2) = -xy
    tcom(2,2) = xx+zz
    tcom(3,2) = -yz
    tcom(1,3) = -xz
    tcom(2,3) = -yz
    tcom(3,3) = xx+yy

    ! invert tensor
    call dgetrf(3,3,tcom,3,indx,info)
    call dgetri(3,tcom,3,indx,work,9,info)

    do m = 1,3
        ocom(m) = 0.d0
        do n = 1,3
            ocom(m) = ocom(m) + tcom(m,n)*acom(n)
        end do
    end do

    ! stop COM translational motion
    do j = 1,natoms
        v(:,j) = v(:,j) - vcom(:)
    end do

    ! stop COM rotational motion
    do j = 1,natoms
        x1 = x(1,j) - xcom(1)
        x2 = x(2,j) - xcom(2)
        x3 = x(3,j) - xcom(3)
        v(1,j) = v(1,j) - ocom(2)*x3 + ocom(3)*x2
        v(2,j) = v(2,j) - ocom(3)*x1 + ocom(1)*x3
        v(3,j) = v(3,j) - ocom(1)*x2 + ocom(2)*x1
    end do

    if( istep .gt. 0 ) write(PMF_OUT,10) istep
    return

10 format('Translational and rotational motion removed in step : ',i8)

end subroutine stop_com

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

! Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122
! rewritten from GSL library 1.3

real(PMFDP) function rand_gaussian(sigma)

    implicit none
    real(PMFDP)       :: sigma
    real(PMFDP)       :: x,y,r2
    ! --------------------------------------------------------------------------

    r2 = 0.0d0

    do while( (r2 .gt. 1.0d0) .or. (r2 .eq. 0d0) )
        ! choose x,y in uniform square (-1,-1) to (+1,+1)
        x = -1 + 2 * rand_uniform()
        y = -1 + 2 * rand_uniform()

        ! see if it is in the unit circle
        r2 = x * x + y * y
    end do

    ! Box-Muller transform
    rand_gaussian = sigma * y * sqrt (-2.0 * log (r2) / r2)

end function rand_gaussian

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

real(PMFDP) function rand_uniform()

    implicit none
    real(PMFSP)   vector(1)
    ! --------------------------------------------------------------------------

    call RANLUX(vector,1)
    rand_uniform = vector(1)

end function rand_uniform

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmfdyn_thermostat

