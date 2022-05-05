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

module pmf_init

implicit none
contains

!===============================================================================
! Subroutine:  pmf_init_dat
!===============================================================================

subroutine pmf_init_dat

    use pmf_dat
    use pmf_timers
    use pmf_unit

    implicit none
    ! --------------------------------------------------------------------------

    ftopology      = ''          ! fake topology
    fdebug         = .false.     ! more verbose output
    frepmpifrag    = .false.
    fprint_inpcrds = .false.
    fprint_masks   = .false.
    fenable_pbc    = .false.
    fmonitor_paths = .false.

    fnatoms     = 0
    fnstlim     = 0
    fdt         = 0.0d0
    ftime       = 0.0d0
    ftemp       = 0.0d0
    fstep       = 0
    fsystype    = SYS_UNK
    fintalg     = IA_LEAP_FROG
    fshake      = .false.

    fdtx        = 0.0d0
    ifdtx       = 0.0d0

    fcanexmdloop = .false.
    fexit_mdloop = 0

    pmf_enabled  = .false.
    cst_enabled  = .false.
    rst_enabled  = .false.
    mtd_enabled  = .false.
    abf_enabled  = .false.
    abp_enabled  = .false.
    mon_enabled  = .false.
    stm_enabled  = .false.
    pdrv_enabled = .false.

    tabf_enabled    = .false.
    usabf_enabled   = .false.

    lng_force_required      = .false.
    shake_force_required    = .false.


    LNG_c_implic    = 1.0d0

    fucell(:,:) = 0.0d0
    frecip(:,:) = 0.0d0
    fbox_volume = 0.0d0
    fbox_sphere = 0.0d0

    ! file names ---------------------------------------------------------------
    fcvsdef     = '{CVS}'
    fpathsdef   = '{PATHS}'

    fcstdef     = '{CST}'
    fcstout     = '_cst.out'
    fcstrst     = '_cst.rst'
    fcstfrst    = '_cst.rst-full'
    fcsttrj     = '_cst.trj'

    frstdef     = '{RST}'
    frstout     = '_rst.out'
    frsthist    = '_rst.hist'

    fmtddef     = '{MTD}'
    fmtdout     = '_mtd.out'
    fmtdrst     = '_mtd.rst'
    fmtdtrj     = '_mtd.trj'
    fmtdhills   = '_mtd.hills'

    fabfdef     = '{ABF}'
    fabfmask    = '_abf.mask'
    fabfout     = '_abf.out'
    fabfrst     = '_abf.rst'
    fabftrj     = '_abf.trj'

    ftabfdef    = '{TABF}'
    ftabfout    = '_tabf.out'
    ftabfrst    = '_tabf.rst'
    ftabftrj    = '_tabf.trj'
    ftabficf    = '_tabf.icf'

    fusabfdef   = '{US-ABF}'
    fusabfout   = '_us-abf.out'
    fusabfrst   = '_us-abf.rst'
    fusabftrj   = '_us-abf.trj'

    fabpdef     = '{ABP}'
    fabpout     = '_abp.out'
    fabprst     = '_abp.rst'
    fabptrj     = '_abp.trj'

    fstmdef     = '{STM}'
    fstmout     = '_stm.out'

    fmondef     = '{MON}'
    fmonout     = '_mon.out'

    fpdrvdef    = '{PDRV}'
    fpdrvout    = '_pdrv.out'

    ! units --------------------------------------
    call pmf_unit_decode_mass_unit('amu',MassUnit)
    call pmf_unit_decode_time_unit('fs',TimeUnit)
    call pmf_unit_decode_energy_unit('kcal/mol',EnergyUnit)
    call pmf_unit_decode_length_unit('A',LengthUnit)
    call pmf_unit_decode_angle_unit('rad',AngleUnit)
    call pmf_unit_decode_temp_unit('K',TemperatureUnit)

    ! timers -------------------------------------
    call pmf_timers_init

end subroutine pmf_init_dat

#ifdef MPI
!===============================================================================
! Subroutine:  pmf_init_taskid_mpi
!===============================================================================

subroutine pmf_init_taskid_mpi(mytaskid,numoftasks)

    use pmf_dat

    implicit none
    integer        :: mytaskid
    integer        :: numoftasks
    ! --------------------------------------------------------------------------

    fmytaskid = mytaskid
    fnumoftasks = numoftasks
    fmaster = fmytaskid .eq. 0


end subroutine pmf_init_taskid_mpi

!===============================================================================
! Subroutine:  pmf_init_bcast_dat_mpi
!===============================================================================

subroutine pmf_init_bcast_dat_mpi

    use pmf_utils
    use pmf_dat
    use mpi

    implicit none
    integer        :: alloc_failed,ierr
    ! --------------------------------------------------------------------------

    ierr = MPI_SUCCESS

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A)') '>> Broadcasting dat section (only master is reporting)'
    end if

    ! integers --------------------------------------
    call mpi_bcast(fdebug, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast fdebug variable!')
    end if

    call mpi_bcast(NumOfLAtoms, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast NumOfLAtoms variable!')
    end if

    call mpi_bcast(fnatoms, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast fnatoms variable!')
    end if

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A)') '>> Broadcasting dat section - end'
        write(PMF_DEBUG+fmytaskid,'(A,L)') '   fdebug      = ',fdebug
        write(PMF_DEBUG+fmytaskid,'(A,I3)') '   NumOfLAtoms = ',NumOfLAtoms
        write(PMF_DEBUG+fmytaskid,'(A,I3)') '   fnatoms      = ',fnatoms
        write(PMF_DEBUG+fmytaskid,*)
    end if

    if( NumOfLAtoms .eq. 0 ) return ! no atoms are there

    ! allocate arrays on slaves ---------------------
    if( .not. fmaster ) then
        allocate(RIndexes(NumOfLAtoms),    &
                stat= alloc_failed )
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to allocate RIndexes!')
        end if
    end if

    ! transfer RIndexes
    call mpi_bcast(RIndexes, NumOfLAtoms, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast RIndexes array!')
    end if

    ! transfer cst_enabled status
    call mpi_bcast(cst_enabled, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast cst_enabled variable!')
    end if

end subroutine pmf_init_bcast_dat_mpi
#endif

!===============================================================================
! Subroutine:  pmf_init_variables
!===============================================================================

subroutine pmf_init_variables(intalg,natom,systype,nstlim,dt,time,temp)

    use pmf_dat
    use pmf_core
    use pmf_utils

    implicit none
    integer        :: intalg
    integer        :: natom    ! total number of atoms
    integer        :: systype  ! system type
    integer        :: nstlim   ! length in simulation in steps
    real(PMFDP)    :: dt       ! dt of step
    real(PMFDP)    :: time     ! actual time
    real(PMFDP)    :: temp     ! simulation temperature
    ! --------------------------------------------------------------------------

    ! MD setup --------------------------------------
    fintalg     = intalg
    fnatoms     = natom
    fnstlim     = nstlim
    fdt         = dt*TimeConv
    ftime       = time*TimeConv
    ftemp       = temp
    fsystype    = systype

    fdtx  = fdt*PMF_DT2VDT

    if( fdtx .le. 0.0d0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] fdt must be greater than zero!')
    end if

    ifdtx = 1.0d0 /  fdtx

end subroutine pmf_init_variables

!===============================================================================
! Subroutine:  pmf_init
!===============================================================================

subroutine pmf_init_all_nocvvalues(amass,ax)

    use pmf_dat
    use pmf_core
    use pmf_cvs
    use pmf_mask
    use pmf_pbc

    implicit none
    real(PMFDP)     :: amass(:)
    real(PMFDP)     :: ax(:,:)
    ! -----------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    ! init all pmf arrays
    call pmf_init_pmf()

    call pmf_init_sys_summary()

    ! if no atoms - exit
    if( NumOfLAtoms .eq. 0 ) return

    ! copy initial crds and masses
    do i=1,NumOfLAtoms
        Mass(i)     = amass(RIndexes(i))*MassConv
        if( Mass(i) .ne. 0 ) then
            MassInv(i) = 1.0d0 / Mass(i)
        else
            MassInv(i) = 0.0
        end if
        Crd(:,i) = ax(:,RIndexes(i))*LengthConv
    end do

end subroutine pmf_init_all_nocvvalues

!===============================================================================
! Subroutine:  pmf_init
!===============================================================================

subroutine pmf_init_all(amass,ax)

    use pmf_dat
    use pmf_core
    use pmf_cvs
    use pmf_mask
    use pmf_pbc
    use pmf_paths

    implicit none
    real(PMFDP)     :: amass(:)
    real(PMFDP)     :: ax(:,:)
    ! -----------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    ! init subsystems
    call pmf_init_all_nocvvalues(amass,ax)

    ! init CVs
    if( NumOfCVs .gt. 0 ) then
        CVContext%CVsValues = 0.0d0
        CVContext%CVsDrvs = 0.0d0
    end if

    do i=1,NumOfCVs
        call CVList(i)%cv%calculate_cv(Crd,CVContext)
    end do

    ! init PATHs
    do i=1,NumOfPATHs
        call pmf_paths_get_path_current_alpha(PathList(i)%path,CVContext)
    end do

    ! try to write imaged coordinates
    if( fdebug ) then
        write(PMF_OUT,'(A)') '# PMF DEBUG: writing imaged coordinates'
        call pmf_pbc_write_xyz(Crd,'imaged.xyz',.true.)
    end if

end subroutine pmf_init_all

!===============================================================================
! Subroutine:  pmf_init_pmf
!===============================================================================

subroutine pmf_init_pmf

    use pmf_dat
    use pmf_utils
    use cst_init

    implicit none
    integer                :: i,j,k,l,itmp,ai,tot_atoms
    logical                :: sorted
    integer                :: alloc_failed
    integer,allocatable    :: tmp_indexes(:)
    ! --------------------------------------------------------------------------

    NumOfLAtoms = 0

    ! update CV and CST according to SHAKE constraints
    call cst_init_add_shake_csts

    !---------------------------------------------------------------------------
    ! calculate total number of atoms including duplicity records
    tot_atoms = 0
    do i=1,NumOfCVs
        tot_atoms = tot_atoms + CVList(i)%cv%natoms
    end do

    if( tot_atoms .eq. 0 ) return

    ! allocate index array
    allocate(tmp_indexes(tot_atoms),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMFLIB] Unable to allocate memory for tmp_indexes array!')
    end if

    ! fill array by indexes ----------------
    ai = 1
    do i=1,NumOfCVs
        do j=1,CVList(i)%cv%natoms
            tmp_indexes(ai) = CVList(i)%cv%rindexes(j)
            ai = ai + 1
        end do
    end do

    ! sort atom indexes by brutal force ----
    sorted = .false.
    do while(.not. sorted)
        sorted = .true.
        do i=1,tot_atoms-1
            if( tmp_indexes(i) .gt. tmp_indexes(i+1) ) then
                itmp = tmp_indexes(i+1)
                tmp_indexes(i+1) = tmp_indexes(i)
                tmp_indexes(i) = itmp
                sorted = .false.
            end if
        end do
    end do

    ! calcuate number of unique atoms
    NumOfLAtoms = 1
    do i=2,tot_atoms
        if( tmp_indexes(i-1) .ne. tmp_indexes(i) ) then
            NumOfLAtoms = NumOfLAtoms + 1
        end if
    end do

    ! init real indexes -----------------------------
    allocate(RIndexes(NumOfLAtoms),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,&
                            '[PMFLIB] Unable to allocate memory for RIndexes array!')
    endif

    do i=1,NumOfCVs
        CVList(i)%cv%idx = i

        do j=1,CVList(i)%cv%natoms
            ! find atom in tmp_indexes
            ai = 1
            do k=2,tot_atoms
                if( CVList(i)%cv%rindexes(j) .eq. tmp_indexes(k-1) ) exit
                if( tmp_indexes(k-1) .ne. tmp_indexes(k) ) then
                    ai = ai + 1
                end if
            end do

            CVList(i)%cv%lindexes(j) = ai
            RIndexes(ai) = CVList(i)%cv%rindexes(j)
        end do
    end do

    ! find individual atoms for each CV lam81
    do i=1,NumOfCVs
        CVList(i)%cv%nindatoms = CVList(i)%cv%natoms
        do j=1,CVList(i)%cv%natoms
           do k=1, j-1
              if( CVList(i)%cv%lindexes(j) .eq. CVList(i)%cv%lindexes(k) ) then
                  CVList(i)%cv%nindatoms = CVList(i)%cv%nindatoms - 1
                  exit
              end if
           end do
        end do
        allocate(CVList(i)%cv%indlindexes(CVList(i)%cv%nindatoms), stat=alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,&
                                '[PMFLIB] Unable to allocate memory for indlindexes array!')
        endif
        l = 0
 outer: do j=1,CVList(i)%cv%natoms
            do k=1, j-1
               if( CVList(i)%cv%lindexes(j) .eq. CVList(i)%cv%lindexes(k) ) then
                   cycle outer
               end if
            end do
            l = l + 1
            CVList(i)%cv%indlindexes(l) = CVList(i)%cv%lindexes(j)
        end do outer
    end do

    ! release temporary array
    deallocate(tmp_indexes)

    ! allocate remaining arrays
    allocate(InitialCrd(3,NumOfLAtoms), &
          Mass(NumOfLAtoms), &
          MassInv(NumOfLAtoms), &
          Crd(3,NumOfLAtoms), &
          Frc(3,NumOfLAtoms), &
          SHAKEFrc(3,NumOfLAtoms), &
          Vel(3,NumOfLAtoms), &
          CVContext%CVsValues(NumOfCVs), &
          CVContext%CVsDrvs(3,NumOfLAtoms,NumOfCVs), &
          stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMFLIB] Unable to allocate memory for common arrays!')
    endif

    if( cst_enabled .or. shake_force_required ) then
        ! allocate arrays used by bluemoon
        allocate( CrdP(3,NumOfLAtoms), &
                  CVContextP%CVsValues(NumOfCVs), &
                  CVContextP%CVsDrvs(3,NumOfLAtoms,NumOfCVs), &
                  CrdBar(3,NumOfLAtoms), &
                  stat=alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[PMFLIB] Unable to allocate memory for CST arrays!')
        endif

        if( fintalg .eq. IA_VEL_VERLET ) then
            ! allocate arrays used by bluemoon
            allocate(VelP(3,NumOfLAtoms), &
                      stat=alloc_failed)

            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT, 1,'[PMFLIB] Unable to allocate memory for CST array!')
            endif
        end if
    end if

end subroutine pmf_init_pmf

!===============================================================================
! Subroutine:  pmf_init_pmf_methods
!===============================================================================

subroutine pmf_init_pmf_methods()

    use pmf_dat
    use mon_init
    use rst_init
    use abf_init
    use abp_init
    use mtd_init
    use cst_init
    use stm_init
    use pdrv_init
    use tabf_init
    use usabf_init

    implicit none
    ! --------------------------------------------------------------------------

    if( pdrv_enabled ) then
        call pdrv_init_method
    end if

    if( rst_enabled ) then
        call rst_init_method
    end if

    if( abf_enabled ) then
        call abf_init_method
    end if

    if( tabf_enabled ) then
        call tabf_init_method
    end if

    if( usabf_enabled ) then
        call usabf_init_method
    end if

    if( abp_enabled ) then
        call abp_init_method
    end if

    if( mtd_enabled ) then
        call mtd_init_method
    end if

    if( stm_enabled ) then
        call stm_init_method
    end if

    if( cst_enabled ) then
        call cst_init_method
    end if

    if( mon_enabled ) then
        call mon_init_method
    end if

    pmf_enabled = abf_enabled .or. abp_enabled .or. mtd_enabled .or. stm_enabled &
               .or. cst_enabled .or. rst_enabled .or. mon_enabled .or. pdrv_enabled &
               .or. tabf_enabled .or. usabf_enabled

end subroutine pmf_init_pmf_methods

!===============================================================================
! Subroutine: pmf_init_title
!===============================================================================

subroutine pmf_init_title(driver_name)

    use pmf_constants
    use pmf_dat
    use pmf_ver

    implicit none
    character(*)       :: driver_name
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,'(A)')   '#==============================================================================#'
    write(PMF_OUT,'(A)')   '# PMFLib - Potential of Mean Force Toolkit                                     #'
    write(PMF_OUT,'(A)')   '# -----------------------------------------------------------------------------#'
    write(PMF_OUT,'(A)')   '# Authors: (c) 2019 - 2021 Petr Kulhanek (NCBR)                                #'
    write(PMF_OUT,'(A)')   '#          (c) 2011 - 2015 Petr Kulhanek (CEITEC)                              #'
    write(PMF_OUT,'(A)')   '#          (c) 2013 - 2015 Letif Mones (EDUC)                                  #'
    write(PMF_OUT,'(A)')   '#          (c) 2009 - 2011 Petr Kulhanek (NCBR)                                #'
    write(PMF_OUT,'(A)')   '#          (c) 2008 Petr Kulhanek (IE), Letif Mones (IE)                       #'
    write(PMF_OUT,'(A)')   '#          (c) 2007 Petr Kulhanek (IE), Martin Petrek (NCBR), Letif Mones (IE) #'
    write(PMF_OUT,'(A)')   '#          (c) 2006 Petr Kulhanek (NCBR), Martin Petrek (NCBR)                 #'
    write(PMF_OUT,'(A)')   '#          (c) 2005 Petr Kulhanek (NCBR)                                       #'
    write(PMF_OUT,'(A)')   '#                                                                              #'
    write(PMF_OUT,'(A)')   '# CEITEC: Central European Institute of Technology, Masaryk University, CZ     #'
    write(PMF_OUT,'(A)')   '# EDUC:   Engineering Department, University of Cambridge, UK                  #'
    write(PMF_OUT,'(A)')   '# IE  :   Institute of Enzymology, Hungarian Academy of Science, HU            #'
    write(PMF_OUT,'(A)')   '# NCBR:   National Centre for Biomolecular Research, Masaryk University, CZ    #'
    write(PMF_OUT,'(A)')   '#                                                                              #'
    write(PMF_OUT,'(A)')   '# Fortran and C++ Cores of PMFLib are license under Lesser GPL v2.1 and above  #'
    write(PMF_OUT,'(A)')   '# PMFLib utilities are license under GPL v2 and above                          #'
    write(PMF_OUT,'(A)')   '#==============================================================================#'
    write(PMF_OUT,'(A,A)') '# PMFLib version : ',trim(PMFLIBVER)
    write(PMF_OUT,'(A,A)') '# Driver         : ',trim(driver_name)

    DriverName = driver_name

end subroutine pmf_init_title

!===============================================================================
! Subroutine: pmf_init_title
!===============================================================================

subroutine pmf_init_sys_summary()

    use pmf_constants
    use pmf_dat
    use pmf_utils
    use pmf_unit
    use prmfile

    implicit none
    character(20)       :: form_str
    type(UnitType)      :: volunit
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,'(A)') '#==============================================================================#'
    write(PMF_OUT,'(A)') '# PMFLib - System summary                                                      #'
    write(PMF_OUT,'(A)') '#==============================================================================#'
    write(PMF_OUT,'(A)') '#'

    write(form_str,'(I7)') fnatoms
    write(PMF_OUT,10) adjustl(form_str)

    write(form_str,'(I7)') NumOfCVs
    write(PMF_OUT,20) adjustl(form_str)

    write(form_str,'(I7)') NumOfLAtoms
    write(PMF_OUT,30) adjustl(form_str)

    write(PMF_OUT,35) prmfile_onoff(fshake)

    write(form_str,'(F10.1)') pmf_unit_get_rvalue(TemperatureUnit,ftemp)
    write(PMF_OUT,40)  trim(adjustl(form_str)), '['//trim(pmf_unit_label(TemperatureUnit))//']'

    write(form_str,'(F10.3)') pmf_unit_get_rvalue(TimeUnit,fdt)
    write(PMF_OUT,50)  trim(adjustl(form_str)), '['//trim(pmf_unit_label(TimeUnit))//']'

    write(form_str,'(I10)') fnstlim
    write(PMF_OUT,60) trim(adjustl(form_str))

    write(form_str,'(F20.2)') pmf_unit_get_rvalue(TimeUnit,fdt*fnstlim)
    write(PMF_OUT,65)  trim(adjustl(form_str)), '['//trim(pmf_unit_label(TimeUnit))//']'

    select case(fsystype)
    case(SYS_NT)
    write(PMF_OUT,70)  'NT'
    case(SYS_NTV)
    write(PMF_OUT,70)  'NTV'
    case(SYS_NTP)
    write(PMF_OUT,70)  'NTp'
    end select

    select case(fbox_type)
    case(BOX_ISOLATED_SYSTEM)
    write(PMF_OUT,80)  'isolated system'
    case(BOX_ORTHOGONAL)
    write(PMF_OUT,80)  'orthogonal'
    case(BOX_GENERAL)
    write(PMF_OUT,80)  'general'
    end select

    volunit = pmf_unit_power_unit(LengthUnit,3)

    select case(fbox_type)
    case(BOX_ORTHOGONAL,BOX_GENERAL)
    write(form_str,'(F15.1)') pmf_unit_get_rvalue(volunit,fbox_volume)
    write(PMF_OUT,90)  trim(adjustl(form_str)),'['//trim(pmf_unit_label(volunit))//']'

    write(PMF_OUT,100) pmf_unit_get_rvalue(LengthUnit,fucell(1,1)), &
                    pmf_unit_get_rvalue(LengthUnit,fucell(2,1)), &
                    pmf_unit_get_rvalue(LengthUnit,fucell(3,1)), &
                    '['//trim(pmf_unit_label(LengthUnit))//']'
    write(PMF_OUT,110) pmf_unit_get_rvalue(LengthUnit,fucell(1,2)), &
                    pmf_unit_get_rvalue(LengthUnit,fucell(2,2)), &
                    pmf_unit_get_rvalue(LengthUnit,fucell(3,2)), &
                    '['//trim(pmf_unit_label(LengthUnit))//']'
    write(PMF_OUT,120) pmf_unit_get_rvalue(LengthUnit,fucell(1,3)), &
                    pmf_unit_get_rvalue(LengthUnit,fucell(2,3)), &
                    pmf_unit_get_rvalue(LengthUnit,fucell(3,3)), &
                    '['//trim(pmf_unit_label(LengthUnit))//']'

    write(form_str,'(F10.2)') pmf_unit_get_rvalue(LengthUnit,fbox_sphere)
    write(PMF_OUT,130) trim(adjustl(form_str)), '['//trim(pmf_unit_label(LengthUnit))//']'
    end select

 10 format('# Total number of atoms    : ',A)
 20 format('# Total number of CVs      : ',A)
 30 format('# Number of atoms involved : ',A)
 35 format('# SHAKE                    : ',A)
 40 format('# Temperature              : ',A,1X,A)
 50 format('# Time step                : ',A,1X,A)
 60 format('# Number of steps          : ',A)
 65 format('# Simulation length        : ',A,1X,A)
 70 format('# System type              : ',A)
 80 format('# Box type                 : ',A)
 90 format('# Initial box volume       : ',A,1X,A)
100 format('# Initial box A vector     : ',F12.4,1X,F12.4,1X,F12.4,1X,A)
110 format('# Initial box B vector     : ',F12.4,1X,F12.4,1X,F12.4,1X,A)
120 format('# Initial box C vector     : ',F12.4,1X,F12.4,1X,F12.4,1X,A)
130 format('# Largest sphere radius    : ',A,1X,A)

end subroutine pmf_init_sys_summary

!===============================================================================

end module pmf_init
