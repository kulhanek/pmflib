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

module pmf_system_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine process_control

    use pmf_system_dat
    use prmfile
    use smf_utils
    use pmf_system
    use pmf_core
    use pmf_utils

    implicit none
    integer    ::  nrst_args
    ! -----------------------------------------------------------------------------

    ! get name of control file from the command line
    nrst_args = command_argument_count()

    ! check if control file was provided
    if( nrst_args .ne. 1 ) then
        call print_usage
        call pmf_utils_exit(PMF_OUT,1,'No control file specified on the command line or too many arguments!')
    end if

    ! process options
    call get_command_argument(nrst_args, ControlFile)   ! argument is name of control file

    ! write header ----------------------------------------------------------
    write(PMF_OUT,*)
    call centered_heading(PMF_OUT,'Reading control file', '-')
    write(PMF_OUT,'(a,a)') 'Control file name : ',trim(ControlFile)

    call prmfile_init(ControlPrmfile)

    if( .not. prmfile_read(ControlPrmfile,ControlFile) ) then
        call pmf_utils_exit(PMF_OUT,1,'Specified control file cannot be opened!')
    end if

    if( .not. prmfile_open_group(ControlPrmfile,'MAIN') ) then
        call pmf_utils_exit(PMF_OUT,1,'Specified control file does not contain {MAIN} group!')
    end if

    ! read user specificaton ------------------------------------------------
    call read_control
    call read_intervals
    call read_files
    call read_path
    call read_fes

    ! now we check if everything was understood from control file
    if( prmfile_count_ulines(ControlPrmfile,'MAIN') .gt. 0 ) then
        write(PMF_OUT,'(/,a,/)') '=== [unprocessed options in {MAIN} group] ======================================'
        call prmfile_dump_group(ControlPrmfile,PMF_OUT,'MAIN',.true.)
        call pmf_utils_exit(PMF_OUT,1,'Some items in the control file were not understood (see above)!')
    end if

    return

end subroutine process_control

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine print_usage

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a,/)') '=== [usage] ===================================================================='
    write(PMF_OUT,10)

 return

10 format('    pmf-path-explore <control>')

end subroutine print_usage

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_control

    use prmfile
    use pmf_system_dat
    use pmf_core
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [control] =================================================================='

    ! open first section
    if( .not. prmfile_open_section(ControlPrmfile,'control') ) then
        call pmf_utils_exit(PMF_OUT,1,'[control] section has to be specified in control file!')
    end if

    ! timings
    if(.not. prmfile_get_integer_by_key(ControlPrmfile,'steps', nsteps)) then
        call pmf_utils_exit(PMF_OUT,1,'Steps are not specified!')
    end if
    write(PMF_OUT,10) nsteps

    if( prmfile_get_real8_by_key(ControlPrmfile,'stepsize', stepsize)) then
        write(PMF_OUT,20) stepsize
    else
        write(PMF_OUT,25) stepsize
    end if

    if( prmfile_get_real8_by_key(ControlPrmfile,'smoothing', smoothing)) then
        write(PMF_OUT,40) smoothing
    else
        write(PMF_OUT,45) smoothing
    end if

    return

 10  format ('Number of steps (steps)                = ',i12)
 20  format ('Step size (stepsize)                   = ',f18.5)
 25  format ('Step size (stepsize)                   = ',f18.5,'            (default)')
 40  format ('Smoothing factor (sfac)                = ',f18.5)
 45  format ('Smoothing factor (sfac)                = ',f18.5,'            (default)')

end subroutine read_control

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_files_input_fes

    use pmf_system_dat
    use pmf_core
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [files] ===================================================================='

    if(.not. prmfile_open_section(ControlPrmfile,'files')) then
        call pmf_utils_exit(PMF_OUT,1,'[files] section was not found!')
    end if

    ! topology file
    if(.not. prmfile_get_string_by_key(ControlPrmfile,'fes', InputFESFile)) then
        call pmf_utils_exit(PMF_OUT,1,'Input free energy surface file is not specified!')
    end if
    write (PMF_OUT,10) trim(InputFESFile)

    return

 10 format('Input FES file (fes)                = ',a)

end subroutine read_files_input_fes

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_files

    use pmf_system_dat
    use pmf_core
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [files] ===================================================================='

    if(.not. prmfile_open_section(ControlPrmfile,'files')) then
        call pmf_utils_exit(PMF_OUT,1,'[files] section not found.')
    end if

    if( .not. prmfile_get_string_by_key(ControlPrmfile,'fes', InputFESFile)) then
        call pmf_utils_exit(PMF_OUT,1,'FES is not specified (fes)!')
    end if

    write (PMF_OUT,10) trim(InputFESFile)

    if( prmfile_get_string_by_key(ControlPrmfile,'input', InputPath)) then
        write (PMF_OUT,20) trim(InputPath)
    else
        write (PMF_OUT,21) trim(InputPath)
    end if

    if( prmfile_get_string_by_key(ControlPrmfile,'output', OutputPathFile)) then
        write (PMF_OUT,30) trim(OutputPathFile)
    else
        write (PMF_OUT,31) trim(OutputPathFile)
    end if

    if( prmfile_get_string_by_key(ControlPrmfile,'summary', OutputPathSummaryFile)) then
        write (PMF_OUT,40) trim(OutputPathSummaryFile)
    else
        write (PMF_OUT,41) trim(OutputPathSummaryFile)
    end if

    if( prmfile_get_string_by_key(ControlPrmfile,'trajectory', TrajectoryFile)) then
        write (PMF_OUT,50) trim(TrajectoryFile)
    else
        write (PMF_OUT,51) trim(TrajectoryFile)
    end if

    return

 10 format('Input FES file (fes)                   = ',a)
 20 format('Input path file (input)                = ',a)
 21 format('Input path file (input)                = ',a12,'                  (default)')
 30 format('Output path file (output)              = ',a)
 31 format('Output path file (output)              = ',a12,'                  (default)')
 40 format('Output path summary file (summary)     = ',a)
 41 format('Output path summary file (summary)     = ',a12,'                  (default)')
 50 format('Trajectory file (trajectory)           = ',a)
 51 format('Trajectory file (trajectory)           = ',a12,'                  (default)')

end subroutine read_files

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_intervals

    use prmfile
    use pmf_system_dat
    use pmf_core

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [intervals] ================================================================'

    if( .not. prmfile_open_section(ControlPrmfile,'intervals') ) then
        write(PMF_OUT,21) print_freq
        write(PMF_OUT,31) output_freq
        write(PMF_OUT,41) traj_freq
        write(PMF_OUT,51) smooth_freq
        write(PMF_OUT,61) reparam_freq
        return
    end if

    if(  prmfile_get_integer_by_key(ControlPrmfile,'summary', print_freq) ) then
        write(PMF_OUT,20) print_freq
    else
        write(PMF_OUT,21) print_freq
    end if

    if( prmfile_get_integer_by_key(ControlPrmfile,'output', output_freq) ) then
        write(PMF_OUT,30) output_freq
    else
        write(PMF_OUT,31) output_freq
    end if

    if( prmfile_get_integer_by_key(ControlPrmfile,'trajectory', traj_freq) ) then
        write(PMF_OUT,40) traj_freq
    else
        write(PMF_OUT,41) traj_freq
    end if

    if( prmfile_get_integer_by_key(ControlPrmfile,'smooth', smooth_freq) ) then
        write(PMF_OUT,50) smooth_freq
    else
        write(PMF_OUT,51) smooth_freq
    end if

    if( prmfile_get_integer_by_key(ControlPrmfile,'reparam', reparam_freq) ) then
        write(PMF_OUT,60) reparam_freq
    else
        write(PMF_OUT,61) reparam_freq
    end if

    ! avoid division by zero if interval is zero

    if( print_freq .le. 0 ) then
        print_freq = -1
    end if

    if( output_freq .le. 0 ) then
        output_freq = -1
    end if

    if( traj_freq .le. 0 ) then
        traj_freq = -1
    end if

    if( smooth_freq .le. 0 ) then
        smooth_freq = -1
    end if

    if( reparam_freq .le. 0 ) then
        reparam_freq = -1
    end if

 return

20 format('Progress interval (summary)            = ',i12)
21 format('Progress interval (summary)            = ',i12,'                  (default)')

30 format('Output interval (output)               = ',i12)
31 format('Output interval (output)               = ',i12,'                  (default)')

40 format('Trajectory write interval (trajectory) = ',i12)
41 format('Trajectory write interval (trajectory) = ',i12,'                  (default)')

50 format('Smoothing interval (smooth)            = ',i12)
51 format('Smoothing interval (smooth)            = ',i12,'                  (default)')

60 format('Reparametrization interval (reparam)   = ',i12)
61 format('Reparametrization interval (reparam)   = ',i12,'                  (default)')

end subroutine read_intervals

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_path

    use prmfile
    use pmf_system_dat
    use pmf_core
    use pmf_control
    use pmf_paths

    implicit none
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{PATHS}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( InputPath(1:1) .eq. '{' ) then
        grpname = InputPath(2:len_trim(InputPath)-1)
        write(PMF_OUT,110) grpname
        ! open goup with name from fpathdef
        if( .not. prmfile_open_group(ControlPrmfile,trim(grpname)) ) then
            write(PMF_OUT,130) trim(grpname)
            return
        end if
        write(PMF_OUT,*)
        call pmf_control_read_paths_gen_fake_cvs(ControlPrmfile)
        write(PMF_OUT,*)
        call pmf_control_read_paths_from_group(ControlPrmfile)
    else
        write(PMF_OUT,120) trim(InputPath)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,InputPath) ) then
            call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to load file: ' // trim(InputPath) // '!')
        end if

        write(PMF_OUT,*)
        call pmf_control_read_paths_gen_fake_cvs(locprmfile)

        write(PMF_OUT,*)
        call pmf_control_read_paths_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Path is read from group: ',A)
120 format('Path is read from file : ',A)
130 format(' >> No ',A,' group was specified.')

end subroutine read_path

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_fes

    use pmf_system_dat
    use pmf_core
    use pmf_utils

    implicit none
    integer            :: alloc_failed,i,it
    character(len=5)   :: sfes
    character(len=5)   :: sver
    ! ------------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{FES}',':')
    write(PMF_OUT,*)

    write(PMF_OUT,5) trim(InputFESFile)

    call pmf_utils_open(IO_FES,InputFESFile,'O')

    ! read header --------------------------
    read(IO_FES,10,end=100,err=100) sfes, sver, ncvs

    if( trim(sfes) .ne. 'FES' ) then
        call pmf_utils_exit(PMF_OUT,1,'Attempt to read non-FES file!')
    end if

    if( trim(sver) .ne. 'V1' ) then
        call pmf_utils_exit(PMF_OUT,1,'Attempt to read FES file of unsupported version!')
    end if

    if( ncvs .le. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'Two or more dimmensional FES is required!')
    end if

    allocate(cvs(ncvs),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for FES cvs!')
    end if

    tot_nbins = 1

    do i=1, ncvs
        read(IO_FES,20,end=100,err=100) it, cvs(i)%ctype, &
                                        cvs(i)%min_value, cvs(i)%max_value, cvs(i)%nbins
        cvs(i)%ctype = adjustl(cvs(i)%ctype)
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'Inccorect CV item in FES file!')
        end if
        read(IO_FES,25,end=100,err=100) it, cvs(i)%name
        cvs(i)%name = adjustl(cvs(i)%name)
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'Inccorect CV item in FES file!')
        end if
        tot_nbins = tot_nbins * cvs(i)%nbins
        cvs(i)%width = (cvs(i)%max_value - cvs(i)%min_value)/cvs(i)%nbins
    end do

    if( tot_nbins .le. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'FES is too small (only one bin or less)!')
    end if

    allocate(fes(tot_nbins), &
             dfes(ncvs,tot_nbins), &
             stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for FES data!')
    end if

    ! read energies
    read(IO_FES,40,end=100,err=100) (fes(i),i=1,tot_nbins)
    dfes(:,:) = 0.0d0

    close(IO_FES)

    call print_fes_info

 return

 5  format('FES is read from file : ',A)
10  format(A3,1X,A2,1X,I3)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
40  format(4(E19.11,1X))

100 call pmf_utils_exit(PMF_OUT,1,'Unable to read from FES file!')

 return

end subroutine read_fes

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine print_fes_info

    use pmf_system_dat

    implicit none
    integer            :: i
    ! -----------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,10) ncvs
    write(PMF_OUT,20) tot_nbins
    write(PMF_OUT,*)
    write(PMF_OUT,30)
    write(PMF_OUT,40)
    do i=1, ncvs
        write(PMF_OUT,50) i,trim(cvs(i)%ctype), &
                          trim(cvs(i)%name), &
                          cvs(i)%min_value,cvs(i)%max_value, &
                          cvs(i)%nbins
    end do

    return

 10 format('Number of CVs  : ',I9)
 20 format('Number of bins : ',I9)

 30 format('ID  Type      Name                          Min value       Max value     NBins ');
 40 format('-- ---------- -------------------------- --------------- --------------- -------');
 50 format(I2,1X,A10,1X,A26,1X,E15.8,1X,E15.8,1X,I7)

end subroutine print_fes_info

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmf_system_control
