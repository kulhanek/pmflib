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

module pmfdyn_restraints_control

use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_restraints

    implicit none
    ! --------------------------------------------------------------------------

    call read_rstseq

end subroutine read_restraints

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_rstseq

    use prmfile
    use pmf_utils
    use pmfdyn_restraints_dat
    use pmfdyn_system_dat

    implicit none
    integer        :: i
    logical        :: yes
    character(80)  :: text
    character(10)  :: hflag,ceflag,chflag
    integer        :: alloc_failed
    ! -----------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [sequence_restraints] ======================================================'

    ! --- nrstr_seq, [rstseq]
    nrstr_seq = prmfile_open_section_and_count(ControlPrmfile,'sequence_restraints')
    write (PMF_OUT,10) nrstr_seq

    if ( nrstr_seq .le. 0 ) return

    ! allocate memory for rstseq
    allocate(rstseq(nrstr_seq), stat=alloc_failed)
        if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for rstseq array!')
    end if

    write (PMF_OUT,20)

    do i=1,nrstr_seq
        rstseq(i)%i = 0
        rstseq(i)%j = 0
        rstseq(i)%fk = 0
        rstseq(i)%ih = 0
        rstseq(i)%to_centre = 0
        rstseq(i)%change = 0
        rstseq(i)%fk_final = 0

        yes = prmfile_get_line(ControlPrmfile,text)
        read(text,*, end=111, err=112) rstseq(i)

    111 select case(rstseq(i)%ih)
            case(0)
                hflag = 'no'
            case(1)
                hflag = 'yes'
            case default
                call pmf_utils_exit(PMF_OUT,1,'Illegal H-flag for sequence restraint!')
        end select

        select case(rstseq(i)%to_centre)
            case(0)
                ceflag = 'atom-centres'
            case(1)
                ceflag = 'geo-centre'
            case(2)
                ceflag = 'mass-centre'
            case default
                call pmf_utils_exit(PMF_OUT,1,'Illegal to_centre flag for sequence restraint!')
        end select

        select case(rstseq(i)%change)
            case(0)
                chflag = 'no'
            case(1)
                chflag = 'yes'
            case default
                call pmf_utils_exit(PMF_OUT,1,'Illegal change flag for sequence restraint!')
        end select

        write(PMF_OUT,30) rstseq(i)%i,rstseq(i)%j,rstseq(i)%fk, &
                           trim(hflag),trim(ceflag),trim(chflag),rstseq(i)%fk_final
    end do

    return

112 call pmf_utils_exit(PMF_OUT,1,'Unable to read sequence restraint line!')

 10 format('Number of sequence restraints         = ',i12)
 20 format(/,'  atom_i  atom_j      fc  H-flag   to_centre  change    final_fc')
 30 format(2i8,f8.2,2X,A6,2X,A10,2X,A6,2X,f10.2)

end subroutine read_rstseq

!===============================================================================

end module pmfdyn_restraints_control
