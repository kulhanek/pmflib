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

module cst_restart

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_restart_read
!===============================================================================

subroutine cst_restart_read

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer :: i,inerr
    integer :: fpmode, fpnitem ! finger print info
    ! --------------------------------------------------------------------------

    ! read restart file - only master read file --------------------------------
    if( .not. frestart ) then
        write(CON_OUT,'(A)') '# RST: frestart = OFF'
        return
    end if

    ! check if restart file exist
    open(CON_RST,FILE=fcstrst,STATUS='OLD',ERR=1000)
    close(CON_RST)

    ! open restart file
    call pmf_utils_open(CON_RST,fcstrst,'O')

    ! read basic info - serves as finger print
    read(CON_RST,'(4I8)') fpmode, fpnitem

    ! check fingerprints
    inerr = 0
    if( fpmode .ne. fmode ) then
        write(PMF_OUT, '(/2x,a)') '[CST] fmode from restart file is not identical with this blue moon setup.'
        inerr = 1
    end if
    if( fpnitem .ne. NumOfCONs ) then
        write(PMF_OUT, '(/2x,a)') '[CST] NumOfCONs from restart file is not identical with this blue moon setup.'
        inerr = 1
    end if
    if (inerr .eq. 1) then
        write(PMF_OUT, '(/,a)') '[CST] restart error(s)!'
        call pmf_utils_exit(PMF_OUT, 1)
    end if

    ! read accumulation info
    read(CON_RST,'(2I8)') faccumulation,faccurst

    ! read current values
    read(CON_RST,'(E25.17)') (CONList(i)%value, i=1,NumOfCONs)

    ! read info about metric tenzor correction and lambda total
    read(CON_RST,'(E25.17)') isrztotal
    read(CON_RST,'(E25.17)') isrztotals
    read(CON_RST,'(E25.17)') (lambdatotal(i), i=1,NumOfCONs)
    read(CON_RST,'(E25.17)') (lambdatotals(i), i=1,NumOfCONs)
    read(CON_RST,'(E25.17)') (CONList(i)%sdevtot, i=1,NumOfCONs)

    close(CON_RST)

    write(CON_OUT,'(A)') '# RST: frestart = ON'

    return

1000 continue
    write(CON_OUT,*) '# WARNING: Unable to load restart information but frestart = on'

    frestart = .false.

    return

end subroutine cst_restart_read

!===============================================================================
! Subroutine:  cst_restart_write
!===============================================================================

subroutine cst_restart_write

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer :: i
    ! --------------------------------------------------------------------------

    ! only master process writes restart information
    call pmf_utils_open(CON_RST,fcstrst,'U')

    ! write basic info - serves as finger print
    write(CON_RST,'(4I8)') fmode, NumOfCONs

    ! write accumulation info
    write(CON_RST,'(2I8)') faccumulation,faccurst

    ! write current values
    write(CON_RST,'(E25.17)') (CONList(i)%value, i=1,NumOfCONs)

    ! write info about Z matrix, up part of BM, and lambda total
    write(CON_RST,'(E25.17)') isrztotal
    write(CON_RST,'(E25.17)') isrztotals
    write(CON_RST,'(E25.17)') (lambdatotal(i), i=1,NumOfCONs)
    write(CON_RST,'(E25.17)') (lambdatotals(i), i=1,NumOfCONs)
    write(CON_RST,'(E25.17)') (CONList(i)%sdevtot, i=1,NumOfCONs)

    close(CON_RST)

    return

end subroutine cst_restart_write

!===============================================================================

end module cst_restart

