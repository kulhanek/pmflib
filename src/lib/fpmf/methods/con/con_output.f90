!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module con_output

use pmf_sizes
use pmf_constants

contains

!===============================================================================
! Subroutine:  con_output_open
!===============================================================================

subroutine con_output_open

    use pmf_utils
    use pmf_dat
    use con_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! open output file
    call pmf_utils_open(CON_OUT,fconout,'R')

    write(CON_OUT,110)
    write(CON_OUT,120)
    write(CON_OUT,130)

    return

110 format('#===============================================================================')
120 format('# Constrained Dynamics')
130 format('#===============================================================================')

end subroutine con_output_open

!===============================================================================
! Subroutine:  con_output_write_header
!===============================================================================

subroutine con_output_write_header

    use pmf_constants
    use pmf_dat
    use con_dat

    implicit none
    integer             :: i,item, nitems
    type(UnitType)      :: lambda_unit
    character(len=15)   :: name
    ! --------------------------------------------------------------------------

    nitems = NumOfCONs
    if( fplevel .eq. 0 ) then
        nitems = NumOfCONs - NumOfSHAKECONs
    end if

    if( faccurst .gt. 0 ) then
        write(CON_OUT,'(A)')    '#-------------------------------------------------------------------------------'
        write(CON_OUT,'(A)')    '# INFO: THIS IS EQUILIBRATION STAGE OF SAMPLING'
        write(CON_OUT,'(A,I7)') '# Number of remaining steps:',faccurst
        write(CON_OUT,'(A)')    '#-------------------------------------------------------------------------------'
    else
        write(CON_OUT,'(A)')    '#-------------------------------------------------------------------------------'
        write(CON_OUT,'(A)')    '# INFO: THIS IS PRODUCTION PHASE OF ACCUMULATION'
        write(CON_OUT,'(A)')    '#-------------------------------------------------------------------------------'

        ! prevent double info about production stage
        if( faccurst .eq. 0 ) faccurst = -1
    end if


    write(CON_OUT,'(A)',advance='NO') '#   NSTEP    NACCU        TIME NI             MTC           <IRZ>        s(<IRZ>)'
    do i=1,nitems
        write(name,11) trim(CONList(i)%cv%name),i
        write(CON_OUT,10,advance='NO') trim(name),i
        write(CON_OUT,20,advance='NO') i,i,i
        if( has_lambdav ) write(CON_OUT,30,advance='NO') i,i,i
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#-------- -------- ----------- -- --------------- --------------- ---------------'
    do i=1,nitems
        write(CON_OUT,'(A)',advance='NO') ' --------------- ---------------'
        write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
        if( has_lambdav ) write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#                         [fs]   '
    write(CON_OUT,39,advance='NO') '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(CON_OUT,'(A)',advance='NO')     '          [i.u.]          [i.u.]'
    do i=1,nitems
        write(CON_OUT,40,advance='NO') '['//trim(CONList(i)%cv%get_ulabel())//']'

        lambda_unit = pmf_unit_div_units( EnergyUnit, CONList(i)%cv%unit )

        write(CON_OUT,45,advance='NO') '['//trim(pmf_unit_label(lambda_unit))//']'

        if( has_lambdav ) then
            write(CON_OUT,45,advance='NO') '['//trim(pmf_unit_label(lambda_unit))//']'
        end if
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#-------- -------- ----------- -- --------------- -------------------------------'
    do i=1,nitems
        write(CON_OUT,'(A)',advance='NO') ' -------------------------------'
        write(CON_OUT,'(A)',advance='NO') ' -----------------------------------------------'
        if( has_lambdav ) write(CON_OUT,'(A)',advance='NO') ' -----------------------------------------------'
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#       1        2           3  4               5               6               7'
    item = 8
    do i=1,nitems
        write(CON_OUT,5,advance='NO') item,item+1
        item = item + 2
        write(CON_OUT,15,advance='NO') item,item+1,item+2
        item = item + 3
        if( has_lambdav ) then
            write(CON_OUT,15,advance='NO') item,item+1,item+2
            item = item + 3
        end if
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#-------- -------- ----------- -- --------------- --------------- ---------------'
    do i=1,nitems
        write(CON_OUT,'(A)',advance='NO') ' --------------- ---------------'
        write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
        if( has_lambdav ) write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
    end do
    write(CON_OUT,*)

    return

 5 format('              ',I2,'              ',I2)
11 format(A,' = x',I2.2)
10 format(' ',A15,'        dev(x',I2.2,')')
15 format('              ',I2,'              ',I2,'              ',I2)
20 format('          L(x',I2.2,')        <L(x',I2.2,')>      s(<L(x',I2.2,'>)')
30 format('          M(x',I2.2,')        <M(x',I2.2,')>      s(<M(x',I2.2,'>)')
39 format(1X,A15)
40 format(1X,A31)
45 format(1X,A47)

end subroutine con_output_write_header

!===============================================================================
! Subroutine: con_output_write
!===============================================================================

subroutine con_output_write

    use pmf_dat
    use con_dat

    implicit none
    integer         :: i,ci,nitems
    real(PMFDP)     :: mtc,aisrz,isrzs,aisrzs
    real(PMFDP)     :: alam,lam,lams,alams
    real(PMFDP)     :: amu,mu,mus,amus
    type(UnitType)  :: cv_unit
    type(UnitType)  :: lambda_unit
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 .and. fstep .ne. fnstlim ) return ! write only every fsample step

    if( faccumulation .ne. 0 ) then
        if( faccumulation*isrztotals - isrztotal**2 .gt. 0 ) then
            isrzs = sqrt( faccumulation*isrztotals - isrztotal**2  ) / faccumulation
        else
            isrzs = 0.0d0
        end if
        aisrz = isrztotal / faccumulation
        mtc = PMF_Rgas*ftemp*log(aisrz)
        aisrzs = isrzs / sqrt(real(faccumulation))
    else
        mtc = 0.0d0     ! actually, this cannot be evaluated
        isrzs = 0.0d0
        aisrz = 0.0d0
        aisrzs = 0.0d0
    end if

    ! write header --------------------------------------------------------------
    write(CON_OUT,170,advance='NO') fstep,faccumulation,pmf_unit_get_rvalue(TimeUnit,ftime),fliter

    ! go to energy units - calculate metric tensor correction
    if( faccumulation .ne. 0 ) then
        write(CON_OUT,180,advance='NO') pmf_unit_get_rvalue(EnergyUnit,mtc), aisrz, aisrzs
    else
        write(CON_OUT,185,advance='NO')
    end if

    nitems = NumOfCONs
    if( fplevel .eq. 0 ) then
        nitems = NumOfCONs - NumOfSHAKECONs
    end if

    ! ---------------------------------------------------------------------------
    do i=1,nitems

        ! units
        cv_unit = CONList(i)%cv%unit
        lambda_unit = pmf_unit_div_units( EnergyUnit, cv_unit )

        ! lambda ----------------------------------
        lam = lambda(i)
        if( faccumulation .gt. 0 ) then
            if( faccumulation*lambdatotals(i) - lambdatotal(i)**2 .gt. 0 ) then
                lams = sqrt( faccumulation*lambdatotals(i) - lambdatotal(i)**2 ) / faccumulation   ! sigma(lambda)
            else
                lams = 0.0
            end if
            alam  = lambdatotal(i) / faccumulation
            alams = lams / sqrt(real(faccumulation))
        else
            alam = 0.0d0
            lams = 0.0d0
            alams = 0.0d0
        end if

        if( has_lambdav ) then
            mu = lambdav(i)
            if( faccumulation .gt. 0 ) then
                if ( faccumulation*lambdavtotals(i) - lambdavtotal(i)**2 .gt. 0 ) then
                    mus = sqrt( faccumulation*lambdavtotals(i) - lambdavtotal(i)**2 ) / faccumulation   ! sigma(lambda)
                else
                    mus = 0.0
                end if
                amu  = lambdavtotal(i) / faccumulation
                amus = mus / sqrt(real(faccumulation))
            else
                amu = 0.0d0
                mus = 0.0d0
                amus = 0.0d0
            end if
        end if

        ! write results ---------------------------------------------------------------
        ci = CONList(i)%cvindx
        write(CON_OUT,190,advance='NO') pmf_unit_get_rvalue(cv_unit,CVContextP%CVsValues(ci)), &
                                        pmf_unit_get_rvalue(cv_unit,CONList(i)%deviation)
        write(CON_OUT,200,advance='NO') pmf_unit_get_rvalue(lambda_unit,lam), &
                                        pmf_unit_get_rvalue(lambda_unit,alam), &
                                        pmf_unit_get_rvalue(lambda_unit,alams)
        if( has_lambdav ) then
            write(CON_OUT,200,advance='NO') pmf_unit_get_rvalue(lambda_unit,mu), &
                                            pmf_unit_get_rvalue(lambda_unit,amu), &
                                            pmf_unit_get_rvalue(lambda_unit,amus)
        end if
    end do
    write(CON_OUT,*)

    return

! output formats -----------------------------
170 format(1X,I8,1X,I8,1X,F11.1,1X,I2)
180 format(1X,E15.8,1X,E15.8,1X,E15.8)
185 format(1X,'---------------',1X,'---------------',1X,'---------------')
190 format(1X,E15.8,1X,E15.8)
200 format(1X,E15.8,1X,E15.8,1X,E15.8)

end subroutine con_output_write

!===============================================================================
! Subroutine:  con_output_close
!===============================================================================

subroutine con_output_close

    use pmf_dat
    use con_dat

    implicit none
    integer             :: i,item, nitems
    real(PMFDP)         :: mtc,aisrz,isrzs,aisrzs
    real(PMFDP)         :: alam,lams,alams
    real(PMFDP)         :: amu,mus,amus
    type(UnitType)      :: cv_unit
    type(UnitType)      :: lambda_unit
    character(len=15)   :: name
    ! -----------------------------------------------------------------------------

    nitems = NumOfCONs
    if( fplevel .eq. 0 ) then
        nitems = NumOfCONs - NumOfSHAKECONs
    end if

    write(CON_OUT,'(A)') '#'
    write(CON_OUT,'(A)',advance='NO') '#   NSTEP    NACCU        TIME             MTC'

    do i=1,nitems
        write(name,11) trim(CONList(i)%cv%name),i
        write(CON_OUT,39,advance='NO') trim(name)
        write(CON_OUT,20,advance='NO') i,i,i
        if( has_lambdav ) write(CON_OUT,30,advance='NO') i,i,i
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#-------- -------- ----------- ---------------'
    do i=1,nitems
        write(CON_OUT,'(A)',advance='NO') ' ---------------'
        write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
        if( has_lambdav ) write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#                         [fs]'
    write(CON_OUT,39,advance='NO') '['//trim(pmf_unit_label(EnergyUnit))//']'
    do i=1,nitems
        write(CON_OUT,39,advance='NO') '['//trim(CONList(i)%cv%get_ulabel())//']'

        lambda_unit = pmf_unit_div_units( EnergyUnit, CONList(i)%cv%unit )

        write(CON_OUT,45,advance='NO') '['//trim(pmf_unit_label(lambda_unit))//']'

        if( has_lambdav ) then
            write(CON_OUT,45,advance='NO') '['//trim(pmf_unit_label(lambda_unit))//']'
        end if
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#-------- -------- ----------- ---------------'
    do i=1,nitems
        write(CON_OUT,'(A)',advance='NO') ' ---------------'
        write(CON_OUT,'(A)',advance='NO') ' -----------------------------------------------'
        if( has_lambdav ) write(CON_OUT,'(A)',advance='NO') ' -----------------------------------------------'
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#       2        3           4               5               6'
    item = 7
    do i=1,nitems
        write(CON_OUT,5,advance='NO') item
        item = item + 1
        write(CON_OUT,15,advance='NO') item,item+1,item+2
        item = item + 3
        if( has_lambdav ) then
            write(CON_OUT,15,advance='NO') item,item+1,item+2
            item = item + 3
        end if
    end do
    write(CON_OUT,*)

    write(CON_OUT,'(A)',advance='NO') '#-------- -------- ----------- ---------------'
    do i=1,nitems
        write(CON_OUT,'(A)',advance='NO') ' ---------------'
        write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
        if( has_lambdav ) write(CON_OUT,'(A)',advance='NO') ' --------------- --------------- ---------------'
    end do
    write(CON_OUT,*)

    if ( faccumulation*isrztotals - isrztotal**2 .gt. 0 ) then
        isrzs = sqrt( faccumulation*isrztotals - isrztotal**2 ) / faccumulation    ! sigma(down)
    else
        isrzs = 0.0
    end if

    if( faccumulation .ge. 1 ) then
        aisrz = isrztotal / faccumulation
        aisrzs = isrzs / sqrt(faccumulation*1.0)
    else
        aisrz = 0.0
        aisrzs = 0.0
    end if

    ! calculate metric tensor correction
    mtc = PMF_Rgas * ftemp * log(aisrz)

    ! write header --------------------------------------------------------------

    write(CON_OUT,170,advance='NO') fnstlim,faccumulation,pmf_unit_get_rvalue(TimeUnit,ftime)
    write(CON_OUT,180,advance='NO') pmf_unit_get_rvalue(EnergyUnit,mtc)

    ! ---------------------------------------------------------------------------
    do i=1,nitems

        ! units
        cv_unit = CONList(i)%cv%unit
        lambda_unit = pmf_unit_div_units( EnergyUnit, cv_unit )

        ! lambda ----------------------------------
        if( faccumulation .gt. 0 ) then
            if ( faccumulation*lambdatotals(i) - lambdatotal(i)**2 .gt. 0 ) then
                lams = sqrt( faccumulation*lambdatotals(i) - lambdatotal(i)**2 ) / faccumulation   ! sigma(lambda)
            else
                lams = 0.0
            end if
            alam     = lambdatotal(i) / faccumulation
            alams    = lams / sqrt(real(faccumulation))
        else
            alam = 0.0d0
            lams = 0.0d0
            alams = 0.0d0
        end if

        if( has_lambdav ) then
            if( faccumulation .gt. 0 ) then
                if ( faccumulation*lambdavtotals(i) - lambdavtotal(i)**2 .gt. 0 ) then
                    mus = sqrt( faccumulation*lambdavtotals(i) - lambdavtotal(i)**2 ) / faccumulation   ! sigma(lambda)
                else
                    mus = 0.0
                end if
                amu     = lambdavtotal(i) / faccumulation
                amus    = mus / sqrt(real(faccumulation))
            else
                amu = 0.0d0
                mus = 0.0d0
                amus = 0.0d0
            end if
        end if

        ! write results ---------------------------------------------------------------
        write(CON_OUT,190,advance='NO') pmf_unit_get_rvalue(cv_unit,CONList(i)%value)
        write(CON_OUT,200,advance='NO') pmf_unit_get_rvalue(lambda_unit,alam), &
                                        pmf_unit_get_rvalue(lambda_unit,lams), &
                                        pmf_unit_get_rvalue(lambda_unit,alams)
        if( has_lambdav ) then
            write(CON_OUT,200,advance='NO') pmf_unit_get_rvalue(lambda_unit,amu), &
                                        pmf_unit_get_rvalue(lambda_unit,mus), &
                                        pmf_unit_get_rvalue(lambda_unit,amus)
        end if
    end do

    write(CON_OUT,*)

    ! close file -----------------------------------

    close(CON_OUT)

    return

 5 format('              ',I2)
11 format(A,' = x',I2.2)
15 format('              ',I2,'              ',I2,'              ',I2)
20 format('        <L(x',I2.2,')>       s(L(x',I2.2,'))      s(<L(x',I2.2,'>)')
30 format('        <M(x',I2.2,')>       s(M(x',I2.2,'))      s(<M(x',I2.2,'>)')

170 format('#',I8,1X,I8,1X,F11.1)
180 format(1X,E15.8)
190 format(1X,E15.8)
200 format(1X,E15.8,1X,E15.8,1X,E15.8)

39 format(1X,A15)
40 format(1X,A31)
45 format(1X,A47)

end subroutine con_output_close

!===============================================================================

end module con_output

