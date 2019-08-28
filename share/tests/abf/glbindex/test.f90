program test_glbindex

    implicit none
    integer,parameter   :: PMFDP            = 8
    integer             :: NumOfABFCVs
    real(PMFDP)         :: min_value(10)
    real(PMFDP)         :: max_value(10)
    integer             :: nbins(10)
    real(PMFDP)         :: width(10)
    real(PMFDP)         :: bin_width(10)
    real(PMFDP)         :: cv_values(10)
    integer             :: tot_nbins
    integer             :: glbindex, old_glbindex, i, j, k
    real(PMFDP)         :: cv1, cv2, cv3
    ! --------------------------------------------------------------------------

    NumOfABFCVs = 1
    min_value(1) = -8.0d0
    max_value(1) =  8.0d0
    nbins(1) = 160
!    call abf_accumulator_init

!    cv1 = -9.0
!    do i=1,200
!        cv_values(1) = cv1
!        glbindex = abf_accumulator_globalindex(cv_values)
!        write(*,*) 'cv1 = ', cv_values(1), ' glbidx = ',  glbindex
!        cv1 = cv1 + 0.01
!    end do

    NumOfABFCVs = 2
    min_value(2) =  0.0d0
    max_value(2) =  3.3d0
    nbins(2) = 33
    call abf_accumulator_init

    cv1 = -7.95
    old_glbindex = 0
    do i=1,160
        cv_values(1) = cv1
        cv2 = 0.05
        do j=1,33
            cv_values(2) = cv2
            glbindex = abf_accumulator_globalindex(cv_values)
            write(*,*) 'cv1 = ', cv_values(1), 'cv2 = ', cv_values(2), ' glbidx = ',  glbindex
            if( old_glbindex + 1 .ne. glbindex ) then
                stop 'glbindex discontinuity'
            end if
            cv2 = cv2 + 0.1
            old_glbindex = glbindex
        end do
        cv1 = cv1 + 0.1
    end do
    write(*,*) 'tot_nbins = ', tot_nbins

    NumOfABFCVs = 3
    min_value(3) =  0.0d0
    max_value(3) =  3.0d0
    nbins(3) = 30
    call abf_accumulator_init

    cv1 = -7.95
    old_glbindex = 0
    do i=1,160
        cv_values(1) = cv1
        cv2 = 0.05
        do j=1,33
            cv_values(2) = cv2
            cv3 = 0.05
            do k=1,30
                cv_values(3) = cv3
                glbindex = abf_accumulator_globalindex(cv_values)
                write(*,*) 'cv1 = ', cv_values(1), 'cv2 = ', cv_values(2), 'cv3 = ', cv_values(3), ' glbidx = ',  glbindex
                if( old_glbindex + 1 .ne. glbindex ) then
                    stop 'glbindex discontinuity'
                end if
                old_glbindex = glbindex
                cv3 = cv3 + 0.1
            end do
            cv2 = cv2 + 0.1
        end do
        cv1 = cv1 + 0.1
    end do
    write(*,*) 'tot_nbins = ', tot_nbins

contains

!===============================================================================
! Subroutine:  abf_accumulator_init
!===============================================================================

subroutine abf_accumulator_init()

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    tot_nbins = 1
    do i=1,NumOfABFCVs
        width(i)      = abs(max_value(i) - min_value(i))
        bin_width(i)  = width(i) / nbins(i)
        tot_nbins = tot_nbins * nbins(i)
    end do

end subroutine abf_accumulator_init

!===============================================================================
! Function:  accumulator_index
! Arguments:
!               idxcoord ... number of ksi coordinate
!               accuvalue    ... value that is used to compute the bin index
! compute index for one accumulator coordinate
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function abf_accumulator_index(idxcoord,accuvalue)

    implicit none
    integer        :: idxcoord
    real(PMFDP)    :: accuvalue
    ! --------------------------------------------------------------------------

    ! we need number from zero - therefore we use floor(x)
    abf_accumulator_index = floor((accuvalue - min_value(idxcoord)) / &
                               bin_width(idxcoord))

    if( abf_accumulator_index .lt. 0 .or. abf_accumulator_index .ge.  nbins(idxcoord)) then
        abf_accumulator_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

    return

end function abf_accumulator_index

!===============================================================================
! Function:  accumulator_globalindex
! Description:  Compute globalindex for accumulator, based on accuvalues of all coordinates
! Arguments:    none
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function abf_accumulator_globalindex(lvalues)

    implicit none
    real(PMFDP)            :: lvalues(:)
    ! -----------------------------------------------
    integer                :: idx_local,i
    ! --------------------------------------------------------------------------

    abf_accumulator_globalindex = 0

    do i=1,NumOfABFCVs
        idx_local = abf_accumulator_index(i,lvalues(i))

        if (idx_local .eq. -1) then
            abf_accumulator_globalindex = -1
            return
        end if

        abf_accumulator_globalindex = abf_accumulator_globalindex*nbins(i) + idx_local
    end do

    abf_accumulator_globalindex = abf_accumulator_globalindex + 1

    return

end function abf_accumulator_globalindex

!===============================================================================

end program test_glbindex
