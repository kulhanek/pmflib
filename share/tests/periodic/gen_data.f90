program gen_data

    implicit none
    integer           :: i
    double precision  :: xvalue, yvalue, noise
    double precision  :: pvalue, pivalue
    double precision  :: ivalue, avalue, dvalue
    double precision  :: mavalue, mdvalue
    ! ---------------------------------------
    double precision  :: minv
    double precision  :: maxv
    double precision  :: interval
    ! ---------------------------------------------------------------------------

    minv =  1.0
    maxv =  7.0

    interval = (maxv-minv)

    pvalue = minv + 0.5*interval
    pivalue = minv + 0.5*interval

    do i=1,1000
        xvalue = 2.0*3.14*real(i)/360.0
        yvalue = sin(xvalue)
        call random_number(noise)
        noise  = 0.2*(noise-0.5)
        yvalue = interval*(yvalue+noise+1.0)/2.0 + minv

        ! image
        ivalue = yvalue
        if( ivalue .gt. maxv ) ivalue = ivalue - interval
        if( ivalue .lt. minv ) ivalue = ivalue + interval

        ! average --------------
        ! correct values
        avalue = 0.5*(yvalue + pvalue)
        if( avalue .gt. maxv ) avalue = avalue - interval
        if( avalue .lt. minv ) avalue = avalue + interval

        ! tested values
        mavalue = get_average_value(ivalue,pivalue)

        if( abs(mavalue-avalue) .gt. 1e-7 ) then
            write(*,10) xvalue, yvalue, ivalue, avalue, mavalue
            stop 1
        end if

        ! difference ----------
        ! correct values
        dvalue = (yvalue - pvalue)
        ! tested values
        mdvalue = get_deviation(ivalue,pivalue)

        if( abs(mdvalue-dvalue) .gt. 1e-7 ) then
            write(*,20) xvalue, yvalue, ivalue, dvalue, mdvalue
            stop 1
        end if

        write(*,*) xvalue, yvalue, ivalue

        pvalue = yvalue
        pivalue = ivalue
    end do

 5 format(3(1x,F10.6))
10 format('Average error',5(1x,F10.6))
20 format('Differe error',5(1x,F10.6))

contains

!===============================================================================

double precision function get_average_value(value1,value2)

    implicit none
    double precision     :: value1
    double precision     :: value2
    ! --------------------------------------------
    double precision     :: vec
    ! --------------------------------------------------------------------------

    if( abs(value2-value1) .lt. 0.5d0*(maxv-minv) ) then
        get_average_value = 0.5d0*(value1+value2)
        return
    else
        ! get vector
        vec = value2 - value1
        ! shift to box center
        vec = vec + 0.5d0*(maxv+minv)
        ! image as point
        vec = vec - (maxv-minv)*floor((vec - minv)/(maxv-minv))
        ! return vector back
        vec = vec - 0.5d0*(maxv+minv)
        ! calculate average
        get_average_value = 0.5d0*vec + value1
        ! image average
        get_average_value = get_average_value - (maxv-minv)*floor((get_average_value - minv)/(maxv-minv))

        ! debug
        !write(*,*) 'get_average_value:',value1,value2,get_average_value
    end if

end function get_average_value

!===============================================================================

double precision function get_deviation(value1,value2)

    implicit none
    double precision     :: value1
    double precision     :: value2
    ! --------------------------------------------
    double precision     :: vec
    ! --------------------------------------------------------------------------

    if( abs(value1-value2) .lt. 0.5d0*(maxv-minv) ) then
        get_deviation = value1 - value2
        return
    else
        ! get vector
        vec = value1 - value2
        ! shift to box center
        vec = vec + 0.5d0*(maxv+minv)
        ! image as point
        vec = vec - (maxv-minv)*floor((vec-minv)/(maxv-minv))
        ! return vector back
        get_deviation = vec - 0.5d0*(maxv+minv)

        ! debug
        ! write(PMF_DEBUG,*) 'get_deviation:',value1,value2,get_deviation
    end if

end function get_deviation

end program
