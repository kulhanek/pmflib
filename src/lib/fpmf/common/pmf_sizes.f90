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

module pmf_sizes

implicit none

!===============================================================================

integer,parameter       :: PMFDP            = 8         ! double precision
integer,parameter       :: PMFSP            = 4         ! float
integer,parameter       :: PMF_MAX_PATH     = 1024      ! max file name
integer,parameter       :: PMF_MAX_SUNIT    = 50
integer,parameter       :: PMF_MAX_CV_NAME  = 50
integer,parameter       :: PMF_MAX_TYPE     = 10
integer,parameter       :: PMF_MAX_MODE     = 1
integer,parameter       :: PMF_MAX_KEY      = 20
integer,parameter       :: PMF_KEYLINE      = 80

! ------------------------------------------------------------------------------

logical                 :: ignored_arg__    = .false.   ! it is used in disabling
                                                        ! unused variable warning
!===============================================================================

end module pmf_sizes
