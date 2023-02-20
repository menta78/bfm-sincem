!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: ModuleGlobalMem
!
! DESCRIPTION
!   This module contains global settings:
!   -general constants for controlling prescision,
!   -parameters defining fle streams and error message numbers
!   -the subroutine for printing the message
!   and aborting the simulation
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! INCLUDE
!
! INTERFACE
  MODULE global_mem

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Global Constants
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! default REALS are double precision
  integer, parameter :: sp = selected_real_kind(6, 37)  ! real 4 (6 digits)
  integer, parameter :: dp = selected_real_kind(12,307) ! real 8 (12 digits)
  integer, parameter :: RLEN = dp                       ! default
  real(RLEN), parameter ::ZERO=0.0_RLEN
  real(RLEN), parameter ::ONE=1.0_RLEN
  real(RLEN), parameter ::PI=3.141592653589793_RLEN
  real(RLEN), parameter ::BASETEMP= 10.0_RLEN
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Control logical units for files I/O
  logical               :: bfm_lwp = .FALSE. ! logical writing for proc 0
  integer, parameter    :: bfm_file_FirstUnit = 1000
  integer               :: LOGUNIT
  integer               :: NMLUNIT
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Control BFM coupling with ocean models (skip BFM is NO ocean points)
  logical               :: SkipBFMCore = .FALSE.
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Next parameters are defined to control types of State variables
  integer,    parameter ::OFF=-100
  integer,    parameter ::SINKSOURCE=-1
  integer,    parameter ::NOTRANSPORT=0
  integer,    parameter ::HORTRANSPORT=10
  integer,    parameter ::ALLTRANSPORT=20
  real(RLEN), parameter :: DONE=1._RLEN
  ! Error codes:
  integer,    parameter ::ALLOC=10
  integer,    parameter ::NML_OPEN=11
  integer,    parameter ::NML_READ=12
  integer,    parameter ::DIM_MISMATCH=13
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  contains

  subroutine error_msg_prn(code,infile,what)

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer :: code
  character(LEN=*) :: what,infile
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   write(LOGUNIT,*) "*********** RUN TIME ERROR BEGIN ***********"
   select case (code)
        case (ALLOC)
          write(LOGUNIT,*) "Unable to allocate ",trim(what)," in ", &
                                trim(infile)
        case (NML_OPEN)
          write(LOGUNIT,*) "Unable to open ",trim(what)," in ",  &
                                trim(infile)
        case (NML_READ)
          write(LOGUNIT,*) "Namelist mismatch in ",trim(what),  &
                         " opened by ",trim(infile)
        case (DIM_MISMATCH)
          write(LOGUNIT,*) "Dimension mismatch while reading ",  &
                trim(what)," in ",trim(infile)
    end select
    write(LOGUNIT,*) "***********  RUN TIME ERROR END  ***********"
    stop "BFM error (see logfile)"

  end subroutine error_msg_prn
  
  end module global_mem

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
