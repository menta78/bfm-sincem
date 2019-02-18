#include "cppdefs.h"

#if defined INCLUDE_PELCO2 || defined INCLUDE_BENCO2
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModulePelCO2
!
! DESCRIPTION
! Module to simulate the Carbonate system dynamics in Pelagic compartment
!
! !INTERFACE
  module mem_CO2
!
! !USES:
  use global_mem
  use SystemForcing, only :ForcingName, ForcingField, FieldInit, FieldClose
  IMPLICIT NONE
!  
!
! !AUTHORS
!   T. Lovato (CMCC) 2017
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2017 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !NAMELIST PelCO2_parameters
  !-------------------------------------------------------------------------!
  ! CARBONATE SYSYEM SETTING
  ! NAME            [UNIT]/KIND             DESCRIPTION
  ! AtmCO20         [ppmv]           Initial atmospheric concentration of CO2
  ! calcAtmpCO2     logical          Compute the partial pressure of Atmospheric CO2
  ! CalcBioAlk      logical          Compute biological processes corrections on total alkalinity
  ! CO2fluxfac      real             Multipling factor for CO2 flux to accelerate air-sea exchange
  !              ---------  SolveSAPHE parameters  -----------
  ! MaxIterPHsolver integer          Maximum number of iterations (default 50)
  !              ---------  Parameters for calcium and calcite ---------
  ! p_kdca          [d-1]            Calcite dissolution rate constant
  ! p_nomega        [-]              Order of the dissolution rate dependence on Omega
  !              ---------  EXTERNAL DATA INPUT STRUCTURES -----------
  ! AtmCO2_N       structure        Read external data for atmospheric CO2 values
  ! AtmSLP_N       structure        Read external data for atmospheric sea level pressure
  ! Example of general input structure for the data structure:
  !          ! Read  !   File                               ! NetCDF  !  Var    !
  !          ! Input !   name                               ! Logical !  name   !
  !AtmCO2_N  =    0  , 'CMIP5_Historical_GHG_1765_2005.dat' , .FALSE.  , 'CO2'  ,
  !          !  RefTime          ! Input      !   Time   !
  !          !  yyyymmdd         ! Frequency  !  interp  !
  !           '1764-07-01 00:00' ,  'yearly'  ,  .TRUE.
  !
  ! Convention for Input reading : 0 = use constant value (default if struct is not initialized)
  !                                1 = read timeseries file ( e.g. CO2 mixing ratios)
  !                                2 = read 2D fields using NEMO fldread 
  ! NOTE: The file "CMIP5_Historical_GHG_1765_2005.dat" is located in "$BFMDIR/tools" folder
  !-----------------------------------------------------------------------------------!
   real(RLEN)           :: AtmCO20 = 365.0_RLEN ! ppm 
   integer              :: MaxIterPHsolver = 50
   real(RLEN)           :: p_kdca
   integer              :: p_nomega
   logical              :: CalcBioAlk = .FALSE.
   real(RLEN)           :: Co2fluxfac = 1.0_RLEN
   type(ForcingName)    :: AtmCO2_N, AtmSLP_N
   type(ForcingField)   :: AtmCO2, AtmSLP
   ! ancillary
   real(RLEN),allocatable,dimension(:) :: patm3d ! atm. pressure over NO_BOXES
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
   public InitCO2, CloseCO2

  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitCO2()

#ifdef NOPOINTERS
  use mem
#else
  use mem,          ONLY: NO_BOXES, NO_BOXES_XY, pH
#endif
  use api_bfm, ONLY: bfm_init
  use mem_Param, ONLY: p_atm0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    namelist /CSYS_parameters/ AtmCO20, MaxIterPHsolver,         &
                            p_kdca, p_nomega,                  &
                            AtmCO2_N, AtmSLP_N, CalcBioAlk, Co2fluxfac
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer            ::error=0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  allocate (patm3d(NO_BOXES))
  patm3d = p_atm0
  !---------------------------------------------------------------------------
  ! Initialize the structured array that defines if a variable is initialized  
  ! with external data.
  !---------------------------------------------------------------------------
                        ! Read  !   File     ! Netcdf  !  Var   ! File    ! Input      !   Time   !
                        ! Input !   name     ! Logical !  name  ! RefTime ! Frequency  !  interp  !
    AtmCO2_N = ForcingName( 0  , "dummy.nc" , .TRUE.  ,"AtmCO2" , "dummy" ,  "dummy"   ,  .TRUE.  )
    AtmSLP_N = ForcingName( 0  , "dummy.nc" , .TRUE.  ,"AtmSLP" , "dummy" ,  "dummy"   ,  .TRUE.  )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    !LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
    !LEVEL1 ' '
    LEVEL1 '     INITIALIZE PELAGIC CARBONATE SYSTEM       ' 
    !LEVEL1 ' '
    LEVEL2 'Namelist content:'
    open(NMLUNIT,file='Carbonate_Dynamics.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=CSYS_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,nml=CSYS_parameters)
    !LEVEL1 ' '
 
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set initial conditions
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Atmospheric CO2 concentration
    AtmCO2%init = AtmCO2_N%init 
    CALL FieldInit(AtmCO2_N, AtmCO2)
    SELECT CASE ( AtmCO2%init )
       CASE (0) ! Costant
          ! the following check is needed to avoid allocation of empty arrays with MPI and land domains
          if (NO_BOXES_XY > 0) then
             AtmCO2%fnow = AtmCO20
             write(LOGUNIT,*) 'Using constant atmospheric CO2 concentration:', AtmCO2%fnow(1)
             write(LOGUNIT,*) ' '
          end if
       CASE (1) ! read external 0-D timeseries
          ! the following check is needed to avoid allocation of empty arrays with MPI and land domains
          if (NO_BOXES_XY > 0) then
             write(LOGUNIT,*) 'BFM Read atmospheric CO2 timeseries. Initial value:', AtmCO2%fnow(1)
             write(LOGUNIT,*) ' '
          end if
       CASE (2) ! read external 2-D fields with NEMO fldread (see envforcing_bfm)
          ! the following check is needed to avoid allocation of empty arrays with MPI and land domains
          if (NO_BOXES_XY > 0) then
             AtmCO2%fnow = AtmCO20
             write(LOGUNIT,*) 'Read CO2 2D fields with NEMO fldread. Initialize with default uniform value', AtmCO2%fnow(1)
             write(LOGUNIT,*) ' '
          end if
    END SELECT

    ! Atmospheric Sea Level Pressure
    AtmSLP%init = AtmSLP_N%init
    CALL FieldInit(AtmSLP_N, AtmSLP)
    if (AtmSLP%init == 0) then
       ! Use constant
       ! the following check is needed to avoid allocation of empty arrays with MPI and land domains
       if (NO_BOXES_XY > 0) then
          AtmSLP%fnow = p_atm0
          write(LOGUNIT,*) 'Using constant atmospheric SLP (see p_atm0 in BFM_General.nml): ', AtmSLP%fnow(1)
          write(LOGUNIT,*) ' '
       end if
    else
      if (AtmSLP%init .eq. 1 ) &
         write(LOGUNIT,*) 'BFM reads atmospheric SLP timeseries from file: ', AtmSLP_N%filename
      if (AtmSLP%init .eq. 2 ) &
         write(LOGUNIT,*) 'Read SLP 2D fields with NEMO fldread. Initialize with default uniform value', AtmSLP%fnow(1)
      write(LOGUNIT,*) ' '
    endif

    ! summary of input parameters
    write(LOGUNIT,*) ' Model uses PH Total Scale '
 
    ! If cold start, CarbonateSystem computes initial pH
    if (bfm_init == 0 ) pH(:) = -ONE

    !LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
    !LEVEL1 ' '

    FLUSH(LOGUNIT)
    return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"ModuleCO2.f90","Carbonate_system.nml")
101 call error_msg_prn(NML_READ,"ModuleCO2.f90","CSYS_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitCO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine CloseCO2()
    implicit none
    
    if (AtmCO2%init == 1 ) then
       ! close external 0-D timeseries
       CALL FieldClose(AtmCO2_N, AtmCO2)
    end if
  end subroutine CloseCO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end module mem_CO2

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
