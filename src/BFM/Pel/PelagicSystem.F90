!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PelagicSystem
!
! DESCRIPTION
!   This is the Pelagic Submodel. 
!   All the pelagic biogeochemical modules are called in sequence
!   according to the logical switches
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
#include "DEBUG.h"
!
! INTERFACE
  subroutine PelagicSystemDynamics
!
! USES
  use global_mem, ONLY:RLEN, ZERO
  use mem, ONLY: iiPhytoPlankton, iiMesoZooPlankton, iiMicroZooPlankton,      &
                 iiPelBacteria, flPTN6r, sediR2, sediR3, sediR6, sediPPY 
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF
#endif
  use mem_Param, ONLY: CalcPhytoPlankton,CalcMicroZooPlankton,                &
    CalcMesoZooPlankton, CalcPelBacteria, CalcPelChemistry, ChlDynamicsFlag
  use global_interface,   ONLY: CalcChlorophylla, CalcOxygenSaturation
  use global_interface, ONLY: PhotoAvailableRadiation, &
    PhytoDynamics, LightAdaptationDynamics, MesoZooDynamics, MicroZooDynamics
  use api_bfm, ONLY: LOGUNIT, BOTindices
  use init_var_bfm_local, only: upd_organic_quotas
  use mem_PelSinkSet, ONLY: p_rR6m, p_rR3m, p_rPIm, p_burvel_R2, p_burvel_R3, p_burvel_R6, p_burvel_PI
#if defined INCLUDE_PELCO2
  use mem,            ONLY: sediO5
  use mem_PelSinkSet, ONLY: p_rO5m, p_burvel_O5
#endif
  use mem_PelBac,     ONLY: p_version

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  flPTN6r(:)  =   ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Diagnostic chlorophyll
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call  CalcChlorophylla( )

  !---------------------------------------------
  ! Update quotas of non- and living organic components
  !---------------------------------------------
  call upd_organic_quotas()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set background sedimentation velocities
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  sediR2(:) = ZERO
  sediR3(:) = p_rR3m
  sediR6(:) = p_rR6m
#if defined INCLUDE_PELCO2
  sediO5(:) = p_rO5m
#endif
  do i = 1 , (iiPhytoPlankton)
    sediPPY(:,i) = p_rPIm( i)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissovled oxygen saturation and air-sea flux
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call PelOxygen

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! This part is executed if Optimal Irradiance is used
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if ( ChlDynamicsFlag== 1) then
     do i =1,iiPhytoPlankton
        if ( CalcPhytoPlankton(i)) call PhotoAvailableRadiation( i )
     end do
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute phytoplankton dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do i =1,iiPhytoPlankton
     if ( CalcPhytoPlankton(i)) call PhytoDynamics( i )
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! This part is executed if Optimal Irradiance is used
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if ( ChlDynamicsFlag== 1) then
     do i =1,iiPhytoPlankton
        if ( CalcPhytoPlankton(i)) call LightAdaptationDynamics( i )
     end do
  end if
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute MesoZooPlankton dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do i =1,iiMesoZooPlankton
     if ( CalcMesoZooPlankton(i)) call MesoZooDynamics( i )
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute MicroZooPlankton dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do i =1,iiMicroZooPlankton
     if ( CalcMicroZooPlankton(i)) call MicroZooDynamics( i )
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute Bacteria dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do i =1,iiPelBacteria
     if ( CalcPelBacteria(i) ) then
        if ( p_version(i) .eq. 4) then
           call PelBacDynamics4( i )
        else
           call PelBacDynamics( i )
        endif
     endif
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Prescribe burial velocities at sediment-water column interface
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if ( BOTindices(1) .NE.0 )  then 
      sediR2(BOTindices) = p_burvel_R2
      sediR3(BOTindices) = p_burvel_R3
      sediR6(BOTindices) = p_burvel_R6
#if defined INCLUDE_PELCO2
      sediO5(BOTindices) = p_burvel_O5
#endif
      do i = 1 , ( iiPhytoPlankton)
          sediPPY(BOTindices,i)  =   p_burvel_PI
      end do
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute HydroChemistry (including CO2 in seawater)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if ( CalcPelChemistry) then
    call PelChemDynamics
  end if

  end subroutine PelagicSystemDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
