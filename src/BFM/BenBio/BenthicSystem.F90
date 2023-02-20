!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BenthicSystem
!
! DESCRIPTION
!   This is the top level of the benthic system :
!   1) initialise benthic system global variables
!   2) solve biological processes affecting the benthic dynamics 
!      in a specified sequence according to the calculation flags.
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
#include "INCLUDE.h"
!
! INTERFACE
  subroutine BenthicSystemDynamics
!
! USES
  use global_mem, ONLY:RLEN, ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY:  ppY1c, ppY1n, ppY1p, ppY2c, ppY2n, ppY2p, ppY4c, ppY4n,  &
    ppY4p, ppY5c, ppY5n, ppY5p, ppH1c, ppH1n, ppH1p, ppH2c, ppH2n, ppH2p,  &
    iiY1, iiY2, iiY3, iiY4, iiY5, iiH1, iiH2, iiBen, iiPel, flux,          &
    rrBTo, rrATo, reBTn, reBTp, reATn, reATp,                              &
    jbotO2o, jbotN1p, jbotN3n, jbotN4n, jbotN5s, jbotN6r
#ifdef INCLUDE_BENCO2
   use mem, ONLY:  jbotO3h,jbotO3c
#endif
#endif
  use mem_Param, ONLY: CalcBenOrganisms, CalcBenBacteria

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Initialize benthic rates 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrBTo(:)  =   ZERO  ! mgO2/m2 # Total Benthic oxic respiration
  reBTn(:)  =   ZERO  ! mmN/m2  # Total Benthic oxic N mineralization
  reBTp(:)  =   ZERO  ! mmP/m2  # Total Benthic oxic P mineralization
  rrATo(:)  =   ZERO  ! mgO2/m2 # Total Benthic anoxic respiration
  reATn(:)  =   ZERO  ! mmN/m2  # Total Benthic anoxic N mineralization
  reATp(:)  =   ZERO  ! mmP/m2  # Total Benthic anoxic P mineralization

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Initialize here inorganic nutrient fluxes that are given back by
  ! Filterfeeders and nutrient regeration model
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#ifdef INCLUDE_BENCO2
  jbotO3h(:)=ZERO
  jbotO3c(:)=ZERO
#endif
  jbotO2o(:)=ZERO
  jbotN1p(:)=ZERO
  jbotN3n(:)=ZERO
  jbotN4n(:)=ZERO
  jbotN5s(:)=ZERO
  jbotN6r(:)=ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of bioturbation (turenh) and bioirrigation (irrenh)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call BioturbationDynamics

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of biological dynamics 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( CalcBenOrganisms(iiY1)) then
    call BenOrganismsDynamics( iiY1, ppY1c, ppY1n, ppY1p)
  end if

  if ( CalcBenOrganisms(iiY2)) then
    call BenOrganismsDynamics( iiY2, ppY2c, ppY2n, ppY2p)
  end if

  if ( CalcBenOrganisms(iiY4)) then
    call BenOrganismsDynamics( iiY4, ppY4c, ppY4n, ppY4p)
  end if

  if ( CalcBenOrganisms(iiY5)) then
    call BenOrganismsDynamics( iiY5, ppY5c, ppY5n, ppY5p)
  end if

  if ( CalcBenOrganisms(iiY3)) then
    call FilterFeederDynamics
  end if

  if ( CalcBenBacteria(iiH1)) then
    call BenBacDynamics( iiH1, ppH1c, ppH1n, ppH1p)
  end if

  if ( CalcBenBacteria(iiH2)) then
    call BenBacDynamics( iiH2, ppH2c, ppH2n, ppH2p)
  end if

  end subroutine BenthicSystemDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
