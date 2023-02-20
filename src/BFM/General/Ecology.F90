!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: Ecology
!
! DESCRIPTION
! This is the main interface to call all the other sub-models of the BFM
! Depending on the choices of macros and parameters this part activates:
! - the sea ice model
! - the pelagic model
! - the benthic model
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
  subroutine EcologyDynamics
!
! USES
  use global_mem, ONLY:RLEN
  use mem,  ONLY: iiBen, iiPel, iiReset, flux
  use mem_Param,  ONLY: CalcPelagicFlag, CalcBenthicFlag, CalcConservationFlag
  use api_bfm, ONLY: LOGUNIT
#ifdef INCLUDE_SEAICE
  use mem_Param,  ONLY: CalcSeaiceFlag
#endif

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#ifdef DEBUG
  call  flux(1,iiReset,1,1,0.0)
#endif

#ifdef INCLUDE_SEAICE
  if ( CalcSeaiceFlag) then

    call SeaiceSystemDynamics

  end if
#endif

  if ( CalcPelagicFlag) then

    call PelagicSystemDynamics

  end if

  if ( CalcBenthicFlag ) then
#if defined BENTHIC_BIO
     ! Intermediate benthic return
     call PelForcingForBenDynamics
     call BenthicSystemDynamics
     call BenthicReminDynamics
     call BenOxygenDynamics

#elif defined BENTHIC_FULL
     ! Full benthic nutrients
     call PelForcingForBenDynamics
     call BenthicSystemDynamics
     call BenthicNutrientDynamics

#else
     ! Simple benthic return
     call BenthicReturnDynamics
#endif
  endif
  
  call PelagicBenthicCoupling

  if (CalcConservationFlag) &
     call CheckMassConservationDynamics

  end subroutine EcologyDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
