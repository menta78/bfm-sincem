#include "DEBUG.h"
#include "INCLUDE.h"

#ifdef INCLUDE_PELCO2
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! DESCRIPTION
! Dissolution of Particulate Inorganic Carbon (calcite/aragonite) in seawater
!
! !INTERFACE
  subroutine PelPICDynamics()
!
! !USES:

  use global_mem, ONLY:RLEN,ONE,ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: NO_BOXES, OCalc, O5c, ppO5c, ppO3c, ppO3h, iiPel, flux_vector
#endif
  use mem_Param,  ONLY: p_small
  use mem_CO2,    ONLY: p_kdca,p_nomega
  use constants,    only: SEC_PER_DAY,C2ALK
!  
!
! !AUTHORS
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2014 BFM System Team (bfm_st@lists.cmcc.it)
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
!
!
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)           :: excess(NO_BOXES),rdiss(NO_BOXES)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute undersaturation
  excess(:) = max(ZERO,ONE - OCalc(:))
  ! Dissolution rate of C in CaCO3 (mg C/m3/d)
  ! from Morse and Berner (1972)
  rdiss(:) = p_kdca * excess(:)**p_nomega * O5c(:)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Inorganic carbon and alkalinity flux
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call flux_vector( iiPel, ppO5c, ppO3c, rdiss(:) )
  call flux_vector( iiPel, ppO3h, ppO3h, -C2ALK*rdiss(:) )

  end subroutine PelPICDynamics
#endif
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
