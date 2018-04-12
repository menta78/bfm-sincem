#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicReturn
!
! DESCRIPTION
!   This process is a very simple parameterisation of benthic 
!   remineralisation.
!   Benthic organism sub-model cannot be used with this model.
!   A constant portion of the organic matter in the sediments is
!   released to the water column as inorganic nutrients.
!   Oxygen consumption is stoichiometrically associated to carbon 
!   remineralisation rates and nitrogen remineralisation is partitioned 
!   into ammonium and nitrate flux with a constant value. 
!   Fluxes are then used as boundary conditions for pelagic variables
!   in PelagicBenthicCoupling.F90
!
!
! !INTERFACE
  subroutine BenthicReturnDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN, ONE, bfm_lwp
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE_BEN, D2DIAGNOS
#else
  use mem,  ONLY: Q6c, Q1c, Q6p, Q1p, Q6n, Q1n, Q6s
#endif
  use mem, ONLY: ppQ6c, ppQ1c, ppQ6p, ppQ1p, ppQ6n, ppQ1n, ppQ6s,  &
    jbotO2o, jbotN1p, jbotN3n, jbotN4n, jbotN5s, jbotO3c, jbotO3h, &
    ppO2o, ppN1p, ppN3n, ppN4n, ppN5s, ppO3c, ppO3h, &
    NO_BOXES_XY, iiBen, iiPel, flux_vector
  use mem_BenthicReturn
  use mem_Param, ONLY: CalcConservationFlag, AssignPelBenFluxesInBFMFlag
  use constants, ONLY: MW_C
!  
!
! !AUTHORS
!   Piet Ruardij 
!
!
!
! !REVISION_HISTORY
!   Created at Fri Apr 30 21:30:27 CEST 2004
!
!
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, the BFM team
!   (rua@nioz.nl, vichi@bo.ingv.it)
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: rate

  if ( .NOT. AssignPelBenFluxesInBFMFlag ) return

  rate  =   p_reminQ6* retfac* Q6c(:)
  call flux_vector( iiBen, ppQ6c,ppQ6c,-( rate) )
  jbotO2o(:)  = - rate / MW_C
  jbotO3c(:)  = rate

  rate  =   p_reminQ1* retfac* Q1c(:)
  call flux_vector( iiBen, ppQ1c,ppQ1c,-( rate) )
  jbotO2o(:)  = jbotO2o(:) - rate / MW_C
  jbotO3c(:)  = jbotO3c(:) + rate

  rate  =   p_reminQ6* retfac* Q6p(:)
  call flux_vector( iiBen, ppQ6p,ppQ6p,-( rate) )
  jbotN1p(:)  = rate

  rate  =   p_reminQ1* retfac* Q1p(:)
  call flux_vector( iiBen, ppQ1p,ppQ1p,-( rate) )
  jbotN1p(:)  = jbotN1p(:)+ rate

  rate  =   p_reminQ6* retfac* Q6n(:)
  call flux_vector( iiBen, ppQ6n,ppQ6n,-( rate) )
  jbotN3n(:)  = rate* p_pQIN3
  jbotN4n(:)  = rate*( ONE - p_pQIN3)

  rate  =   p_reminQ1* retfac* Q1n(:)
  call flux_vector( iiBen, ppQ1n,ppQ1n,-( rate) )
  jbotN3n(:)  = jbotN3n(:)+ rate* p_pQIN3
  jbotN4n(:)  = jbotN4n(:)+ rate*( ONE - p_pQIN3)

  rate  =   p_reminQ6* retfac* Q6s(:)
  call flux_vector( iiBen, ppQ6s,ppQ6s,-( rate) )
  jbotN5s(:)  = rate

  end subroutine BenthicReturnDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
