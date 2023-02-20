!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BenthicReturn
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
  subroutine BenthicReturnDynamics
!
! USES
  use global_mem, ONLY:RLEN, ONE, bfm_lwp
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE_BEN, D2DIAGNOS, D3STATE, D3DIAGNOS
#else
  use mem,  ONLY: Q6c, Q1c, Q6p, Q1p, Q6n, Q1n, Q6s
#endif
  use mem, ONLY: ppQ6c, ppQ1c, ppQ6p, ppQ1p, ppQ6n, ppQ1n, ppQ6s,  &
    jbotO2o, jbotN1p, jbotN3n, jbotN4n, jbotN5s, jbotO3c, jbotO3h, &
    ppO2o, ppN1p, ppN3n, ppN4n, ppN5s, ppO3c, ppO3h, ppETW, &
    NO_BOXES_XY, iiBen, iiPel, flux_vector, ETW, O2o
  use mem_BenthicReturn
  use mem_Param, ONLY: CalcConservationFlag, AssignPelBenFluxesInBFMFlag
  use constants, ONLY: MW_C
  use mem_globalfun, ONLY: MM_POWER
  use api_bfm, ONLY: BOTindices

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer   :: box, kbot
  real(RLEN),dimension(NO_BOXES_XY)  :: rate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if ( .NOT. AssignPelBenFluxesInBFMFlag ) return

  ! remineralization fluxes
  rate  =   p_rmnQ6c* Q6c(:)
  call flux_vector( iiBen, ppQ6c,ppQ6c,-( rate) )
  jbotO2o(:)  = - rate / MW_C
  jbotO3c(:)  = rate

  rate  =   p_rmnQ1c* Q1c(:)
  call flux_vector( iiBen, ppQ1c,ppQ1c,-( rate) )
  jbotO2o(:)  = jbotO2o(:) - rate / MW_C
  jbotO3c(:)  = jbotO3c(:) + rate

  rate  =   p_rmnQ6p* Q6p(:)
  call flux_vector( iiBen, ppQ6p,ppQ6p,-( rate) )
  jbotN1p(:)  = rate

  rate  =   p_rmnQ1p* Q1p(:)
  call flux_vector( iiBen, ppQ1p,ppQ1p,-( rate) )
  jbotN1p(:)  = jbotN1p(:)+ rate

  rate  =   p_rmnQ6n* Q6n(:)
  call flux_vector( iiBen, ppQ6n,ppQ6n,-( rate) )
  jbotN3n(:)  = rate* p_pQIN3
  jbotN4n(:)  = rate*( ONE - p_pQIN3)

  rate  =   p_rmnQ1n* Q1n(:)
  call flux_vector( iiBen, ppQ1n,ppQ1n,-( rate) )
  jbotN3n(:)  = jbotN3n(:)+ rate* p_pQIN3
  jbotN4n(:)  = jbotN4n(:)+ rate*( ONE - p_pQIN3)

  rate  =   p_rmnQ6s* Q6s(:)
  call flux_vector( iiBen, ppQ6s,ppQ6s,-( rate) )
  jbotN5s(:)  = rate

  end subroutine BenthicReturnDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
