!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: Bioturbation
!
! DESCRIPTION
!   Describes sediment reworking by benthic organisms and their effect on
!   vertical transport of dissolved (bioirrigation) and particulate 
!   (bioturbation) matter
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
  subroutine BioturbationDynamics
!
! USES
  use global_mem, ONLY:RLEN,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: D6m, D7m, D8m, D9m, Y2c, Y5c, Y1c, Y3c, Y4c
  use mem, ONLY: ppD6m, ppD7m, ppD8m, ppD9m, ppY2c, ppY5c, ppY1c, &
    ppY4c, turenh, irrenh, &
    ETW_Ben, NO_BOXES_XY, iiBen, iiPel, flux_vector
#endif
  use mem_Bioturbation
  use mem_globalfun,   ONLY: MM

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: Ytur
  real(RLEN),dimension(NO_BOXES_XY)  :: Yirr
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature Response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   (p_q10)**(( ETW_Ben(:)- 10.0_RLEN)* 0.1_RLEN)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of bioturbation. ''turenh'' is a global box variable!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  Ytur  =   Y2c(:)+ Y5c(:)+ p_turY1* Y1c(:)

  turenh(:)  =   ONE+ p_cmtur* MM( Ytur,  p_chtur)* et

  call flux_vector(iiBen, ppD6m,ppD6m, p_Etur* turenh(:)*( ONE- exp( - &
    p_cturm/ D6m(:)))/ D6m(:))
  call flux_vector(iiBen, ppD7m,ppD7m, p_Etur* turenh(:)*( ONE- exp( - &
    p_cturm/ D7m(:)))/ D7m(:))
  call flux_vector(iiBen, ppD8m,ppD8m, p_Etur* turenh(:)*( ONE- exp( - &
    p_cturm/ D8m(:)))/ D8m(:))
  call flux_vector(iiBen, ppD9m,ppD9m, p_Etur* turenh(:)*( ONE- exp( - &
    p_cturm/ D9m(:)))/ D9m(:))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of bioirrigation. ''irrenh'' is a global box variable!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  Yirr  =   Y2c(:)+ Y5c(:)+ p_irrY4* Y4c(:) +Y3c(:)

  irrenh(:)  =   ONE+ p_cmirr* MM( Yirr,  p_chirr)* et

  end subroutine BioturbationDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
