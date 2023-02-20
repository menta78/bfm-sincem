!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BFM1D_NO_BOXES
!
! DESCRIPTION
!   Initializing procedure
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

#include "cppdefs.h"
!
subroutine BFM1D_NO_BOXES(N,X,Y,Z,XY)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN, bfm_lwp, LOGUNIT
  use mem, ONLY:NO_BOXES,NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z,NO_BOXES_XY, NO_STATES, &
       NO_D3_BOX_STATES
  use api_bfm, ONLY: SRFindices,BOTindices

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, intent(in) :: N,X,Y,Z,XY


  NO_BOXES   = N  
  NO_BOXES_X = X
  NO_BOXES_Y = Y
  NO_BOXES_Z = Z
  NO_BOXES_XY = XY
  NO_STATES   = NO_D3_BOX_STATES * NO_BOXES
  LEVEL1 'NO_BOXES_TOT', N
  LEVEL1 'NO_BOXES_X', X
  LEVEL1 'NO_BOXES_Y', Y
  LEVEL1 'NO_BOXES_Z', Z
  LEVEL1 'NO_BOXES_XY', XY
  LEVEL1 'allocating NO_BOXES_XY'
  allocate(SRFindices(NO_BOXES_XY))
  allocate(BOTindices(NO_BOXES_XY))
  
  
  LEVEL1 'NO_STATES', NO_STATES  
  LEVEL1 'NO_D3_BOX_STATES',NO_D3_BOX_STATES 

end subroutine BFM1D_NO_BOXES
