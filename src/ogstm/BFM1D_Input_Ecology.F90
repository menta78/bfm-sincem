!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BFM1D_Input_EcologyDynamics
!
! DESCRIPTION
!   Interface to receive inputs from ogstm
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


#include"cppdefs.h"
!
subroutine BFM1D_Input_EcologyDynamics(bot,BFM1D_er)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
  use api_bfm, ONLY: SRFindices, BOTindices
  use mem
  use mem_CO2

  IMPLICIT NONE

  integer   , intent(in) :: bot
  integer dim_BFM1D_trn
  real(RLEN), intent(in) :: BFM1D_er(NO_BOXES,11)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer i,ib
  real(RLEN) bottom


 !      BOTindices(:) = 0
 !      SRFindices(:) = 0
 !      SRFindices(1) = 1
 !      BOTindices(bot)=1
        SRFindices(1) = 1
        BOTindices(1) = bot

!  DO ib=bot+1,NO_BOXES
!     BOTindices(ib)=1
!  ENDDO


!  DO ib=1,NO_BOXES
!  DO i=1,dim_BFM1D_trn
!      D3STATE(ib,i)=BFM1D_trn(ib,i)
!  END DO
!       LEVEL1
!       'BFM1D_Input_EcologyDynamics:D3STATE(',i,')=',D3STATE(i,1) 
!  END DO

  ! Environmental regulating factors
  ETW(:)    = BFM1D_er(:,1)
  ESW(:)    = BFM1D_er(:,2)
  ERHO(:)   = BFM1D_er(:,3)
  EICE(:)   = BFM1D_er(1,4)
#ifdef INCLUDE_PELCO2
  AtmCO2%fnow(:) = BFM1D_er(1,5)
#endif
  EIR(:)    = BFM1D_er(:,6)
  SUNQ(:)   = BFM1D_er(1,7)
  DEPTH(:)  = BFM1D_er(:,8)
  EWIND(:)  = BFM1D_er(1,9)
  ph(:)     = BFM1D_er(:,10)
  BAC_ACT_FACT(:) = BFM1D_er(:,11)

  bottom=0
  DO ib=1,NO_BOXES
     bottom = bottom + DEPTH(ib)
     EPR(ib) = bottom - DEPTH(ib)/2
  ENDDO
  ! LEVEL1 'BFM1D_Input_EcologyDynamics:ETW', ETW
end subroutine BFM1D_Input_EcologyDynamics



subroutine BFMmit_Input_EcologyDynamics(bot,BFM1D_trn,dim_BFM1D_trn,BFM1D_er)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
  use api_bfm, ONLY: SRFindices, BOTindices
  use mem
  use mem_CO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  integer   , intent(in) :: bot
  integer dim_BFM1D_trn
!  real(RLEN), intent(in) :: BFM1D_trn(dim_BFM1D_trn,NO_BOXES),BFM1D_er(NO_BOXES,10)
  real(RLEN), intent(in) :: BFM1D_trn(NO_BOXES,dim_BFM1D_trn),BFM1D_er(NO_BOXES,11)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer i,ib
  real(RLEN) bottom


 !      BOTindices(:) = 0
 !      SRFindices(:) = 0
 !      SRFindices(1) = 1
 !      BOTindices(bot)=1
        SRFindices(1) = 1
        BOTindices(1) = bot

!  DO ib=bot+1,NO_BOXES
!     BOTindices(ib)=1
!  ENDDO


  DO ib=1,NO_BOXES
  DO i=1,dim_BFM1D_trn
      D3STATE(ib,i)=BFM1D_trn(ib,i)
  END DO
!       LEVEL1
!       'BFM1D_Input_EcologyDynamics:D3STATE(',i,')=',D3STATE(i,1) 
  END DO

  ! Environmental regulating factors
  ETW(:)    = BFM1D_er(:,1)
! LEVEL1 'BFM1D_Input_EcologyDynamics:ETW', ETW
  ESW(:)   = BFM1D_er(:,2)
! LEVEL1 'BFM1D_Input_EcologyDynamics:ESW', ESW
  ERHO(:)   = BFM1D_er(:,3)
! LEVEL1 'BFM1D_Input_EcologyDynamics:ERHO', ERHO
  EICE(:)   = BFM1D_er(1,4)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EICE', EICE
#ifdef INCLUDE_PELCO2
  AtmCO2%fnow(:) = BFM1D_er(1,5)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EPCO2air', EPCO2air
#endif
  EIR(:)    = BFM1D_er(:,6)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EIR', EIR
  SUNQ(:)   = BFM1D_er(1,7)
! LEVEL1 'BFM1D_Input_EcologyDynamics:DL', SUNQ
  DEPTH(:)  = BFM1D_er(:,8)
! LEVEL1 'BFM1D_Input_EcologyDynamics:DEPTH', DEPTH
  EWIND(:)  = BFM1D_er(1,9)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EWIND,', EWIND
  ph(:)     = BFM1D_er(:,10)
! LEVEL1 'BFM1D_Input_EcologyDynamics:PH,', ph
  BAC_ACT_FACT(:) = BFM1D_er(:,11)
! LEVEL1 'BFM1D_Input_EcologyDynamics:BAC_ACT_FACT',BAC_ACT_FACT

  bottom=0
  DO ib=1,NO_BOXES
     bottom = bottom + DEPTH(ib)
     EPR(ib) = bottom - DEPTH(ib)/2
  ENDDO
end subroutine BFMmit_Input_EcologyDynamics


