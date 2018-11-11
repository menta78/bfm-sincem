#include"cppdefs.h"
!
subroutine BFM1D_Input_EcologyDynamics(bot,BFM1D_trn,dim_BFM1D_trn,BFM1D_er)
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
  real(RLEN), intent(in) :: BFM1D_trn(NO_BOXES,dim_BFM1D_trn),BFM1D_er(NO_BOXES,10)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer i,ib


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
      D3STATE(i,ib)=BFM1D_trn(ib,i)
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
  EICE(:)   = BFM1D_er(:,4)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EICE', EICE
#ifdef INCLUDE_PELCO2
  AtmCO2%fnow(:) = BFM1D_er(:,5)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EPCO2air', EPCO2air
#endif
  EIR(:)    = BFM1D_er(:,6)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EIR', EIR
  SUNQ(:)   = BFM1D_er(:,7)
! LEVEL1 'BFM1D_Input_EcologyDynamics:DL', SUNQ
  DEPTH(:)  = BFM1D_er(:,8)
! LEVEL1 'BFM1D_Input_EcologyDynamics:DEPTH', DEPTH
  EWIND(:)  = BFM1D_er(:,9)
! LEVEL1 'BFM1D_Input_EcologyDynamics:EWIND,', EWIND
  ph(:)     = BFM1D_er(:,10)
! LEVEL1 'BFM1D_Input_EcologyDynamics:PH,', ph

end subroutine BFM1D_Input_EcologyDynamics
