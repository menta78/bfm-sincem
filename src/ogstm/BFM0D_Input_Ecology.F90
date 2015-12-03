#include"cppdefs.h"
!
subroutine BFM0D_Input_EcologyDynamics(sur,bot,BFM0D_trn,dim_BFM0D_trn,BFM0D_er)
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

  logical   , intent(in) :: sur, bot
  integer dim_BFM0D_trn
  real(RLEN), intent(in) :: BFM0D_trn(dim_BFM0D_trn), BFM0D_er(10)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer i

    IF(sur) then
       SRFindices(1) = 1
    ELSE
       SRFindices(1) = 0
    ENDIF


    IF(bot) then
       BOTindices(1) = 1
    ELSE
       BOTindices(1) = 0
    ENDIF


  DO i=1,dim_BFM0D_trn
      D3STATE(i,1)=BFM0D_trn(i)
!       LEVEL1 'BFM0D_Input_EcologyDynamics:D3STATE(',i,')=',D3STATE(i,1) 
  END DO

  ! Environmental regulating factors
  ETW    = BFM0D_er(1)
! LEVEL1 'BFM0D_Input_EcologyDynamics:ETW', ETW
  ESW    = BFM0D_er(2)
! LEVEL1 'BFM0D_Input_EcologyDynamics:ESW', ESW
  ERHO   = BFM0D_er(3)
! LEVEL1 'BFM0D_Input_EcologyDynamics:ERHO', ERHO
  EICE   = BFM0D_er(4)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EICE', EICE
#ifdef INCLUDE_PELCO2
  AtmCO2%fnow = BFM0D_er(5)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EPCO2air', EPCO2air
#endif
  EIR    = BFM0D_er(6)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EIR', EIR
  SUNQ   = BFM0D_er(7)
! LEVEL1 'BFM0D_Input_EcologyDynamics:DL', SUNQ
  DEPTH  = BFM0D_er(8)
! LEVEL1 'BFM0D_Input_EcologyDynamics:DEPTH', DEPTH
  EWIND  = BFM0D_er(9)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EWIND,', EWIND
  ph     = BFM0D_er(10)
! LEVEL1 'BFM0D_Input_EcologyDynamics:PH,', ph

end subroutine BFM0D_Input_EcologyDynamics
