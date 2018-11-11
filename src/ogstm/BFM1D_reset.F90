#include "cppdefs.h"
!
subroutine BFM1D_reset
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  use (import) other modules
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, only:ZERO
  use mem

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Start the allocation of pelagic state global
  ! matrix and pointers
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer i, dim_BFM1D_tra

  call ResetFluxes()

  SiPToutr = ZERO
  SiPToutc = ZERO
  SiPTouto = ZERO
  SiPToutn = ZERO
  SiPToutp = ZERO
  SiPTouti = ZERO
  SoO2Airo = ZERO
  SoPTinr  = ZERO
  SoPTinc  = ZERO
  SoPTino  = ZERO
  SoPTinn  = ZERO
  SoPTinp  = ZERO
  SoPTini  = ZERO
  SoRIc    = ZERO
  SoRIn    = ZERO
  SoRIp    = ZERO
  SoRIs    = ZERO

#ifdef INCLUDE_PELCO2
  jsurO3c=0
#endif

end subroutine BFM1D_reset
