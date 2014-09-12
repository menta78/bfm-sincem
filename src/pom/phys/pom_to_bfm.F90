#include "INCLUDE.h"
!
!     PASS PHYSICAL VARIABLES TO BFM
!
      subroutine pom_to_bfm
#ifdef NOPOINTERS
      use mem
#else
      use Mem, ONLY: ETW, ESW, EIR, ESS,ERHO, ewind 
#endif
      use POM
      use Service, ONLY: ISM
      implicit none
      real rcp,dep
      real tauw
      parameter (rcp=4.187E6)
      integer :: k,ll
!      REAL(RLEN) :: 
!
!     -----1D ARRAYS FOR BFM-----
!
      do k = 1 , KB - 1
             ETW(k) = tb(k)
             ESW(k) = sb(k)
             ERHO(k) = (rho(k)*1.E3)+1000.
             ESS(k) = ISM(k)
             Depth(k) = dz(k)*h   !   (G)
      end do
!  Surface radiation in deg.C transformed in W.m-2
!
      EIR(1) = (-1)*SWRAD*rcp
!
!    Wind approximate velocity calculation

       tauw=sqrt(wusurf**2+wvsurf**2)*1.e3
       ewind=sqrt(tauw/(1.25*.0014))
!
      return
      end

