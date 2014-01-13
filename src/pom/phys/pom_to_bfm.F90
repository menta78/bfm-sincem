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
!      REAL(RLEN) :: ISM
!
!     -----1D ARRAYS FOR BFM-----
!
      do k = 1 , KB - 1
             ETW(k) = tb(k)
             ESW(k) = sb(k)
             ERHO(k) = (rho(k)*1.E3)+1025.
             ESS(k) = ISM(k)
      end do

      Depth(1) = dz(1)
      do k = 2, KB-1
             Depth(k) = dz(k)+depth(k-1)
      end do
      Depth=Depth*h
!
!
!    Surface radiation in deg.C transformed in W.m-2
!
      EIR(1) = (-1)*SWRAD*rcp
!
!    Wind approximate velocity calculation

       tauw=sqrt(wusurf**2+wvsurf**2)*1.e3
       ewind=sqrt(tauw/(1.25*.0014))
!
      return
      end

