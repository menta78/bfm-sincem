#include "INCLUDE.h"
      subroutine env_forcing_pom_bfm_1d 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   BFM1D-Version 0.9                !
! Marco & Giulia
!   Momme Butenschoen, January 2005  !
!   Dipartimento di Fisica           !
!   Universita' di Bologna           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------------------
!  Calculation of physical forcing functions
!-----------------------------------------------------------------------
      use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
      use Mem, ONLY: EIR
#endif
!      use param_mem, ONLY:CalcEcologyFlag
!      use TimeModule, ONLY: time,outdelt,start
       use POM
      implicit none
!
!    PASS PHYSICAL VARIABLES INTO BFM
!
         call pom_to_bfm
!
!  Calculation of vertical extinction coefficient
!  (input is ESSNO, PHYTOPLANKTON, DETRITUS AND BACTERIA BIOMASS)
!
         call CalcVerticalExtinction
!
!  Calculation of the irradiation (forcing function)
!  (input is xEPSNO: vertical extinction)
!         write(6,*) 'Entering CalcLightDistribution...' 
!
         call CalcLightDistribution
 
!      end if
      return
      end

