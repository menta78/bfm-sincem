!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: CorrectConcNearBed
!
! DESCRIPTION
!   Near bottom correction 
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
      subroutine CorrectConcNearBed(depthlayer, sedi, fto, p_max, wf, correction)
! USES
     use constants, only: RLEN, SEC_PER_DAY
     use mem_Param, only: p_small
     use global_mem, only: ONE,ZERO
#ifdef NOPOINTERS
     use mem
#else
     use mem,         ONLY: ETAUB,NO_BOXES_XY
#endif

      IMPLICIT NONE
 
     ! INPUT
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      real(RLEN), intent(IN),dimension(NO_BOXES_XY)   ::DepthLayer
      real(RLEN), intent(IN),dimension(NO_BOXES_XY)   ::Sedi
      real(RLEN), intent(IN)                          ::fto
      real(RLEN), intent(IN)                          ::p_max
      real(RLEN), intent(IN),dimension(NO_BOXES_XY)   ::wf          ! volumefiltered*Y3c (m/d)
     ! OUTPUT
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      real(RLEN),dimension(NO_BOXES_XY),intent(OUT)   ::correction

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Local Variables
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      real(RLEN),dimension(NO_BOXES_XY)                ::b
      real(RLEN),dimension(NO_BOXES_XY)                ::f
      real(RLEN),dimension(NO_BOXES_XY)                ::r
      real(RLEN),parameter                             ::kappa=1.0E-5_RLEN
    
      f=  min(fto,DepthLayer)
      b = min(100.0_RLEN,Sedi/(p_small+SEC_PER_DAY*kappa*ETAUB(:))) 

      where (b < ONE)
        r=ONE/(ONE-b)*(0.5*DepthLayer)**b
        correction=r*(f**(ONE-b) -p_small**(ONE-b))/(f-p_small)
      elsewhere
        correction=ONE
      endwhere

      return

      end subroutine CorrectConcNearBed

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
