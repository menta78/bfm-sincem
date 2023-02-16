#include "INCLUDE.h"
!
! **************************************************************
! **************************************************************
! **                                                          **
! ** ONE-DIMENSIONAL BFM-POM  MODELING SYSTEM (BFM-POM1D)     **
! **                                                          **
! ** The modeling system originate from the direct on-line    **
! ** coupling of the 1D Version of the Princeton Ocean model  **
! ** "POM" and the Biological Flux Model "BFM".               **
! **                                                          **
! ** The whole modelling system and its documentation are     **
! ** available for download from the BFM web site:            **
! **                                                          **
! **                  bfm-community.eu                        **
! **                                                          **
! ** For questions and/or information please address to the   **
! ** BFM system team:                                         **
! **                                                          **
! **                 (bfm_st@lists.cmcc.it)                   **
! **                                                          **
! ** Version 1.0 2016                                         **
! **                                                          **
! ** This release has been finalised by Marco Zavatarelli,    **
! ** Giulia Mussap and Nadia Pinardi. However, previous       **
! ** significant contributions were provided also by          **
! ** Momme Butenschoen and Marcello Vichi.                    **
! ** Thanks are due to Prof. George L. Mellor that allowed us **
! ** to modify, use and distribute the one dimensional        **
! ** version of the Princeton Ocean Model.                    **
! **                                                          **
! **                            Marco.Zavatarelli@unibo.it    **
! **                                                          **
! ** This program is free software; you can redistribute it   **
! ** and/or modify it under the terms of the GNU General      **
! ** Public License as published by the Free Software         **
! ** Foundation.                                              **
! ** This program is distributed in the hope that it will be  **
! ** useful,but WITHOUT ANY WARRANTY; without even the        **
! ** implied warranty of  MERCHANTEABILITY or FITNESS FOR A   **
! ** PARTICULAR PURPOSE.  See the GNU General Public License  **
! ** for more details.                                        **
! ** A copy of the GNU General Public License is available at **
! ** http://www.gnu.org/copyleft/gpl.html or by writing to    **
! ** the Free Software Foundation, Inc. 59 Temple Place,      **
! ** Suite 330, Boston, MA 02111, USA.                        **
! **                                                          **
! **************************************************************
! **************************************************************
!
! !ROUTINE: pom_to_bfm
!
! DESCRIPTION    
!    This subroutine passes the physical variables to the BFM
!
! !INTERFACE
!
   subroutine pom_to_bfm
!
! DESCRIPTION
!
!    This subroutine transfers information about the physical anvironment 
!    to BFM. 
!    Data transferred are:
!    Temperature [TB(k)-->ETW(k)]
!    Salinity [SB(k)-->ESW(k)]
!    Density [(RHO(k)*1000)+1000-->ERHO]
!    Water Pressure rho(k)*grav*zz(k)*H*1.e-5 --> EPR(k)
!    Inorganic suspended matter [ISM(k)-->ESS(k)]
!    Surface short wave radiation Recomputed to W/m2)
!    Wind velocity [EWIND, m/s, recomputed from wind stress WUSURF, WVSURF]
!    Bottom stress [ETAUB computed from WUBOT, WVBOT]
!
!*********************************************************************
!
!  -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem,ONLY: RLEN, ONE
!
#ifdef NOPOINTERS
!
      use mem
!
#else
!
      use Mem, ONLY: ETW,  &
                     ESW,  &
                     EIR,  &
                     ESS,  &
                     ERHO, &
                     EWIND,&
                     EPR
!
#endif
!
      use POM, ONLY: KB,     &
                     rcp,    &
                     TB,     &
                     SB,     &
                     RHO,    &
                     SWRAD,  &
                     WUSURF, &
                     WVSURF, &
                     WUBOT,  &
                     WVBOT,  &
                     RHO0,   &
                     RHOSEA, &
                     GRAV,   &
                     H,      &
                     ZZ
!
      use CPL_VARIABLES, ONLY: ISM
!
!  -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!  -----APPROXIMATE DRAG COEFFICIENT FOR WIND VELOCITY COMPUTATION-----
!
   real(RLEN), parameter :: Cd=1.40E-3
!
!  -----REFERENCE AIR DENSITY FOR WIND VELOCITY COMPUTATION-----
!
   real(RLEN), parameter :: rhoa=1.25
!
!  -----LOCAL SCALARS-----
!
      real(RLEN) :: tauw, taub
!
!  -----LOOP COUNTER-----
!
      integer :: k
!
!  ******************************
!  ******************************
!  **                          **
!  ** START DATA TRANSFER      **
!  **                          **
!  ******************************
!  ******************************
!
!  -----LOOP OVER GRID POINT-----
!
      do k = 1 , KB - 1
!
!            --TRANSFER TEMPERATURE--
!
             ETW(k)   = tb(k)
!
!            --TRANSFER SALINITY--
!
             ESW(k)   = sb(k)
!
!            --TRANSFER DENSITY--
!
             ERHO(k)  = (rho(k)*RHO0)+RHO0
!           
!            --COMPUTE PRESSURE--
!
             EPR(k)=-ERHO(k)*GRAV*ZZ(k)*H*1.e-5_RLEN
!
!            ---TRANSFER INORGANIC SUSPENDED MATTER---
!
             ESS(k)   = ISM(k)
!
      end do
!
!  -----TRANSFER URFACE SHORTWAVE RADIATION (RECOMPUTED TO W/M2)-----
!
      EIR(1) = (-ONE)*SWRAD*rcp
!
!  -----COMPUTE WIND VELOCITY (APPROXIMATE CALCULATION FROM WIND STRESS)-----
!
      tauw = sqrt(wusurf**2+wvsurf**2)*rho0
      EWIND = sqrt(tauw/(rhoa*Cd))
!
!  -----TRANSFER BOTTOM STRESS (Nm/s)-----
!
      taub=sqrt(WUBOT**2+WVBOT**2)*rho0
      ETAUB = taub
!
      return
!
      end subroutine pom_to_bfm
!
