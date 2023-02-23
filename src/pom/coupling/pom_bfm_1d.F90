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
! ** Subsequent maintance by Lorenzo Mentaschi.               **
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
! !ROUTINE: pom_bfm
!
! !INTERFACE:

       subroutine pom_bfm_1d
!
! !DESCRIPTION:
!  This subroutine handles the BFM-POM1D coupling by:
!    -Passing to BFM Information about the physical environment
!    -Computing the Biogeochemical rates of change (BFM core)
!    -Computing the Physical rates of change for The BFM State Var's
!    -Integrating forward in time BFM state Var's with Source splitting
!     method and leapFrog numerical scheme.
!    -Handling the model output
!    -Reset the BFM state Var's trend arrays at the end of each iteration
!
!*******************************************************************
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
       use api_bfm,ONLY           :out_delta
       use CPL_VARIABLES,ONLY           :savef
       use constants,ONLY         :SEC_PER_DAY
       use POM,ONLY               :time,time0,dti,intt,ilong
!
       use global_mem, ONLY       :RLEN
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
       IMPLICIT NONE
!
!      -----MODEL TIME IN DAYS-----
!
       real(RLEN),save       :: TT
!
!      -----FLAG FOR TT INITIALISATION-----
!
       logical, save          ::first
       data first /.true./
!
!      ********************************************
!      ********************************************
!      **                                        **
!      ** Pass physical variables into BFM core  **
!      ** Compute extinction coefficients        **
!      ** compute vertical light distribution    **
!      **                                        **
!      ********************************************
!      ********************************************

       call pom_env_forcing_1d
!
!      -----EXECUTE BFM CORE-----
!
       call EcologyDynamics
!
!      *******************************************************************
!      *******************************************************************
!      **                                                               **
!      ** Physical transport (Vertical diffusion and sinking) of        **
!      ** BFM state Var's.                                              **
!      ** Integration of BFM pelagic state Var's with Source Splitting  **
!      ** method and Leapfrog numerical scheme                          **
!      **                                                               **
!      *******************************************************************
!      *******************************************************************
!
       call VERT_INTEGRATION
!
       call BENTHIC_TIME_INTEGRATION
!
!
!      -----DEFINE AND UPDATE TIME FOR OUTPUT WRITING-----
!
       if(first) then
!
           TT=time-time0-(dti/SEC_PER_DAY)
           out_delta=savef
           first=.false.
!
       endif
!
       TT = TT + dti/SEC_PER_DAY
!
!      -----MANAGE OUTPUT-----
!
       call pom_dia_bfm(TT)
!
!      -----RESET BFM STATE VAR'S TREND ARRAYS-----
!
       call ResetFluxes
!
       return
!
       end subroutine pom_bfm_1d

!EOC


