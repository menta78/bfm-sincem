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
! ** Subsequent maintance by Lorenzo Mentaschi.               **
! ** Thanks are due to Prof. George L. Mellor that allowed us **
! ** to modify use and distribute the one dimensional         **
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
!DESCRIPTION    
! 
!  This subroutine is called from subroutine pom_bfm_1d  and 
!  calls the routines that provides BFM with all the needed information 
!  about the physical environment:
!
!   pom_to_bfm: Transfers the  POM state variables into the correspondent
!               BFM's.
! 
!   CalcVerticalExtinction: computes the light vertical extinction coefficients.
!
!   CalcLightDistribution: define the irradiance vertical profile.
!
! !INTERFACE
!
    SUBROUTINE pom_env_forcing_1d 
!
      IMPLICIT NONE
!
!        -----TRANSFER PHYSICAL VARIABLES INTO BFM-----
!
         call pom_to_bfm
!
!        -----COMPUTE VERTICAL EXTINCTION COEFFICIENTS-----
!
         call CalcVerticalExtinction
!
!        -----COMPUTE THE IRRADIANCE VERTICAL PROFILE-----
!
         call CalcLightDistribution
!
      return
!
      end subroutine pom_env_forcing_1d


