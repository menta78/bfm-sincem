! **************************************************************
! **************************************************************
! **                                                          **
! ** ONE-DIMENSIONAL BFM-POM  MODELING SYSTEM (BFM-POM1D)     **
! **                                                          **
! ** The modeling system originate from the direct on-line    **
! ** coupling of the 1D Version of the Princeton Ocean model  **
! ** "POM" and the Biogeochemical Flux Model "BFM".           **
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
! !ROUTINE: Calcdepth
!
!DESCRIPTION    
!
!  This subroutine establishes the vertical resolution with log distributions
!  at the top and bottom, and a linear distribution in between 
!  The number of surface layers with log distribution is given by
!
!  KL1-2
!
!  while the number of bottom layrr with log distribution is given by
!
!  KB-KL2-1
!
!  The values: KL1 = .3*KB AND KL2 = KB-2.
!  Yields a log distribution at the top and none at the bottom.
!
!  All the data needed for the vertical coordinate system deinition reach
!  the rsubroutine via the module POM
!
!*************************************************************************
!
! !INTERFACE
!
      Subroutine CALCDEPTH
!
! -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----

      use global_mem,ONLY: RLEN
!
      use POM,ONLY: KB, Z, ZZ, DZ, DZZ, KL1, KL2
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----LOCAL SCALARS-----
!
      REAL(RLEN) :: BB,CC,DEL1,DEL2
!
!     -----LOOP COUNTER-----
!
      INTEGER :: K
!
!     -----INTRINSIC FUNCTIONS-----
!
      INTRINSIC EXP,FLOAT
!
!     -----START COMPUTATION OF Z AND ZZ-----
!
      BB = FLOAT(KL2-KL1) + 4.
      CC = FLOAT(KL1) - 2.
      DEL1 = 2./BB/EXP(.693147*FLOAT(KL1-2))
      DEL2 = 2./BB/EXP(.693147*FLOAT(KB-KL2-1))
      Z(1) = 0.
      ZZ(1) = -DEL1/2.
!
!     -----SURFACE LOG DISTRIBUTION-----
!
      DO K = 2,KL1 - 2
          Z(K) = -DEL1*EXP(.693147*FLOAT(K-2))
          ZZ(K) = -DEL1*EXP(.693147* (FLOAT(K)-1.5))
      END DO
!
!     -----LINEAR DISTRIBUTION IN THE MIDDLE-----
!
      DO K = KL1 - 1,KL2 + 1
          Z(K) = - (FLOAT(K)-CC)/BB
          ZZ(K) = - (FLOAT(K)-CC+0.5)/BB
      END DO
!
!     -----BOTTOM LOG DISTRIBUTION-----
!
      DO K = KL2 + 1,KB - 1
          Z(K) = (1.0-DEL2*EXP(.693147*FLOAT(KB-K-1)))* (-1.)
          ZZ(K) = (1.0-DEL2*EXP(.693147* (FLOAT(KB-K)-1.5)))* (-1.)
      END DO
!
      Z(KB) = -1.0
      ZZ(KB-1) = -1.* (1.-DEL2/2.)
      ZZ(KB) = -1.* (1.+DEL2/2.)
!
!     -----DEFINE DZ AND DZZ-----
!
      DO K = 1,KB - 1
          DZ(K) = Z(K) - Z(K+1)
          DZZ(K) = ZZ(K) - ZZ(K+1)
      END DO
      RETURN
      
      end subroutine CALCDEPTH
