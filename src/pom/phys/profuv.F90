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
! !ROUTINE: PROFUV 
!
! !INTERFACE
!
      Subroutine PROFUV(DT2,VELF,SSTRESS,BSTRESS)
!
!DESCRIPTION
!
! This subroutine solves for the vertical profile of the velocity
! zonal and meridional components.
! It handles the surface and bottom boundary conditions.
!
! The routine dummy arguments are:
! DT2: twice the time step
! VELF: The velocity component to be computed
! SSTRESS: The surface stress
! Bstress: The bottom stress.
!
!*****************************************************************
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem,ONLY: RLEN, ZERO
      use POM,ONLY: H,A,C,KM,DZ,DZZ,VH,VHP,UB,VB,UMOL,CBC,KB
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----SCALAR ARGUMENTS-----
!
      REAL(RLEN) :: DT2, SSTRESS,BSTRESS
!
!    -----ARRAY ARGUMENT-----
!
     REAL(RLEN), dimension (KB) :: VELF
!
!     -----LOOP COUNTERS (DOWNWARD AND UPWARD)-----
!
      INTEGER :: K,KI
!
!     -----INTRINSIC FUNCTION-----
!
      INTRINSIC SQRT
!
!    *******************************************************
!    *******************************************************
!    **                                                   **
!    ** THE FOLLOWING SECTION SOLVES THE EQUATION         **
!    ** DT2*(KM*U')' - U= -UB                             **
!    **                                                   **
!    *******************************************************
!    *******************************************************
!
      DO K = 2,KB - 1
!
          A(K-1) = -DT2* (KM(K)+UMOL)/ (DZ(K-1)*DZZ(K-1)*H*H)
          C(K) = -DT2* (KM(K)+UMOL)/ (DZ(K)*DZZ(K-1)*H*H)
!
      END DO
!
      VH(1) = A(1)/ (A(1)-1.)
!
      VHP(1) = (-DT2*SSTRESS/ (-DZ(1)*H)-VELF(1))/ (A(1)-1.)
!
      DO K = 2,KB - 2
!
          VHP(K) = 1./ (A(K)+C(K)* (1.-VH(K-1))-1.)
          VH(K) = A(K)*VHP(K)
          VHP(K) = (C(K)*VHP(K-1)-VELF(K))*VHP(K)
!
      END DO
!
      BSTRESS = CBC*SQRT(UB(KB-1)**2+VB(KB-1)**2)
!
      VELF(KB-1) = (C(KB-1)*VHP(KB-2)-VELF(KB-1))/ &
                   (BSTRESS*DT2/                   &
                   (-DZ(KB-1)*H)-1.- (VH(KB-2)-1.)*C(KB-1))
!
      DO K = 2,KB - 1
!
          KI = KB - K
          VELF(KI) = VH(KI)*VELF(KI+1) + VHP(KI)
!
      END DO
!
!     -----COMPUTE NEW BOTTOM STRESS-----
!
      BSTRESS = -BSTRESS*VELF(KB-1)
!
!     -----HOUSE CLEANING-----
!
      DO K = 1,KB
!
          VH(K)  = ZERO
          VHP(K) = ZERO
          A(K)   = ZERO
          C(K)   = ZERO
      END DO

      RETURN
!
      end subroutine PROFUV
