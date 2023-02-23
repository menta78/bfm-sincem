!   SCHEMES FOR VERTICAL ADVECTION
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!
! !SCHEMES FOR VERTICAL ADVECTION
!
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

!   --- W(K-1)   ---
!    |            |
!    |  F(K-1)    |  H*DZ(K-1) ---
!    |            |             |
!   --- W(K)     ---            |  H*DZZ(K-1)
!    |            |             |
!    |  F(K)      |  H*DZ(K)   ---
!    |            |             |
!   --- W(K+1)   ---            |  H*DZZ(K)
!                               |
!                              ---
!  de F/de t = - grad(W*F) = - W*grad(F) - F*grad(W)



SUBROUTINE VERT_ADV_UPWIND_SCHEME(F, W, DZ, KB, H, ADVTND)
        ! 
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-= 
        ! Calculate vertical advection. Mind downward velocities are negative! 
        ! Upwind scheme: VELOCITIES CAN BE POSITIVE OR NEGATIVE
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-= 
!    THE DUMMY ARGUMENTS ARE:
!    F:      CONCENTRATION
!    W:      VELOCITY (<0 IF SINKING, IF PROFILES ARE GIVEN CAN BE POSITIVE)
!    DZ:     THICKNESS OF THE CELLS
!    KB:     NUMBER OF LAYERS (THE CELLS ARE KB-1)
!    H:      WATER DEPTH
!    ADVTND: OUTPUT ADVECTION TENDENCY
        !

        IMPLICIT NONE

        INTEGER           :: KB
        REAL, INTENT(IN)  :: F(KB), W(KB), DZ(KB)
        REAL, INTENT(IN)  :: H
        REAL, INTENT(OUT) :: ADVTND(KB)

        INTEGER           :: K
        REAL              :: ZFLUX(KB)

        ZFLUX(1) = 0
        ZFLUX(KB) = 0.e0

        DO K=2,KB-1
           ZFLUX(K)=0.5E0*( (W(K)+ABS(W(K)))*F(K) + (W(K)-ABS(W(K)))*F(K-1) )
        END DO

        DO K=1,KB-1
           ADVTND(K) = 1/(H*DZ(K)) * (ZFLUX(K+1) - ZFLUX(K))
        END DO
        ADVTND(KB) = 0

END SUBROUTINE



SUBROUTINE VERT_ADV_SINK_UPWIND_SCHEME(F, W, DZ, KB, H, ADVTND)
        ! 
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-= 
        ! Calculate vertical advection. Mind downward velocities are negative! 
        ! Upwind scheme: VELOCITY MUST BE NEGATIVE, ONLY FOR SINKING STUFF!!! 
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-= 
        !
!    THE DUMMY ARGUMENTS ARE:
!    F:      CONCENTRATION
!    W:      VELOCITY (HERE CAN BE ONLY A SKINK<0!!!)
!    DZ:     THICKNESS OF THE CELLS
!    KB:     NUMBER OF LAYERS (THE CELLS ARE KB-1)
!    H:      WATER DEPTH
!    ADVTND: OUTPUT ADVECTION TENDENCY

        IMPLICIT NONE

        INTEGER           :: KB
        REAL, INTENT(IN)  :: F(KB), W(KB), DZ(KB)
        REAL, INTENT(IN)  :: H
        REAL, INTENT(OUT) :: ADVTND(KB)

        INTEGER           :: K

        ADVTND(1) = 1/(H*DZ(1))*F(1)*W(2)
        DO K=2,KB-1
           ADVTND(K) = 1/(H*DZ(K)) * (F(K)*W(K+1) - F(K-1)*W(K))
        END DO
        ADVTND(KB) = 0

END SUBROUTINE



SUBROUTINE VERT_ADV_CENTERED_SCHEME(F, W, DZ, DZZ, KB, H, ADVTND)
        !  ! centered scheme: Gordon & Stern 1982
        !  TO BE VERIFIED, DO NOT USE!!
        !   --- W(K-1)   ---
        !    |            |
        !    |  F(K-1)    |  H*DZ(K-1) ---
        !    |            |             |
        !   --- W(K)     ---            |  H*DZZ(K-1)
        !    |            |             |
        !    |  F(K)      |  H*DZ(K)   ---
        !    |            |             |
        !   --- W(K+1)   ---            |  H*DZZ(K)
        !                               |
        !                              ---
        !  de F/de t = - grad(W*F) = - W*grad(F) - F*grad(W)
!    THE DUMMY ARGUMENTS ARE:
!    F:      CONCENTRATION
!    W:      VELOCITY (HERE CAN BE ONLY A SKINK<0!!!)
!    DZ:     THICKNESS OF THE CELLS
!    KB:     NUMBER OF LAYERS (THE CELLS ARE KB-1)
!    H:      WATER DEPTH
!    ADVTND: OUTPUT ADVECTION TENDENCY

        IMPLICIT NONE

        INTEGER           :: KB
        REAL, INTENT(IN)  :: F(KB), W(KB), DZ(KB), DZZ(KB)
        REAL, INTENT(IN)  :: H
        REAL, INTENT(OUT) :: ADVTND(KB)
        INTEGER           :: K
        REAL              :: GRDFTRM, GRDWTRM

        GRDFTRM = 0.5*(H*DZZ(1)) * W(2) * (F(2)-F(1))
        GRDWTRM = F(1) * 0.5*W(2)/(H*DZ(1))
        ADVTND(1) = GRDFTRM + GRDWTRM
        DO K=2,KB-1
           GRDFTRM = 0.5/H*( W(K)*(F(K)-F(K-1))/DZZ(K-1) + W(K+1)*(F(K+1)-F(K))/DZZ(K) )
           GRDWTRM = F(K) * (W(K+1) - W(K)) / (H*DZ(K))
           ADVTND(K) = GRDFTRM + GRDWTRM
        END DO
        ADVTND(KB) = 0

END SUBROUTINE
