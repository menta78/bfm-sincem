#include "INCLUDE.h"
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
!
! !ROUTINE: SCHEME_BENTHIC_LF1D
!
!
! !INTERFACE
!
  SUBROUTINE SCHEME_BENTHIC_LF1D
!
!DESCRIPTION
!
!  This routine calculates and solves for the leap-frog time step
!  applied to the BFM scalar state variables (the benthic state variables)
!
!***************************************************************************
!
!   -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
    use mem_Param, ONLY: p_small
!
    use global_mem, ONLY:RLEN,ZERO
!
    use Pom, ONLY: smoth,dti, ilong
!
    use Mem, ONLY:D2STATE_BEN,NO_D2_BOX_STATES_BEN,D2SOURCE_BEN,NO_BOXES_BEN
!
    use api_bfm, ONLY:D2STATEB_BEN
!
!    -----IMPLICIT TYPING IS NEVER ALLOWED----
!
     IMPLICIT NONE
!
!    ----- LOOP COUNTER-----
!
     integer(ilong) ::  n
!
!    -----TWICE THE TIME STEP-----
!
     real(RLEN) :: dti2
!
!    -----TEMPORARY STORAGE FOR STATE VARIABLES COMPUTED @ t+dt-----
!
     real(RLEN) :: tempo(NO_D2_BOX_STATES_BEN)
!
!    -----ZEROING-----
!
     tempo = ZERO
!
!    -----TWICE THE TIME STEP-----
!
     dti2 = dti*2.0_RLEN
!
!    -----COMPUTE SOURCE/SINKS TREND-----
!
     do n = 1, NO_D2_BOX_STATES_BEN
        tempo(n) = tempo(n) + D2SOURCE_BEN(1,n)
     end do 
!
!-----LEAP FROG INTEGRATION-----
!
     tempo=D2STATEB_BEN(1,:)+(tempo*dti2)
!         
!    -----CLIPPING (IF NEEDED....)-----
!
    do  n = 1,NO_D2_BOX_STATES_BEN
!
            tempo(n)=max(p_small,tempo(n))
!
    end do
!
!    -----ASSELIN FILTER-----
!
     D2STATE_BEN(1,:)=D2STATE_BEN(1,:)          + &
                      0.5_RLEN*smoth            * &
                      (                           &
                       tempo(:)                 + &
                       D2STATEB_BEN(1,:)        - &
                       2.0_RLEN*D2STATE_BEN(1,:)  &
                      )
!
!    -----RESTORE TIME SEQUENCE-----
!
     D2STATEB_BEN(1,:)=D2STATE_BEN(1,:)
     D2STATE_BEN(1,:)=tempo(:)
!
    return
!
      end subroutine SCHEME_BENTHIC_LF1D




SUBROUTINE SCHEME_BENTHIC_EFW
    ! super simple Euler-Forward scheme
    USE Pom, ONLY: DTI
    USE Mem, ONLY: D2STATE_BEN, D2SOURCE_BEN

    D2STATE_BEN = D2STATE_BEN + DTI*D2SOURCE_BEN

END SUBROUTINE SCHEME_BENTHIC_EFW



SUBROUTINE BENTHIC_TIME_INTEGRATION
   USE CPL_VARIABLES, ONLY: BENT_INT_SCHEME

   SELECT CASE (BENT_INT_SCHEME)
       CASE (1)
           CALL SCHEME_BENTHIC_EFW
       CASE DEFAULT
           CALL SCHEME_BENTHIC_LF1D
   END SELECT

END SUBROUTINE BENTHIC_TIME_INTEGRATION
