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
! !ROUTINE: DENS
!
!DESCRIPTION
!
! !INTERFACE
!
      Subroutine DENS
!
!DESCRIPTION
!
!  This subroutine computes density.
!
!  The algorithm computing density is the one proposed by:
!  Mellor G.L. (1991)
!  An equation of state for numerical models of ocean and estuaries.
!  Journal of atmospheric and oceanic technology, 8, 609-611.
!
!  The input fields are:
!  Potential temperature
!  Salinity
!  Layers depth (ZZ*H)
!
!  Pressure is computed from the hydrostatic eq. assuming the
!  ocean average density.
!
!  The output is density (the reference density is subtracted and then
!  The resulting value is then  divided by 1000).
!  
!  Both input and output are handled via module POM
!
!  **************************************************************************
!
! -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!

      use global_mem,ONLY: RLEN, ZERO, ONE
!
      use POM,ONLY: GRAV, RHO0, H, T, S, ZZ, KB, RHO, RHOSEA
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!    -----LOOP COUNTER-----
!
     INTEGER :: K
!
!    -----LOCAL SCALARS-----
!
     REAL(RLEN) :: CR, P, RHOR, SR, TR, TEMPO
!
!    -----INTRINSIC FUNCTION-----
!
     INTRINSIC ABS
!
!    -----ZEROING-----
!    
    RHOR=ZERO
!
!    -----START COMPUTATION-----

     DO K = 1,KB - 1
!
!        -----TEMPORARY STORAGE FOR T AND S-----
!
!
         TR = T(K)
         SR = S(K)
!
!        -----APPROXIMATE PRESSURE IN UNITS OF BARS-----
!
          P = -GRAV              *              &
               RHOSEA            *              &
               1.0e-5_RLEN       *              &
               ZZ(K)             *              &
               H 
!
!         -----GO FOR DENSITY----
!
          RHOR =  999.842594_RLEN +         &
                 6.793952E-2_RLEN * TR    - &
                 9.095290E-3_RLEN * TR**2 + &
                 1.001685E-4_RLEN * TR**3 - &
                 1.120083E-6_RLEN * TR**4 + &
                 6.536332E-9_RLEN * TR**5
!
          RHOR = RHOR                     + &
                (                           &
                 0.824493_RLEN            - &
                 4.0899E-3_RLEN*TR        + &
                 7.6438E-5_RLEN*TR**2     - &
                 8.2467E-7_RLEN*TR**3     + &
                 5.3875E-9_RLEN*TR**4       &
                )              *SR        + &
                (                           &
                 -5.72466E-3_RLEN         + &
                 1.0227E-4     *TR        - &
                 1.6546E-6_RLEN*TR**2       &
                )                         * &
                 ABS(SR)**1.5             + &
                 4.8314E-4_RLEN*SR**2
!

          CR = 1449.1_RLEN                + &
               0.0821_RLEN      *P        + &
               4.55_RLEN        *TR       - &
               0.045_RLEN       *TR**2    + &
               1.34_RLEN        * (SR-35.)
!
!         **************************************************
!         **************************************************
!         **                                              **
!         ** THE FOLLOWING STATEMENTS SHOULD BE           **
!         ** UNCOMMENTED FOR "DEEP" OCEAN IMPLEMENTATIONS **
!         ** OF THE MODEL (SEE ALSO SUBROUTINE PROFQ).    **
!         **                                              **
!         **************************************************
!         **************************************************
!
!         TEMPO=P/(CR*CR)
!
!         RHOR= RHOR        + &
!                1.0e5_RLEN * &
!                TEMPO      * &
!                (            &
!                ONE        - &
!                 (ONE+ONE) * &
!                 TEMPO       &
!                 )
!
!         -----DENSITY READY FOR USE IN POM-----
!
          RHO(K) = ( RHOR - RHO0 ) / RHO0
 
      END DO
!    
!     -----COSMETIC: T, S, RHO VALUES AT ZZ(KB) ARE NOT USED-----
!
      RHO(KB) = RHO(KB-1)

      RETURN

      end subroutine DENS

