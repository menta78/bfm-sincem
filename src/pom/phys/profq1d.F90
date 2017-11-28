!
! !ROUTINE: Profq
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
! !INTERFACE
!
      SUBROUTINE PROFQ(DT2)
!
!DESCRIPTION
!
! This this is the
! Mellor G.L. and Yamada T (1982)
! Development of a turbulence closure model for geophysical fluid problems.
! Review of Geophysics and Space Physics, 20, 851-875.
! Turbulence closure model.
! It solves for:.
! Turbulent kinetic energy (Q2/2)
! Turbulent length scale (Q2l)
! and provides the turbulent diffusion coefficients for tracers and momentum.
!
! Other useful references are:
! Galperin, B., L. H. Kantha, S. Hassid, and A. Rosati (1988) A quasi-equilibrium
! turbulent energy model for geophysical flows. 
! Journal of Atmospheric Science, 45, 55-62.
!
! Mellor, G. L. (1989)
! Retrospect on oceanic boundary layer modeling and second moment closure.
! Hawaiian Winter Workshop on "Parameterization of Small-Scale Processes",
! University of Hawaii, Honolulu, Hawaii, 1989.
!
!*****************************************************************************
!
!  -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem,ONLY: RLEN, ONE, ZERO
!
      use POM,ONLY: h,      &
                    a,      &
                    c,      &
                    kb,     &
                    kq,     &
                    dzz,    &
                    dz,     &
                    vh,     &
                    vhp,    &
                    wusurf, &
                    wvsurf, &
                    wubot,  &
                    wvbot,  &
                    q2f,    &
                    s,      &
                    t,      &
                    q2lb,   &
                    rho,    &
                    dtef,   &
                    sprod,  &
                    km,     &
                    u,      &
                    v,      &
                    bprod,  &
                    prod,   &
                    q2lf,   &
                    z,      &
                    l,      &
                    sh,     &
                    sm,     &
                    kn,     &
                    kh,     &
                    gm,     &
                    gh,     &
                    zz,     &
                    q2b,    &
                    q2,     &
                    UMOL,   &
                    GRAV,   &
                    H,      &
                    RHOSEA

!
!
!  -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
! -----SCALAR INPUT-----
!
      REAL(RLEN) :: DT2
!
! -----LOCAL SCALARS-----
!
      REAL(RLEN) :: A1,     &
                    A2,     &
                    B1,     &
                    B2,     &
                    C1,     &
                    CIWC,   &
                    COEF1,  &
                    COEF2,  &
                    COEF3,  &
                    COEF4,  &
                    COEF5,  &
                    CONST1, &
                    E1,     &
                    E2,     &
                    KAPPA,  &
                    P,      &
                    SMALL,  &
                    zco,    &
                    coef6
!
!     -----LOOP COUNTERS-----
!
      INTEGER :: K,KI
!
!     -----LOCAL ARRAYS-----
!
      REAL(RLEN) :: BOYGR(KB), &
                    CC(KB)

!
!     -----INTRINSIC FUNCTIONS-----
!
      INTRINSIC ABS,MIN,SQRT
!
!     -----TURBULENCE CLOSURE PARAMETERS-----
!
      DATA A1/0.92_RLEN/
      DATA B1/16.6_RLEN/
      DATA A2/0.74_RLEN/
      DATA B2/10.1_RLEN/
      DATA C1/0.08/
      DATA E1/1.8_RLEN/
      DATA E2/1.33_RLEN/
!
!     -----VON KARMAN CONSTANT-----
!
      DATA KAPPA/0.40_RLEN/
!
!     -----THIS IS TO "MODULATE" WIND AND BOTTOM STRESS (0<CIWC<1)-----
!
      DATA CIWC/ONE/
!
!     -----A SMALL NUMBER...-----
!
      DATA SMALL/1.E-8_RLEN/

      DO K = 2,KB - 1
!
          A(K) = -DT2* (KQ(K+1)+KQ(K)+2.*UMOL)*.5/ (DZZ(K-1)*DZ(K)*H*H)
          C(K) = -DT2* (KQ(K-1)+KQ(K)+2.*UMOL)*.5/ (DZZ(K-1)*DZ(K-1)*H*H)
!
      END DO
!
!    *******************************************************
!    *******************************************************
!    **                                                   **
!    ** THE FOLLOWING SECTION SOLVES THE EQUATION:        **
!    ** DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B        **
!    **                                                   **
!    *******************************************************
!    *******************************************************
!
      CONST1 = B1**(2.0_RLEN/3.0_RLEN)*CIWC
!
!    -----BOUNDARY CONDITION-----
!
      VH(1)   = ZERO
      VHP(1)  = SQRT(WUSURF**2+WVSURF**2)*CONST1
      Q2F(KB) = SQRT(WUBOT**2 +WVBOT**2) *CONST1
!
      DO K = 1,KB - 1
!
!        -----CALCULATE PRESSURE IN UNITS OF DECIBARS-----
!
         P = -GRAV*(RHOSEA)*ZZ(k)*H*1.E-4_RLEN
!
!        -----CALCULATE SPEED OF SOUND SQUARED-----
!
         cc(k) = 1449.1_RLEN +                    &
                 0.00821_RLEN*p   +               &
                 4.55_RLEN*t(k)   -               &
                 0.045_RLEN*t(k)**2 +             &
                 1.34*(S(k)-35.)
!
         cc(k) = cc(k) /                          &
                 sqrt(                            &
                      (ONE-0.01642_RLEN*p/cc(k))* &
                      (ONE-0.40_RLEN*p/cc(k)**2)  &
                     )
      END DO
!
!     -----CALCULATE BUOYANCY GRADIENT-----
!
      DO K = 2,KB - 1
!
          Q2B(K)  =  ABS(Q2B(K))
          Q2LB(K) = ABS(Q2LB(K))
!
          BOYGR(K) = GRAV* &
                     (RHO(K-1)-RHO(K))/ (DZZ(K-1)*H)
!
!         ***********************************************************
!         ***********************************************************
!         **                                                       **
!         ** THIS LINE SHOULD BE COMMENTED IF DENSITY CALCULATION  **
!         ** DOES NOT INCLUDE THE EFFECT OF PRESSURE               **
!         **                                                       **
!         ***********************************************************
!         ***********************************************************
!
!                     +GRAV**2*2.0_RLEN*(RHOSEA)*1.0e-3RLEN/ (CC(K-1)**2+CC(K)**2)
!
          DTEF(K) = Q2B(K)*SQRT(Q2B(K))/ (B1*Q2LB(K)+SMALL)
!
!         -----T.K.E. PRODUCTION-----
!
          SPROD(K) = KM(K)            * &
                     (                  &
                     (U(K)-U(K-1))**2 + &
                     (V(K)-V(K-1))**2   &
                     )                / &
                     (DZZ(K-1)*H)**2  * &
                     CIWC**2
!
          BPROD(K) = KH(K)*BOYGR(K)
!
          PROD(K)  = SPROD(K) + BPROD(K)
!
      END DO
!
!     -----SWEEP DOWNWARD-----
!

      DO K = 2,KB - 1
!
          VHP(K) = ONE/ (A(K)+C(K)* (1.-VH(K-1))- ((ONE+ONE)*DT2*DTEF(K)+ONE))
!
          VH(K) = A(K)*VHP(K)
!
          VHP(K) = (-(ONE+ONE)*DT2*PROD(K)+C(K)*VHP(K-1)-Q2B(K))*VHP(K)
!
      END DO
!
!    -----SWEEP UPWARD-----
!
      DO K = 1,KB - 1
!
          KI = KB - K
!
          Q2F(KI) = VH(KI)*Q2F(KI+1) + VHP(KI)
!
      END DO
!
!    *******************************************************
!    *******************************************************
!    **                                                   **
!    ** THE FOLLOWING SECTION SOLVES THE EQUATION:        **
!    ** DT2*(KQ*Q2L')' - Q2L*(2.*DT2*DTEF+1.) = -Q2BLB    **
!    **                                                   **
!    *******************************************************
!    *******************************************************
!
!    -----BOUNDARY CONDITION-----
!
      VH(1)    = ZERO
      VHP(1)   = ZERO
      Q2LF(KB) = ZERO
!
!    -----SWEEP DOWNWARD-----
!
      DO K = 2,KB - 1
!
          DTEF(K) = DTEF(K)*                   &
                    (                          &
                      ONE+E2*                  &
                       (                       &
                        (                      &
                         ONE/ABS(Z(K)-Z(1))+   &
                         ONE/ABS(Z(K)- Z(KB))  &
                        )*                     &
                        L(K)/(H*KAPPA)         &
                       )**2                    &
                     )
!
          VHP(K) = ONE/ (A(K)+C(K)* (ONE-VH(K-1))- (DT2*DTEF(K)+ONE))
!
          VH(K) = A(K)*VHP(K)
!
          VHP(K) = (DT2* (- (SPROD(K)+BPROD(K))*L(K)*E1)+ &
                  C(K)*VHP(K-1)-Q2LB(K))*VHP(K)
!
      END DO
!
!     -----SWEEP UPWARD-----
!
      DO K = 1,KB - 1
!
          KI = KB - K
!
          Q2LF(KI) = VH(KI)*Q2LF(KI+1) + VHP(KI)
!
      END DO

      DO K = 2,KB - 1
!
          IF (Q2F(K).LE.SMALL .OR. Q2LF(K).LE.SMALL) THEN
              Q2F(K)  = SMALL
              Q2LF(K) = SMALL
          ENDIF
!
      END DO
!
!    *******************************************************
!    *******************************************************
!    **                                                   **
!    ** THE FOLLOWING SECTION SOLVES FOR KM AND KH        **
!    **                                                   **
!    *******************************************************
!    *******************************************************
!
      COEF1 = A2*(ONE-6.0_RLEN*A1/B1)
      COEF2 = 3.0_RLEN*A2*B2 + 18.0_RLEN*A1*A2
      COEF3 = A1*(ONE-3.0_RLEN*C1-6.0_RLEN*A1/B1)
      COEF4 = 18.0_RLEN*A1*A1 + 9.0_RLEN*A1*A2
      COEF5 = 9.0_RLEN*A1*A2
!
      L(1)   = ZERO
      L(KB)  = ZERO
      GH(1)  = ZERO
      GH(KB) = ZERO
!
      DO K = 2,KB - 1
!
          L(K) = Q2LF(K)/Q2F(K)
!
          GH(K) = L(K)**2/Q2F(K)*BOYGR(K)
!
      END DO
!
      DO K = 1,KB
!
!         ---SM AND SH LIMIT TO INFINITY WHEN GH APPROACHES 0.0288---
!
          GH(K) = MIN(GH(K),.028)
!
          SH(K) = COEF1/ (1.-COEF2*GH(K))
!
          SM(K) = COEF3 + SH(K)*COEF4*GH(K)
          SM(K) = SM(K)/ (1.-COEF5*GH(K))
!
      END DO
!
!     -----FINAL COMPUTATION FOR KQ, KN, KM, KH-----
!
      DO K = 1,KB
!
          KN(K) = L(K)*SQRT(ABS(Q2(K)))
!
          KQ(K) = (KN(K)*0.41_RLEN*SH(K)+KQ(K))*0.5_RLEN
!
          KM(K) = (KN(K)*SM(K)+KM(K))*0.5_RLEN
!
          KH(K) = (KN(K)*SH(K)+KH(K))*0.5_RLEN
!
      END DO
!
      RETURN
!
      end subroutine PROFQ

