#include"cppdefs.h"
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
! DESCRIPTION
!
! This is the "main" program of the coupled numerical model originating by
! the direct on-line coupling of the 1D Version of the Princeton Ocean model
! "POM" and the Biological Flux Model "BFM".
!
!****************************************************************************
!
      PROGRAM MAIN
!
! -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem,ONLY: RLEN,ZERO,PI,ONE
!
      use constants, ONLY: SEC_PER_DAY
!
      use POM, ONLY: IDIAGN,                                              &
                     KL1, KL2,                                            &
                     IHOTST,                                              &
                     DTI,                                                 &
                     IDAYS,                                               &
                     IEND,                                                &
                     INTT,                                                &
                     TIME,                                                &
                     TIME0,                                               &
                     H,                                                   &
                     ALAT,                                                &
                     COR,                                                 &
                     UMOL,                                                &
                     UMOLT,UMOLS,UMOLBFM,                                 &
                     SMOTH,                                               &
                     DAYI,                                                &
                     KB,                                                  &
                     NBCT,NBCS,NBCBFM,                                    &
                     NTP,                                                 &
                     TRT, SRT,                                            &
                     upperH,                                              &
                     SSRT,                                                &
                     Z,ZZ,DZ,DZZ,DZR,                                     &
                     TF,T,TB,                                             &
                     SF,S,SB,                                             &
                     RHO,                                                 &
                     UF,U,UB,VF,V,VB,                                     &
                     Q2F,Q2,Q2B,                                          &
                     L,                                                   &
                     Q2LF,Q2L,Q2LB,                                       &
                     KM,KH,KQ,                                            &
                     GM,GH,SM,SH,KN,SPROD,BPROD, A, C, VH, VHP,PROD,DTEF, &
                     WUSURF,WVSURF,                                       &
                     WUBOT,WVBOT,                                         &
                     SWRAD,                                               &
                     TSURF, SSURF,                                        &
                     TSTAR, SSTAR,                                        &
                     WSSURF,                                              &
                     WTADV, WSADV,                                        &
                     NRT,                                                 &
                     Z0B, CBCMIN, CBC,                                    &
                     D, DT,                                               &
                     vonkarmann,                                          &
                     NC_OUT_STARTTIME
!
      use CPL_VARIABLES,ONLY: ilong, savef, USE_KH_EXT
!
      use Forcing,ONLY: Forcing_manager
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----NAMELIST READING UNIT------
!
      integer(ilong),parameter            :: namlst=10
!
!     -----LOOP COUNTER-----
!
      integer(ilong)               :: K
!
!     -----TWICE THE TIME STEP-----
!
      REAL(RLEN)                   :: DT2
!
!     -----EARTH ROTATION ANGULAR VELOCITY-----
!
      REAL(RLEN),parameter         ::OMEGA=7.29E-5
!
!     -----EXTERNAL SUBROUTINES-----
!
      EXTERNAL DENS,            &
               PROFQ,           &
               CALCDEPTH,       &
               get_init_TS_IC,       &
               load_restart,    &
               pom_ini_bfm_1d,  &
               pom_bfm_1d,      &
               save_restart
!
!
!     -----INTRINSIC FUNCTIONS-----
!
      INTRINSIC FLOAT,MOD,NINT,SIN,IFIX,MAX,LOG
!
!     -----NAMELIST READING-----
!
      NAMELIST /Params_POMBFM/ H,DTI,ALAT,IDIAGN,IDAYS,SMOTH,&
                               ihotst,UMOL,KL1,KL2,savef,NRT,NBCT,NBCS,NBCBFM,&
                               UMOL,UMOLT,UMOLS,UMOLBFM,NTP,TRT,SRT,UPPERH,SSRT, &
                               CBCMIN,Z0B,NC_OUT_STARTTIME
!
      OPEN(namlst,file='pom_bfm_settings.nml',status='old',action='read')
      READ(namlst,nml=Params_POMBFM)
      CLOSE(namlst)
!
!     -----GENERAL INITIALISATION (ZEROING)-----
!
          Z(:)           = ZERO
          ZZ(:)          = ZERO
          DZ(:)          = ZERO
          DZZ(:)         = ZERO
          TF(:)          = ZERO
          T(:)           = ZERO
          TB(:)          = ZERO
          SF(:)          = ZERO
          S(:)           = ZERO
          SB(:)          = ZERO
          RHO(:)         = ZERO
          UF(:)          = ZERO
          U(:)           = ZERO
          UB(:)          = ZERO
          VF(:)          = ZERO
          V(:)           = ZERO
          VB(:)          = ZERO 
          Q2F(:)         = 1.E-7_RLEN
          Q2(:)          = Q2F(:)
          Q2B(:)         = Q2(:)
          Q2LF(:)        = 1.E-7_RLEN
          Q2L(:)         = Q2LF(:)
          Q2LB(:)        = Q2L(:)
          L(:)           = 1.0_RLEN
          L(1)           = ZERO
          L(KB)          = ZERO
          KM(:)          = ZERO
          KH(:)          = ZERO
          KQ(:)          = ZERO 
          GM(:)          = ZERO
          GH(:)          = ZERO
          SM(:)          = ZERO
          SH(:)          = ZERO
          KN(:)          = ZERO
          SPROD(:)       = ZERO
          BPROD(:)       = ZERO
          PROD(:)        = ZERO    
          VH(:)          = ZERO
          VHP(:)         = ZERO
          DTEF(:)        = ZERO
          D(:)           = ZERO
          DT(:)          = ZERO
          A(:)           = ZERO 
          C(:)           = ZERO
          TSTAR(:)       = ZERO
          SSTAR(:)       = ZERO
          WTADV(:)       = ZERO
          WSADV(:)       = ZERO
          WUSURF         = ZERO
          WVSURF         = ZERO
          WUBOT          = ZERO
          WVBOT          = ZERO
          SWRAD          = ZERO 
          WSSURF         = ZERO
          TSURF          = ZERO
          SSURF          = ZERO
          CBC            = ZERO
!
!     -----DEFINE VERTICAL COORDINATE SYSTEM-----
!
      call calcdepth
!     
      dz(kb)=1.e-6_RLEN
      dzr(:)=ONE/dz(:)
      dz(kb)=ZERO
!
!     -----CORIOLIS PARAMETER-----
!
      COR = 2.0_RLEN*OMEGA*SIN(ALAT*2.0_RLEN*PI/360.0_RLEN)
!
!     -----TWICE THE TIME STEP-----
!
      DT2  = 2.0_RLEN*DTI
!
!     -----COMPUTE BOTTOM DRAG COEFFICIENT-----
!
      CBC=LOG((ZZ(KB-1)-Z(KB))*H/Z0B)
      CBC=MAX(CBCMIN,vonkarmann**2/CBC**2)
!
!     -----ITERATIONS NEEDED TO CARRY OUT AN"IDAYS" SIMULATION-----
!
      iend = idays*IFIX(SEC_PER_DAY)/IFIX(dti)

!     opening the input files
      call opendat
!
!     -----READ  T&S INITIAL CONDITIONS (IHOTST=0) OR RESTART FILE (IHOTST=1)-----
!

      select case (ihotst)
!
             case (0)
!
                  time0=ZERO
!
                  call get_init_TS_IC
!       
!                 -----DEFINE INITIAL DENSITY FIELD----
!    
                  call DENS
!
             case (1)
!
!                 -----READ RESTART-----
!
                  call load_restart
! 
       end select 
!
#ifndef POM_only
!
! -----BFM SET-UP, INITIALISATION AND RESTART READING (IF REQUIRED)-----
!
      call pom_ini_bfm_1d
!
#endif
!
!     ***************************
!     ***************************
!     **                       **
!     ** BEGIN THE TIME MARCH  **
!     **                       **
!     ***************************
!     ***************************
!
      DO intt = 1, IEND                        
!
!        -----COMPUTE TIME IN DAYS-----
!
         time = time0+(DTI*float(intt)*dayi)
!
!        -----DEFINE ALL FORCINGS-----
!
         call forcing_manager
!
         IF (.NOT. USE_KH_EXT) THEN
!            ----- MELLOR-YAMADA (1982) TURBULENCE CLOSURE-----
             Q2F(:)  = Q2B(:)
             Q2LF(:) = Q2LB(:)
!
             call PROFQ(DT2)
!
!            -----COMPUTE VELOCITY-----
!
             UF(:) = UB(:) + DT2 * COR * V(:) ! CORIOLIS TERM & EXPLICIT LEAPFROG
!
             call PROFUV(DT2,UF,WUSURF,WUBOT)
!
             VF(:) = VB(:) - DT2 * COR * U(:) ! CORIOLIS TERM & EXPLICIT LEAPFROG
!
             call PROFUV(DT2,VF,WVSURF,WVBOT)
!
!            -----MIXING TIME STEP (ASSELIN FILTER)-----
!
             Q2(:)   = Q2(:)  + 0.5_RLEN * SMOTH * (Q2F(:)  + Q2B(:)  - 2.0_RLEN * Q2(:))
             Q2L(:)  = Q2L(:) + 0.5_RLEN * SMOTH * (Q2LF(:) + Q2LB(:) - 2.0_RLEN * Q2L(:))
!
             U(:)    = U(:)   + 0.5_RLEN * SMOTH * (UF(:)   + UB(:)   - 2.0_RLEN * U(:))
             V(:)    = V(:)   + 0.5_RLEN * SMOTH * (VF(:)   + VB(:)   - 2.0_RLEN * V(:))
!
!            -----RESTORE TIME SEQUENCE-----
!
             Q2B(:)  = Q2(:)
             Q2(:)   = Q2F(:)
             Q2LB(:) = Q2L(:)
             Q2L(:)  = Q2LF(:)
!
             UB(:) = U(:)
             U(:)  = UF(:)
             VB(:) = V(:)
             V(:)  = VF(:)
         END IF
!
         TB(:) = T(:)
         T(:)  = TF(:)
         SB(:) = S(:)
         S(:)  = SF(:)
!
!        -----UPDATE DENSITY-----
!
         call DENS
!
!        **********************************************
!        **********************************************
!        **                                          **
!        ** PHYSICS HAS BEEN COMPUTED.....GO BFM!!!  **
!        **                                          **
!        **********************************************
!        **********************************************
!
#ifndef POM_only
!
         call pom_bfm_1d
!
#endif
!
      END DO
!
!     *******************************************
!     *******************************************
!     **                                       **
!     **           END TIME MARCH              **
!     **                                       **
!     *******************************************
!     *******************************************
!
!     ----- WRITING POM RESTART-----
!
      WRITE (71) time,                      &
                 u,ub,v,vb,                 &
                 t,tb,s,sb,                 &
                 q2,q2b,q2l,q2lb,           &
                 kh,km,kq,                  &
                 l,                         &
                 wubot,wvbot,               &
                 RHO
!
!     ----- WRITING BFM RESTART-----
!
#ifndef POM_only

!
      call save_restart
!
#endif
!
       PRINT *,'Main done'
!
     end program MAIN 
