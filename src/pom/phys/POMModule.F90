
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
! INTERFACE:
!
   MODULE POM
!
!  DESCRIPTION
!
! This module contains definition and allocation of parameter, scalars and arrays
! used (mostly) by the physical component od the BFM-POM1D system.
!
!***************************************************************************
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem,ONLY: RLEN,ONE
      use constants, ONLY: SEC_PER_DAY
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----SET INTEGER PRECISION------
!
      integer,parameter                 :: ilong=selected_int_kind(12)
!
!     -----VERTICAL LAYERS-----
!
      integer(ilong),parameter :: KB=31
!
!     -----DEFINE SURFACE AND BOTTOM LAYERS LOG DISTRIBUTION-----
!
!     *****************************************************
!     *****************************************************
!     **                                                 **
!     **  INTEGERS DEFINING THE NUMBER OF LOGARITHMIC    **
!     **  LAYERS AT THE SURFACE AND BOTTOM (USED IN      **
!     **  SUBROUTINE CALCDEPTH.                          **
!     **  THE NUMBER OF SURFACE LOG. LAYERS IS KL1-2     **
!     **  THE NUMBER OF BOTTOM  LOG. LAYERS IS KB-KL2-1  **
!     **  FOR NO LOG. DISTRIBUTION EVERYWHERE SET:       **
!     **  KL1=2 AND KL2=KB-1                             **
!     **                                                 **
!     *****************************************************
!     *****************************************************
!
      integer(ilong)          :: KL1, Kl2
!
!     -----SWITCH FOR COLD/HOT START-----
!
!     *****************************************************
!     *****************************************************
!     **                                                 **
!     ** IHOTST=0 "COLD" START FROM INITIAL CONDITION    **
!     ** IHOTST=1 "HOT" START FROM RESTART FILE          **
!     **                                                 **
!     *****************************************************
!     *****************************************************
!
      integer(ilong)           ::IHOTST
!
!     -----MODEL TIME STEP (SECONDS)-----
 ! 
      REAL(RLEN)               :: DTI
!
!     -----LENGTH OF THE RUN (DAYS)-----
!
      integer(ilong)           :: IDAYS
!
!     ----ITERATIONS NEEDED FOR AN "IDAYS" RUN-----
!  
      integer(ilong)           :: IEND
!
!     -----COUNTER FOR THE TIME MARCHING LOOP-----
!
      integer(ilong)           :: intt
!     
!     -----RUNNING TIME (IN DAYS)-----
!
      REAL(RLEN)               :: TIME
!
!     -----TIME AT RESTART-----
!
      REAL(RLEN)               :: TIME0
!
!     -----BOTTOM DEPTH (METRES)-----
!
      REAL(RLEN)               :: H
!   
!     -----LATITUDE & LONGITUDE-----
!  
      REAL(RLEN)               :: ALAT, ALON
!
!     -----GRAVITY-----
!
      real(RLEN),parameter     :: GRAV=9.806
!    
!     -----CORIOLIS PARAMETER-----
!
      REAL(RLEN)               :: COR
!
!     -----BACKGROUND DIFFUSION FOR U, V, Q2, Q2L-----
!
      REAL(RLEN)               :: UMOL
!
!     -----BACKGROUND DIFFUSION FOR T,S & BFM TRACERS
!
      REAL(RLEN)               :: UMOLT,UMOLS,UMOLBFM
!
!     -----NUTRIENT RELAXATION TIME-----
!    
      REAL(RLEN)               :: NRT
!
!     -----PARAMETER FOR THE ASSELIN FILTER----
!            -----(TIME STEP MIXING)-----
!
      REAL(RLEN)               :: SMOTH
!
!     -----REFERENCE DENSITY OF FRESH WATER-----
!
      real(RLEN),parameter         :: RHO0=1000.0
!
!     -----REFERENCE DENSITY OF SEAWATER-----
!
      real(RLEN),parameter         :: RHOSEA=RHO0 + 25.0
!
!     -----SPECIFIC HEAT TIMES RHO0-----
!
      real(RLEN),parameter         :: RCP=4.186*RHOSEA
!
!     -----1 DAY IN SECONDS (RECIPROCAL)-----
!
      real(RLEN), Parameter        :: DAYI=ONE/SEC_PER_DAY
!
!     -----FLAGS TO CHOOSE T,S AND BFM TRACERS SURFACE B.C. IN PROFTS-----
!
      integer(ilong)           :: NBCT, NBCS, NBCBFM
!
!     -----FLAG TO CHOOSE JERLOV WATER TYPE IN PROFTS-----
!
      integer(ilong)           :: NTP
!
!     -----T&S RELAXATION TIME (DAYS) FOR LATERAL ADVECTION-----
!
      integer(ilong)           :: TRT, SRT
!
!     -----DEPTH (m) AT WHICH LATERAL ADVECTION STARTS-----
!
      REAL(RLEN)               :: upperH 
!
!     -----RELAXATION TIME (DAYS) FOR SURFACE SALINITY FLUX-----
!
      REAL(RLEN)               :: SSRT
!
!     -----VERTICAL COORDINATE SYSTEM-----
!     
      real(RLEN),dimension(KB) :: Z,ZZ,DZ,DZZ,DZR     
!
!     -----TEMPERATURE----
!
      real(RLEN),dimension(KB) :: TF,T,TB
!
!     -----SALINITY-----
!
      real(RLEN),dimension(KB) :: SF,S,SB
!
!     -----DENSITY-----
!
      real(RLEN),dimension(KB) :: RHO
!
!     -----VELOCITY-----
!
      real(RLEN),dimension(KB) :: UF,U,UB,VF,V,VB
!
!     -----TURBULENT KINETIC ENERGY (T.K.E.X2)-----
!
      real(RLEN),dimension(KB) :: Q2F,Q2,Q2B
!
!     -----LENGTH SCALE-----
!
      real(RLEN),dimension(KB) :: L
!
!     -----(T.K.E.X2) TIMES LENGTH SCALE-----
!
      real(RLEN),dimension(KB) :: Q2LF,Q2L,Q2LB
!
!     -----VERTICAL DIFFUSION COEFFICIENTS-----
!
      real(RLEN),dimension(KB) ::KM,KH,KQ
!
!     -----SERVICE ARRAYS USED IN POM ROUTINES-----
!
      real(RLEN),public,dimension(KB) :: GM,GH,SM,SH,KN,SPROD,BPROD, A, C, VH, VHP, & 
                                         PROD,DTEF,D, DT
!
!     -----WIND STRESS----
!
      real(RLEN)               :: WUSURF,WVSURF
!
!     -----BOTTOM STRESS----
!
      real(RLEN)               :: WUBOT,WVBOT
!
!     -----SHORT WAVE  RADIATION-----
!
      real(RLEN)               :: SWRAD
!
!     -----BOTTOM ROUGHNESS LENGTH-----
!
real(RLEN)                :: Z0B
!
!     -----MINIMUM ACHIEVABLE BOTTOM DRAG COEFF.-----
!
real(RLEN)                :: CBCMIN
!
!     -----VON KARMANN CONSTANT-----
!
real(RLEN), Parameter    :: vonkarmann=0.40_RLEN
!
!     -----BOTTOM DRAG COEFFICIENT-----
!
real(RLEN)               :: CBC
!
!     *****************************************
!     *****************************************
!     **                                     **
!     ** N.B.                                **
!     ** THE FOLLOWING SCALAR IS USED ONLY   **
!     ** WHEN THE MODEL IS RUN IN PROGNOSTIC **
!     ** MODE                                **
!     **                                     **
!     *****************************************
!     *****************************************
! 
!     ------PRESCRIBED T & S VERTICAL PROFILES-----
!
       real (RLEN),dimension (KB) :: TSTAR, SSTAR

!     ------PRESCRIBED SURFACE T & S-----
!
!     *****************************************
!     *****************************************
!     **                                     **
!     ** N.B.                                **
!     ** TO BE USED (IF DESIRED) AS SURFACE  **
!     ** T & S SURFACE BOUNDARY CONDITION    **
!     ** (IDIAGN=1)                          **
!     ** INTERPOLATED, CLIMATOLOGICAL        **
!     ** AND TIME VARYING.                   **
!     **                                     **
!     *****************************************
!     *****************************************
!
      real(RLEN)               :: TSURF, SSURF
!
!     -----SURFACE SALINITY FLUX-----
!
      real(RLEN)               :: WSSURF
!
!     -----LATERAL ADVECTION FLUX FOR T & S-----
!
      real(RLEN), dimension(KB):: WTADV, WSADV
!
!
!
      CHARACTER(LEN=20)        :: NC_OUT_STARTTIME = '01-01-0000'
!
!     ----- old flag to set prognostic/diagnostic runs. Not used anymore, only diagnostic runs are supported.
      integer(ilong)           ::IDIAGN
      end module POM

! EOC
