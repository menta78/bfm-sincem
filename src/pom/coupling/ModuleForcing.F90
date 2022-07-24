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
  MODULE  Forcing
!
! DESCRIPTION
!
! DEFINITION AND ALLOCATION OF PARAMETERS, SCALAR AND ARRAYS USED---
! TO DEFINE THE PHYSICAL AND  BIOGEOCHEMICAL FORCING
!
!**********************************************************************
!
! -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
  use POM,ONLY: KB, ilong, NUTSBC_MODE, DZ
!
  use global_mem, ONLY:RLEN
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
  IMPLICIT NONE
!
!     -----SHORTWAVE RADIATION DATA-----
! 
!     ***************************************
!     ***************************************
!     ** N.B.!                             **
!     ** ALWAYS NEEDED: WHEN THE MODEL IS  **
!     ** RUN IN DIAGNOSTIC MODE PROVIDES   **
!     ** ONLY PAR TO BFM. IN PROGNOSTIC    **
!     ** CONTRIBUTES TO THE DEFINITION OF  **
!     ** THE TEMPERATURE SURFACE BOUNDARY  **
!     ** CONDITION.                        **
!     ***************************************
!     ***************************************
!
      REAL(RLEN),SAVE                     :: SWRAD1,SWRAD2
!
!     ***************************************
!     ***************************************
!     ** N.B.!                             **
!     ** THE FOLLOWING 2 SCALARS ARE USED  **
!     ** ONLY WHEN THE MODEL IS RUN IN     **
!     ** PROGNOSTIC MODE.                  **
!     ***************************************
!     ***************************************
!    
!     -----LOSS TERM OF THE SURFACE HEAT FLUX DATA-----

       REAL (RLEN),SAVE                     :: WTSURF1,WTSURF2
!
!     -----HEAT FLUX CORRECTION TERM (NO MORE IN USE)------
!
      REAL (RLEN),SAVE                     :: QCORR1,QCORR2
!     ***************************************
!     ***************************************
!     ** N.B.!                             **
!     ** THE FOLLOWING SCALARS ARE ALWAYS  **
!     ** USED.                             **
!     **                                   **
!     ***************************************
!     ***************************************
!
!     -----WIND STRESS DATA-----
!
      REAL (RLEN),SAVE                       :: WSU1,WSU2,WSV1,WSV2
!
!     -----SURFACE NITRATE DATA-----
!
      REAL (RLEN),SAVE                       :: NO3_1,NO3_2
!
!     -----SURFACE PHOSPHATE DATA-----
!
      REAL (RLEN),SAVE                       :: PO4_1,PO4_2
!
!     -----MONTHLY SURFACE AMMONIA DATA-----
!
      REAL (RLEN),SAVE                       :: NH4_1,NH4_2
!
!     -----MONTHLY SURFACE SILICATE-----
!
      REAL (RLEN),SAVE                       :: SIO4_1,SIO4_2
!
!     -----MONTHLY SURFACE DISCHARGE, used only if NUTSBC_MODE == 1 -----
!     -------- units: kg/m2/s, as in nemo
      REAL (RLEN),SAVE                       :: DIS_1,DIS_2
!
!     -----MONTHLY PROFILES OF CURRENTS SPEED, used only if NUTSBC_MODE == 1
      REAL (RLEN), DIMENSION(KB-1)             :: CUR_1,CUR_2
      REAL (RLEN)                              :: CUR_MEAN
!
!     -----VERTICAL PROFILES OF INORGANIC SUSPENDED MATTER DATA-----
!
      real(RLEN),public,dimension(KB-1),SAVE :: ISM1,ISM2
!
!     -----VERTICAL PROFILES OF T & S-----
!
      real(RLEN),public,dimension(KB),SAVE :: TCLIM1,TCLIM2
      real(RLEN),public,dimension(KB),SAVE :: SCLIM1,SCLIM2

!
!     -----INTERPOLATORS AND COUNTERS-----
!
      INTEGER(ilong),SAVE                   ::  ICOUNTF,        &
                                                IFCHGE,         &
                                                IFINT
!
      REAL(RLEN),SAVE                       ::  RATIOF
!
contains
!
!
      SUBROUTINE FORCING_MANAGER
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
!This subroutine reads the physical and biogeochemical forcing data and carries
! out the time linear interpolation.
!
! N.B. The subroutine is currently arranged to work with monthly climatological
! data series. A more general version is going to be prepared.
!
!*******************************************************************************
!
! -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem,ONLY: RLEN,ZERO,ONE,NML_OPEN,NML_READ,error_msg_prn
!
      use constants, ONLY: SEC_PER_DAY
!
      use POM, ONLY: IDIAGN,                                              &
                     DTI,                                                 &
                     INTT,                                                &
                     RCP,                                                 &
                     KB,                                                  &
                     TF,                                                  &
                     SF ,                                                 &
                     WUSURF,WVSURF,                                       &
                     SWRAD,                                               &
                     WTSURF,WSSURF,                                       &
                     TSURF,SSURF,                                         &
                     TSTAR,SSTAR,                                         &
                     ilong,                                               &
                     RHO0
!
      use Service,ONLY: ISM,PO4SURF,NO3SURF,NH4SURF,SIO4SURF,DISSURF,CURRENTS_SPEED_PROF,Cprofile_input_exist
!
! -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----LOOP COUNTER-----
!
      integer(ilong)               :: K
!
!     -------INITIALISATION AND FIRST FORCING READING-----
!
      if(intt==INT(ONE)) then
!
!       -----INITIALISE (ZEROING)-----
!   
        SWRAD1       =ZERO
        SWRAD2       =ZERO
        WTSURF1      =ZERO
        WTSURF2      =ZERO
        WSU1         =ZERO
        WSU2         =ZERO
        WSV1         =ZERO
        WSV2         =ZERO
        NO3_1        =ZERO
        NO3_2        =ZERO
        PO4_1        =ZERO
        PO4_2        =ZERO
        NH4_1        =ZERO
        NH4_2        =ZERO
        SIO4_1       =ZERO
        SIO4_2       =ZERO
        DIS_1        =ZERO
        DIS_2        =ZERO
        TCLIM1(:)    =ZERO
        TCLIM2(:)    =ZERO
        SCLIM1(:)    =ZERO
        SCLIM2(:)    =ZERO
        ISM1(:)      =ZERO
        ISM2(:)      =ZERO
!
!         -----DATA READING COUNTER-----
!
          ICOUNTF =  1
!
!         -----TIME STEPS TO COVER ONE MONTH----
!
          IFCHGE  = 30_ilong*INT(SEC_PER_DAY)/INT(DTI)
!
!         -----MONTH INTERPOLATOR---
!
!         ******************************************************************
!         ******************************************************************
!         **                                                              **
!         ** THE MONTHLY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE    **
!         ** CENTERED AT DAY 15 OF EACH CLIMATOLOGICAL MONTH. THEREFORE   **
!         ** THE MONTH INTERPOLATOR (IFINT) IS INITIALISED AT THE VALUE   **
!         ** (IFCHGE/2)-1 CORRESPONDING TO DAY 15 MINUS 1 TIMESTEP.       **
!         **                                                              **
!         ******************************************************************
!         ******************************************************************
!
          IFINT   = (IFCHGE/2_ilong)-INT(ONE)
!
!         -----OPEN DATA FILES-----
!
          CALL opendat
!
!         *************************************************
!         *************************************************
!         **                                             **
!         **  INITIAL READING OF THE MONTHLY FORCING     **
!         **                                             **
!         *************************************************
!         *************************************************
!
!         -----WIND STRESS-----
!
          READ (11,REC=ICOUNTF)   WSU1,WSV1
          READ (11,REC=ICOUNTF+1) WSU2,WSV2
!
!         -----SURFACE HEAT FLUX-----
!
!         **************************************************
!         **************************************************
!         **                                              **
!         ** N.B.: IF THE MODEL IS RUN IN DIAGNOSTIC MODE **
!         ** ONLY THE SOLAR RADIATION  (SWRAD#) IS USED   **
!         **                                              **
!         **************************************************
!         **************************************************
!
          select case (IDIAGN)
!
                 case(0)
!
                    READ(21,REC=ICOUNTF)   SWRAD1,WTSURF1
                    READ(21,REC=ICOUNTF+1) SWRAD2,WTSURF2
!
                 case(1)
!
                    READ (21,REC=ICOUNTF)   SWRAD1
                    READ (21,REC=ICOUNTF+1) SWRAD2
!
           end select
!
!         ***********************************************************
!         ***********************************************************
!         **                                                       **
!         ** ELIMINATE SOLAR RADIATION from TOTAL HEAT FLUX        **
!         ** N.B.: THE FOLLOWING TWO LINES SHOULD BE ACTIVE IF     **
!         **       (AND ONLY IF!) IN WTSURF# THERE IS THE TOTAL    **
!         **       HEAT FLUX AND NOT ONLY THE LOSS TERMS           **
!         **                                                       **
!         ***********************************************************
!         ***********************************************************
!
!         WTSURF1 = WTSURF1-SWRAD1
!         WTSURF2 = WTSURF2-SWRAD2
!
!
!         -----CLIMATOLOGICAL T&S PROFILES-------
!
!         ***********************************************************
!         ***********************************************************
!         **                                                       **
!         ** N.B.:IF THE MODEL IS RUN IN DIAGNOSTIC MODE  THE      **
!         **      WHOLE PROFILES ARE USED. IN PROGNOSTIC MODE      **
!         **      ONLY THE SURFACE T&S VALUES ARE RETAINED AND     **
!         **      STORED IN TSURF AND SURF AS OPTION FOR THE       **
!         **      T&S SURFACE BOUNDARY CONDITION                   **
!         **                                                       **
!         ***********************************************************
!         ***********************************************************
!
          DO K = 1,KB
             READ (20,REC=(ICOUNTF-1)*KB+K) SCLIM1(K)
             READ (15,REC=(ICOUNTF-1)*KB+K) TCLIM1(K)
             READ (20,REC=ICOUNTF*KB+K)     SCLIM2(K)
             READ (15,REC=ICOUNTF*KB+K)     TCLIM2(K)
          END DO
!
!         -----SUSPENDED INORGANIC MATTER PROFILES-----
!
          DO K = 1,KB-1
             READ (19,REC=(ICOUNTF-1)*(KB-1)+K) ISM1(K)
             READ (19,REC=ICOUNTF*(KB-1)+K)     ISM2(K)
          END DO
!
!         -----SURFACE NUTRIENTS-----
!
          SELECT CASE (NUTSBC_MODE)
             CASE (1)
                  ! reading the surface concentration and the discharge
                  READ(18,REC=ICOUNTF)   NO3_1,NH4_1,PO4_1,SIO4_1,DIS_1
                  READ(18,REC=ICOUNTF+1) NO3_2,NH4_2,PO4_2,SIO4_2,DIS_2
             CASE DEFAULT
                  ! reading the surface concentration
                  READ(18,REC=ICOUNTF)   NO3_1,NH4_1,PO4_1,SIO4_1
                  READ(18,REC=ICOUNTF+1) NO3_2,NH4_2,PO4_2,SIO4_2
          END SELECT
!
!         -----PROFILES OF CURRENTS SPEED----
!
          IF ((NUTSBC_MODE .EQ. 1) .AND. Cprofile_input_exist) THEN
             DO K = 1,KB-1
                READ (35,REC=(ICOUNTF-1)*(KB-1)+K) CUR_1(K)
                READ (35,REC=ICOUNTF*(KB-1)+K)     CUR_2(K)
                IF ((CUR_1(K) .LT. ZERO) .OR. (CUR_2(K) .LT. ZERO)) THEN
                   WRITE(6,*) 'FATAL ERROR!!!! THE FILE WITH CURRENTS SPEED SHOULD ONLY CONTAIN POSITIVE VALUES!!!'
                   STOP 999
                END IF
             END DO
             ! normalizing to 1
             CUR_MEAN = SUM(CUR_1*DZ(:KB-1))/SUM(DZ(:KB-1)) ! mean weighted on DZ
             CUR_1 = CUR_1/CUR_MEAN
             CUR_MEAN = SUM(CUR_2*DZ(:KB-1))/SUM(DZ(:KB-1)) ! mean weighted on DZ
             CUR_2 = CUR_2/CUR_MEAN
          ELSE
             CUR_1 = 1
             CUR_2 = 1
          END IF

#ifdef SAVEFORCING
          write(400,'(1i8, 8e18.8)') ICOUNTF, WSU1,WSV1,SWRAD1,NO3_1,NH4_1,PO4_1,SIO4_1,DIS_1
          write(401,'(1i8,40e18.8)') ICOUNTF, ISM1
          write(402,'(1i8,40e18.8)') ICOUNTF, SCLIM1
          write(403,'(1i8,40e18.8)') ICOUNTF, TCLIM1
          write(400,'(1i8, 8e18.8)') ICOUNTF+1, WSU2,WSV2,SWRAD2,NO3_2,NH4_2,PO4_2,SIO4_2,DIS_2
          write(401,'(1i8,40e18.8)') ICOUNTF+1, ISM2
          write(402,'(1i8,40e18.8)') ICOUNTF+1, SCLIM2
          write(403,'(1i8,40e18.8)') ICOUNTF+1, TCLIM2
#endif
!
!         -----WIND STRESS CONVERTED TO POM UNITS (N/m2-->m2/s2)-----
!
          WSU1 = WSU1*(-ONE)/RHO0
          WSU2 = WSU2*(-ONE)/RHO0
          WSV1 = WSV1*(-ONE)/RHO0
          WSV2 = WSV2*(-ONE)/RHO0
!
!         -----HEAT FLUX CONVERTED TO POM UNITS(W/m2-->deg.C*m/s)-----
!
          SWRAD1  = SWRAD1*(-ONE)/rcp
          SWRAD2  = SWRAD2*(-ONE)/rcp
!   
          IF(IDIAGN==INT(ZERO)) THEN
             WTSURF1 = WTSURF1*(-ONE)/rcp
             WTSURF2 = WTSURF2*(-ONE)/rcp
          ENDIF
!
!         -----UPDATE THE MONTH COUNTER-----
!
          ICOUNTF = ICOUNTF + INT(ONE)
!
      endif
!
!     -----END OF INITIALISATION AND FIRST DATA READING (intt=1)-----
!
!     -----BEGIN TIME INTERPOLATION-----
!
!     -----UPDATE INTERPOLATION COUNTER-----
!
      IFINT  = IFINT + INT(ONE)
!
!     -----UPDATE INTERPOLATOR-----
!
      RATIOF = FLOAT(IFINT)/FLOAT(IFCHGE)
!
!     -----INTERPOLATE WIND STRESS-----
!
      WUSURF = WSU1 + RATIOF * (WSU2-WSU1)
      WVSURF = WSV1 + RATIOF * (WSV2-WSV1)
!
      select case (IDIAGN)
!
             case (INT(ZERO))
!
!                 -----PROGNOSTIC RUN: INTERPOLATE TOTAL HEAT FLUX-----
!
                  WTSURF = WTSURF1 + RATIOF * (WTSURF2-WTSURF1)
                  SWRAD  = SWRAD1  + RATIOF * (SWRAD2-SWRAD1)
!
             case (INT(ONE))
!
!                  -----DIAGNOSTIC RUN: INTERPOLATE SOLAR RADIATION ONLY-----
!
                  SWRAD  = SWRAD1  + RATIOF * (SWRAD2-SWRAD1)
!
      end select
!
!     -----INTERPOLATE T&S PROFILES-----
!
      TSTAR(:) = TCLIM1(:) + RATIOF * (TCLIM2(:)-TCLIM1(:))
      SSTAR(:) = SCLIM1(:) + RATIOF * (SCLIM2(:)-SCLIM1(:))
!
      select case (IDIAGN)
!
             case (INT(ZERO))
!
!                 -----PROGNOSTIC RUN: SAVE SST AND SSS FOR SURF. B.C.-----
!
                  TSURF = TSTAR(1)
                  SSURF = SSTAR(1)
!
             case (INT(ONE))
!
!                 -----DIAGNOSTIC RUN: THE WHOLE T&S PROFILES ARE PRESCRIBED-----
!
                  TF(:) = TSTAR(:)
                  SF(:) = SSTAR(:)
!
      end select
!
!     -----INTERPOLATE SUSPENDED INORGANIC MATTER-----
!
      ISM(:) = ISM1(:) + RATIOF * (ISM2(:)-ISM1(:))
!
!     -----INTERPOLATE SURFACE NUTRIENTS-----

      NO3SURF  = NO3_1  + RATIOF * (NO3_2-NO3_1)
      NH4SURF  = NH4_1  + RATIOF * (NH4_2-NH4_1)
      PO4SURF  = PO4_1  + RATIOF * (PO4_2-PO4_1)
      SIO4SURF = SIO4_1 + RATIOF * (SIO4_2-SIO4_1)
      DISSURF  = DIS_1  + RATIOF * (DIS_2-DIS_1)
      CURRENTS_SPEED_PROF = CUR_1 + RATIOF * (CUR_2-CUR_1)
!
!     -----END OF INTERPOLATION SECTION-----
!
!     -----BEGIN DATA UPDATE SECTION-----
!
      IF (IFINT==IFCHGE) THEN
!
!        -----A MONTH HAS GONE...IT IS NECESSARY TO....-----
!
!        -----....UPDATE MONTH COUNTER....-----
!
         ICOUNTF = ICOUNTF + 1
         PRINT *, 'ICOUNTF', ICOUNTF
!
!        -----....RESET INTERPOLATION COUNTER....-----
!
         IFINT   = INT(ZERO)
!
!        -----....SHIFT THE MONTHLY DATA....-----
!
         SWRAD1     = SWRAD2
         if (IDIAGN==INT(ZERO)) WTSURF1 = WTSURF2
         WSU1       = WSU2
         WSV1       = WSV2
         NO3_1      = NO3_2
         NH4_1      = NH4_2
         PO4_1      = PO4_2
         SIO4_1     = SIO4_2
         DIS_1      = DIS_2
         CUR_1      = CUR_2      
         ISM1(:)    = ISM2(:)
         TCLIM1(:)  = TCLIM2(:)
         SCLIM1(:)  = SCLIM2(:)
!
!            -----IF 12 MONTHS HAVE GONE.... -----
!
         IF (ICOUNTF.GT.13) THEN
!
!            -----RESTART THE READING SEQUENCE-----
!
             ICOUNTF = 2
!
             READ (11,REC=1) WSU1,WSV1
!
!            -----POM UNITS-----
!
             WSU1 = WSU1*(-ONE)/RHO0
             WSV1 = WSV1*(-ONE)/RHO0
!
             select case (IDIAGN)
!                   
                    case(INT(ZERO))
!
                         READ(21,REC=1) SWRAD1,WTSURF1
!          
!                        -----POM UNITS-----
!
                         WTSURF1=WTSURF1*(-ONE)/rcp
!
                    case(INT(ONE))
!
                         READ(21,REC=1) SWRAD1
!
             end select
!
!            -----POM UNITS-----
!
             SWRAD1=SWRAD1*(-ONE)/rcp
!
             SELECT CASE (NUTSBC_MODE)
                CASE (1)
                     ! reading the surface concentration and the discharge
                     READ(18,REC=1) NO3_1,NH4_1,PO4_1,SIO4_1, DIS_1
                CASE DEFAULT
                     ! reading the surface concentration
                     READ(18,REC=1) NO3_1,NH4_1,PO4_1,SIO4_1
             END SELECT
!
             DO K = 1,KB
                READ (20,REC=K) SCLIM1(K)
                READ (15,REC=K) TCLIM1(K)
             END DO
!
             DO K = 1, KB-1
                READ (19,REC=K) ISM1(K)
             END DO

             IF ((NUTSBC_MODE .EQ. 1) .AND. Cprofile_input_exist) THEN
                DO K = 1,KB-1
                   READ (35,REC=K) CUR_1(K)
                   IF (CUR_1(K) .LT. ZERO) THEN
                      WRITE(6,*) 'FATAL ERROR!!!! THE FILE WITH CURRENTS SPEED SHOULD ONLY CONTAIN POSITIVE VALUES!!!'
                      STOP 999
                   END IF
                END DO
                ! normalizing to 1
                CUR_MEAN = SUM(CUR_1*DZ(:KB-1))/SUM(DZ(:KB-1)) ! mean weighted on DZ
                CUR_1 = CUR_1/CUR_MEAN
             ELSE
                CUR_1 = 1
             END IF
!
         END IF
!
!        -----READ FOLLOWING MONTH-----
!
         SELECT CASE (NUTSBC_MODE)
            CASE (1)
                 ! reading the surface concentration and the discharge
                 READ(18,REC=ICOUNTF) NO3_2,NH4_2,PO4_2,SIO4_2, DIS_2
            CASE DEFAULT
                 ! reading the surface concentration
                 READ(18,REC=ICOUNTF) NO3_2,NH4_2,PO4_2,SIO4_2
         END SELECT

         READ (11,REC=ICOUNTF) WSU2,WSV2
!
         select case (IDIAGN)
!
               case (INT(ZERO))
!
                     READ (21,REC=ICOUNTF) SWRAD2,WTSURF2
!
!                    -----POM UNITS-----
!
                     WTSURF2=WTSURF2*(-ONE)/rcp
!
               case (INT(ONE))
!
                     READ (21,REC=ICOUNTF) SWRAD2
!
        end select
!
         IF ((NUTSBC_MODE .EQ. 1) .AND. Cprofile_input_exist) THEN
            DO K = 1,KB-1
               READ (35,REC=(ICOUNTF-1)*(KB-1)+K)     CUR_2(K)
               IF (CUR_2(K) .LT. ZERO) THEN
                  WRITE(6,*) 'FATAL ERROR!!!! THE FILE WITH CURRENTS SPEED SHOULD ONLY CONTAIN POSITIVE VALUES!!!'
                  STOP 999
               END IF
            END DO
            CUR_MEAN = SUM(CUR_2*DZ(:KB-1))/SUM(DZ(:KB-1)) ! mean weighted on DZ
            CUR_2 = CUR_2/CUR_MEAN
         ELSE
            CUR_1 = 1
         END IF
!
#ifdef SAVEFORCING
          write(400,'(1i8, 8e18.8)') ICOUNTF, WSU2,WSV2,SWRAD2,NO3_2,NH4_2,PO4_2,SIO4_2,DIS_2
          flush(400)
#endif
!       -----POM UNITS-----
!
        WSU2 = WSU2*(-ONE)/RHO0
        WSV2 = WSV2*(-ONE)/RHO0
        SWRAD2=SWRAD2*(-ONE)/rcp
!
         DO K = 1,KB
             READ (20,REC=(ICOUNTF-1)*KB+K) SCLIM2(K)
             READ (15,REC=(ICOUNTF-1)*KB+K) TCLIM2(K)
         END DO

         DO K = 1,KB-1
             READ (19,REC=(ICOUNTF-1)*(KB-1)+K) ISM2(K)
         END DO
!
#ifdef SAVEFORCING
          write(401,'(1i8,40e18.8)') ICOUNTF, ISM2
          write(402,'(1i8,40e18.8)') ICOUNTF, SCLIM2
          write(403,'(1i8,40e18.8)') ICOUNTF, TCLIM2
          flush(401) ; flush(402) ; flush(403)
#endif
      END IF
!
      return
!
      end subroutine FORCING_MANAGER
 end module Forcing
  
