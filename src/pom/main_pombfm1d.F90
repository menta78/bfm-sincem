#include"cppdefs.h"

      PROGRAM MAIN
!
      use global_mem,ONLY: RLEN,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
!      use TimeModule,ONLY: start,outdelt,dtphys,ihotst
      use POM,ONLY: KB,Z,ZZ,DZ,DZZ,H,COR,ALAT,ALON,TINIZ,SINIZ,UMOL, &
          UF,U,UB,VF,VB,V,WT,TB,T,TF,SF,S,SB,Q2F,Q2,Q2B,Q2LF,Q2L &
          ,Q2LB,KM,KH,KQ,L,A,C,VH,VHP,RHO,WUSURF,WVSURF,WUBOT,WVBOT,&
          SWRAD,ISHIFT,INSHIFT,IDSHIFT,TSURF,SSURF,WTSURF,WSSURF,DTI,&
          SMOTH,iint,dzr,iend,time,time0,idays 
      use Service,ONLY: ISM, savef,pon1p,pon3n,pon4n,pon5s 
      use Mem
      implicit none
!     .. Parameters ..
      integer,parameter            :: namlst=10
      real, parameter              :: PI=3.1415927
      INTEGER,parameter            :: ilong=selected_int_kind(12)
      INTEGER,PARAMETER            :: N_COMP=KB-1 !
      integer                      :: intt
      INTEGER(ilong)               :: istart
      INTEGER(ilong)               :: icm1
      real(RLEN),parameter         :: rcp=4.187E6
!     .. Local Scalars ..
      REAL(RLEN)                   :: DT2,KAPPA,RATIODF,RATIOF,RATIONF,&
                                               DAYI
      REAL(RLEN),SAVE              :: SWRAD1,SWRAD2,WTSURF1,WTSURF2,&
                                      Qcorr1,QCORR2,SSS1,SSS2,TSS1,TSS2&
                                      ,WSU1,WSU2,WSV1,WSV2,SLUX1,SLUX2
      REAL(RLEN),SAVE              :: NO3_1,PO4_1,NH4_1,SIO4_1, &
                                      NO3_2,PO4_2,NH4_2,SIO4_2
      INTEGER(ilong)               :: IMAXDELT,IDIAGN
      INTEGER(ilong)               :: K,M,RLENGTH,ihotst,KL1,KL2
      INTEGER(ilong)               :: ICOUNTF,IDOUNTF,IFCHGE,IFDCHGE,&
                                      IFDINT,IFINT,IFNCHGE, IFNINT
!     .. Local Arrays ..
      REAL(RLEN),SAVE              :: SINIZ1(KB),SINIZ2(KB),TINIZ1(KB),&
                                      TINIZ2(KB),UG(KB),VG(KB)
      REAL(RLEN)                   :: SSTAR(KB),WSADV(KB),TSTAR(KB),&
                                      WTADV(KB)
      REAL(RLEN),dimension(N_COMP) :: ISM1,ISM2
      logical                      :: ex,op
      character (len=11)           :: acc
!     .. External Subroutines ..
      EXTERNAL DENS,PROFQ,PROFT,PROFU,PROFV
!     .. Intrinsic Functions ..
      INTRINSIC FLOAT,MOD,NINT,SIN
!     .. Data statements ..
      DATA KAPPA/0.40/
!
      NAMELIST /Params_POMBFM/ H,DTI,ALAT,IDIAGN,IDAYS,SMOTH,&
                               ihotst,UMOL,KL1,KL2,savef
!       open(namlst,file='params_POMBFM.nml',status='old',action='read',err=100)
!       read(namlst,nml=Params_POMBFM,err=102)
      OPEN(namlst,file='params_POMBFM.nml',status='old',action='read')
      READ(namlst,nml=Params_POMBFM)
      CLOSE(namlst)
!======================================

      iend = idays*24*3600/IFIX(dti)
      DT2  = 2.0*DTI
      dayi = 1.e0/86400.e0
      WRITE(6,*) 'dayi =',dayi
!
!=======================================================
!    OPEN forcing files from subroutine:
!===============================================================
      CALL opendat(WSU1,WSV1,SSS1,SWRAD1,WTSURF1,QCORR1,SINIZ1, &
                     ISM1,NO3_1,NH4_1,PO4_1,SIO4_1,TINIZ,SINIZ,Tiniz1)
!===============================================================
      CALL calcdepth(z,zz,dz,dzz,kb,KL1,KL2)
!
      DO k=1,N_COMP         
        dzr(K)=1./dz(k)
      END DO
!==================================================================
!
      IF (ihotst.eq.1) THEN !restart from fileinquire(IOLENGTH=rlength) SB

       WRITE(6,*) 'reading POM-restart-File... '
       READ(70) u,v,t,s,q2,q2l,ub,vb,tb,sb,q2b,q2lb,kh,km,kq,L,WUBOT &
                ,WVBOT,RHO 
!====================================================================
! WE SHOULD ADD BFM RESTART
!====================================================================
      ELSE IF (ihotst.eq.0) THEN ! generic restart

!====================================================================
!READING OF T AND S INITIAL CONDOTION SHOULD BE PUT INTO A SUBROUTINE
!===================================================================
       DO K = 1,KB
           READ (29,REC=K) SINIZ(K)
           READ (10,REC=K) TINIZ(K)
       END DO
!======================================================================

       DO K = 1,KB
          SB(K)   = SINIZ(K)
          S(K)    = SB(K)
          TB(K)   = TINIZ(K)
          T(K)    = TB(K)
          UB(k)   = 0.
          VB(k)   = 0.
          U(K)    = 0.
          V(K)    = 0.
          Q2B(K)  = 1.E-7
          Q2(K)   = Q2B(K)
          Q2LB(K) = 1.E-7
          Q2L(K)  = Q2LB(K)
          KM(K)   = 0.
          KH(K)   = 0.
          KQ(K)   = 0.
          L(K)    = Q2L(K)/Q2(K)
          L(1)    = 0.
          L(KB)   = 0.
          UF(K)   = U(K)
          VF(K)   = V(K)
          UG(K)   = 0.
          VG(K)   = 0.
          WT(K)   = 0.
          TF(K)   = T(K)
          SF(K)   = S(K)
          Q2F(k)  = Q2(K)
          Q2LF(K) = Q2L(K)
          VH(K)   = 0.
          VHP(K)  = 0.
          A(K)    = 0.
          C(K)    = 0.
       ENDDO
          WUBOT   = 0.
          WVBOT   = 0.
          WTSURF  = 0.
          SWRAD   = 0.

      ENDIF
      COR = 1.e-4*SIN(ALAT*2*PI/360.)
      IDOUNTF = 0 ! Day counter                                                         !    (G)

!
!
!==============================================================
! COMPUTE # OF TIME STEPS NEEDED TO COVER ONE DAY AND ONE MONTH
!==============================================================
!      IFCHGE AND IFINT REFER TO MONTHLY FILES                         
!      IFDCHGE AND IFDINT REFER TO DAYLY FILES                         
!==============================================================
      IFDCHGE = (3600*24)/NINT(DTI) ! No. of timesteps per day
      IFDINT  = -1
      IFCHGE  = (3600*24*30)/NINT(DTI) ! No. of timesteps per month

      IF(mod(IDOUNTF,30).le.14) THEN
      IFINT   = mod(IDOUNTF,30)*(3600*24)/NINT(DTI)+(3600*24*15)/&
                NINT(DTI)-1
      ELSE
      IFINT   = mod(IDOUNTF,30)*(3600*24)/NINT(DTI)-(3600*24*15)/&
                NINT(DTI)-1
      ENDIF
!
!========================================================================
!     Define initial density field                c
      CALL DENS(T,S,ZZ,H,RHO,KB)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     INITIAL READING OF SURFACE FORCING                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     SHIFT THE FORCING FILES TO START AT THE DESIRED MONTH            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ICOUNTF = (IDOUNTF+15)/30 + 1 ! month counter 
      IDOUNTF = IDOUNTF + 1         ! day counter + 1 
      IF(ICOUNTF==13) ICOUNTF=1
      IF(IDOUNTF==361) IDOUNTF=1

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      INITIAL  WIND STRESS READING                                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      READ (11,REC=ICOUNTF) WSU1,WSV1
      READ (11,REC=ICOUNTF+1) WSU2,WSV2
!     ------CONVERT TO POM UNITS---------
!           from N/m2 to m2/s2
      WSU1 = -WSU1*1.e-3
      WSU2 = -WSU2*1.e-3
      WSV1 = -WSV1*1.e-3
      WSV2 = -WSV2*1.e-3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      READ INITIAL SST-SSS                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      READ (13,REC=ICOUNTF) SSS1
      READ (13,REC=ICOUNTF+1) SSS2
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      READ INITIAL HEAT FLUX
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc`
!
      READ(21,REC=ICOUNTF) SWRAD1,WTSURF1,QCORR1
      READ(21,REC=ICOUNTF+1) SWRAD2,WTSURF2,QCORR2
!     -----ELIMINATE SW RADIATION FROM TOTAL HEAT FLUX-------
      WTSURF1 = WTSURF1-SWRAD1+QCORR1
      WTSURF2 = WTSURF2-SWRAD2+QCORR2
!      ------ CONVERT TO POM UNITS---------
      SWRAD1  = -SWRAD1/rcp
      SWRAD2  = -SWRAD2/rcp
      WTSURF1 = -WTSURF1/rcp
      WTSURF2 = -WTSURF2/rcp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      READ INITIAL TINIZ-SINIZ PROFILES                                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO K = 1,KB
          READ (20,REC=(ICOUNTF-1)*KB+K) SINIZ1(K)
          READ (15,REC=(ICOUNTF-1)*KB+K) TINIZ1(K)
      END DO
      DO K = 1,KB
          READ (20,REC=ICOUNTF*KB+K) SINIZ2(K)
          READ (15,REC=ICOUNTF*KB+K) TINIZ2(K)
      END DO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      READ INITIAL SUSPENDED SEDIMENT  PROFILES                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO K = 1,N_COMP
          READ (19,REC=(ICOUNTF-1)*N_COMP+K) ISM1(K)
      END DO
      DO K = 1,N_COMP
          READ (19,REC=ICOUNTF*N_COMP+K) ISM2(K)
      END DO 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      READ INITIAL SURFACE NUTRIENTS                                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      READ(18,REC=ICOUNTF) NO3_1,NH4_1,PO4_1,SIO4_1
      READ(18,REC=ICOUNTF+1) NO3_2,NH4_2,PO4_2,SIO4_2
!     
      ICOUNTF = ICOUNTF + 1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      READ SOLAR RADIATION AT SURFACE                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      READ (17,REC=IDOUNTF) SLUX1
      READ (17,REC=IDOUNTF+1) SLUX2
!
      IDOUNTF = IDOUNTF + 1
!
!***********************************************************************
!       BFM SET UP (include initial conditions definition)
!***********************************************************************
#ifndef POM_only
!-----------------------------------!
!   initialization of BFM           !
!-----------------------------------!
      CALL pom_ini_bfm_1d
#endif
      time0=0.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      BEGIN TIME MARCH
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      WRITE(6,*) 'Icount before time march loop = ',ICOUNTF
!      print *,'params',H,DTI,ALAT,IDIAGN,IDAYS,SMOTH,IHOTST,UMOL,KL1,KL2,SAVEF
      istart = 1  
      
      DO intt=ISTART, IEND                       ! (G)
       IF (intt.EQ.IEND) GO TO 9000              ! (G)
!      DO 9000 intt = ISTART, IEND               ! (G)
!      iint = intt                               ! (G)
!
      time = time0+(DTI*float(intt)*dayi)       ! (G)
 
!     -----------INTERPOLATION COUNTERS---------
      IFINT  = IFINT + 1
      RATIOF = FLOAT(IFINT)/FLOAT(IFCHGE)
!
      IFDINT  = IFDINT + 1
      RATIODF = FLOAT(IFDINT)/FLOAT(IFDCHGE)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     INTERPOLATE  MONTHLY WIND STRESS                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      WUSURF = WSU1 + RATIOF* (WSU2-WSU1)
      WVSURF = WSV1 + RATIOF* (WSV2-WSV1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      INTERPOLATE MONTHLY SST-SSS                                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      TSURF = TINIZ1(1) + RATIOF* (TINIZ2(1)-TINIZ1(1))
      SSURF = SINIZ1(1) + RATIOF* (SINIZ2(1)-SINIZ1(1))
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      INTERPOLATE HEAT FLUX
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      SWRAD  = SWRAD1+RATIOF*(SWRAD2-SWRAD1)
      WTSURF = WTSURF1+RATIOF*(WTSURF2-WTSURF1)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      INTERPOLATE SALINITY PROFILES
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SSTAR = SINIZ1+RATIOF*(SINIZ2-SINIZ1)
      TSTAR = TINIZ1+RATIOF*(TINIZ2-TINIZ1)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      INTERPOLATE ISM PROFILES
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ISM = ISM1+RATIOF*(ISM2-ISM1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      INTERPOLATE MONTHLY TINIZ-SINIZ PROFILES                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (IDIAGN.EQ.1) THEN
        DO K = 1,KB
             TF(K) = TINIZ1(K) + RATIOF* (TINIZ2(K)-TINIZ1(K))
             SF(K) = SINIZ1(K) + RATIOF* (SINIZ2(K)-SINIZ1(K))
        END DO
      END IF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      INTERPOLATE SURFACE NUTS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pon3n = NO3_1+RATIOF*(NO3_2-NO3_1)
      pon4n = NH4_1+RATIOF*(NH4_2-NH4_1)
      pon1p = PO4_1+RATIOF*(PO4_2-PO4_1)
      pon5s = SIO4_1+RATIOF*(SIO4_2-SIO4_1)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            MY Turbulence closure     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO K = 1,KB
          Q2F(K)  = Q2B(K)
          Q2LF(K) = Q2LB(K)
      END DO

      CALL PROFQ(DT2)
!
!     ----MIXING TIME STEP (ASSELIN)-----
!
      DO K = 1,KB
          Q2(K)   = Q2(K) + .5*SMOTH* (Q2F(K)+Q2B(K)-2.*Q2(K))
          Q2B(K)  = Q2(K)
          Q2(K)   = Q2F(K)
          Q2L(K)  = Q2L(K) + .5*SMOTH* (Q2LF(K)+Q2LB(K)-2.*Q2L(K))
          Q2LB(K) = Q2L(K)
          Q2L(K)  = Q2LF(K)
      END DO

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     SKIP THE T AND S CALC> IF PROFILES PROVIDED BY DATA              c
!     (IDIAGN=1)                                                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      IF(IDIAGN.EQ.0 ) THEN
         DO K = 1,KB
            TF(K) = TB(K)
            SF(K) = SB(K)
         END DO
!
!         ------TEMPERATURE------
!
!         CALL PROFT1D(TF,TB,WTSURF,SWRAD,TSURF,1,DT2)
         CALL PROFT1D(TF,TB,WTSURF,SWRAD,TSURF,3,DT2)

!
!         ----SALINITY INCLUDING LATERAL ADVECTION TERM-----
!
         wssurf = -(ssurf-S(1))*5.68/(86400.)
         wsadv = 0.
         DO k=1,KB
           IF(-ZZ(k)*H.ge.1.) THEN
             wsadv(k) = -(sstar(k)-s(k))/(1.*86400.)
           ENDIF
         ENDDO
!         CALL PROFS(SF,SB,WSSURF,WSADV,0.,SSURF,1,DT2)
        CALL PROFS(SF,SB,WSSURF,WSADV,0.,SSURF,3,DT2)
!
!     ----MIXING TIME STEP (ASSELIN)-----
!
         DO K = 1,KB
             T(K)  = T(K) + .5*SMOTH* (TF(K)+TB(K)-2.*T(K))
             TB(K) = T(K)
             T(K)  = TF(K)
             S(K)  = S(K) + .5*SMOTH* (SF(K)+SB(K)-2.*S(K))
             SB(K) = S(K)
             S(K)  = SF(K)
         END DO
      ELSE
         DO K = 1,KB
            TB(K)=T(K)
            SB(K)=S(K)
            T(k)=TF(K)
            S(K)=SF(K)
         END DO
      END IF
!
!     -----NEW DENSITY-----
!
      CALL DENS(T,S,ZZ,H,RHO,KB)
!
!     ----------CALL THE BFM-------
!
#ifdef  POM_only
        if(mod(intt,864).eq.0) then
        do k= 1, kb
        WRITE(102,*) time,zz(k),TB(k),SB(k),RHO(k)
        enddo
        endif
#endif

#ifndef POM_only
!---------------------------------------------!
!  Call BFM after physics has been calculated ! 
!---------------------------------------------!
      CALL pom_bfm_1d

#endif
!     -----VELOCITY-----
!
      DO K = 1,KB - 1
          UF(K) = UB(K) + DT2*COR* (V(K)-VG(K))
          VF(K) = VB(K) - DT2*COR* (U(K)-UG(K))
      END DO
!
      CALL PROFU(DT2)
      CALL PROFV(DT2)
!     ----MIXING TIME STEP (ASSELIN)-----
!
      DO K = 1,KB
          U(K)  = U(K) + .5*SMOTH* (UF(K)+UB(K)-2.*U(K))
          V(K)  = V(K) + .5*SMOTH* (VF(K)+VB(K)-2.*V(K))
          UB(K) = U(K)
          U(K)  = UF(K)
          VB(K) = V(K)
          V(K)  = VF(K)
      END DO
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     READ NEW WIND STRESS, SURFACE FLUXES  AND RIVER RUNOFF EVERY MONTH
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (IFINT.EQ.IFCHGE) THEN
          PRINT *,'inside loop: read new wind stress, surface fluxes &
                   and river runoff every month'
!         print*,'PO4_1,PO4_2',PO4_1,PO4_2   
          ICOUNTF = ICOUNTF + 1
          IFINT   = 0
          PRINT *,'icountf =',ICOUNTF
          WSU1    = WSU2
          WSV1    = WSV2
          SWRAD1  = SWRAD2
          QCORR1  = QCORR2
          WTSURF1 = WTSURF2
          NO3_1   = NO3_2
          NH4_1   = NH4_2
          PO4_1   = PO4_2
          SIO4_1  = SIO4_2
!          TSS1 = TSS2
          SSS1    = SSS2
          SINIZ1  = SINIZ2
          ISM1    = ISM2
         DO K = 1,KB
             TINIZ1(K) = TINIZ2(K)
             SINIZ1(K) = SINIZ2(K)
         END DO

         IF (ICOUNTF.GT.13) THEN
             ICOUNTF = 2
             READ (11,REC=1) WSU1,WSV1
             WSU1 = -WSU1*1.e-3
             WSV1 = -WSV1*1.e-3
             READ (13,REC=1) SSS1
             READ(21,REC=1) SWRAD1,WTSURF1,QCORR1
             WTSURF1=WTSURF1-SWRAD1+QCORR1
             WTSURF1=-WTSURF1/rcp
             SWRAD1=-SWRAD1/rcp

             READ(18,REC=1) NO3_1,NH4_1,PO4_1,SIO4_1
 
             DO K = 1,KB
                 READ (20,REC=K) SINIZ1(K)
                 READ (15,REC=K) TINIZ1(K)
             END DO
             DO K = 1,N_COMP
                  READ (19,REC=K) ISM1(K)
             END DO
         END IF

         READ(18,REC=ICOUNTF) NO3_2,NH4_2,PO4_2,SIO4_2
         READ (11,REC=ICOUNTF) WSU2,WSV2
         WSU2 = -WSU2*1.e-3
         WSV2 = -WSV2*1.e-3
         READ (13,REC=ICOUNTF) SSS2
         READ (21,REC=ICOUNTF) SWRAD2,WTSURF2,QCORR2
         WTSURF2=WTSURF2-SWRAD2+QCORR2
         WTSURF2=-WTSURF2/rcp
         SWRAD2=-SWRAD2/rcp
         icm1=ICOUNTF-1
          
         DO K = 1,KB
             READ (20,REC=icm1*KB+K) SINIZ2(K)
             READ (15,REC=icm1*KB+K) TINIZ2(K)
         END DO
         
         DO K = 1,N_COMP
             READ (19,REC=icm1*N_COMP+K) ISM2(K)
         END DO
      END IF
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     READ NEW SOLAR RADIATION EVERY DAY                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      IF (IFDINT.EQ.IFDCHGE) THEN
          IDOUNTF = IDOUNTF + 1
          IFDINT = 0.
          SLUX1 = SLUX2
          IF (IDOUNTF.GT.361) THEN
              IDOUNTF = 2
              READ (17,REC=1) SLUX1
          END IF
          READ (17,REC=IDOUNTF) SLUX2
      END IF

      ENDDO                         !   (G)
!
 9000 CONTINUE
!
!100   call error_msg_prn(NML_OPEN,"main_POMBFM.F90","params_POMBFM.nml")
!102   call error_msg_prn(NML_READ,"main_POMBFM.F90","Params_POMBFM")

       PRINT *,'main done'
! 105 FORMAT(f6.1,4f15.2)

      END
