C $Header: MITgcm/pkg/BFMcoupler/BFMcoupler_calc_tendency.F,v 1.01
C 2014/04/10  cgp
C modified 4 sept 2015 by Valeria, introducing N1p, N3n, N5s, O3c, O3h
C contribution by rivers
C in sept 2017 for the bottom contributions
C in feb 2019 for the daily PAR modulation
CBOP
C !ROUTINE: BFMcoupler_calc_tendency

#include "BFMcoupler_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_BFMCOUPLER



C !INTERFACE: ==========================================================
      SUBROUTINE BFMCOUPLER_CALC_TENDENCY(
     I                bi,bj,imin,imax,jmin,jmax,
     I                myIter,myTime,myThid)
C !DESCRIPTION:
C passes Tr and environmental variables to the BFM library and get the gchemTendency
C due to biological activity

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "BFMcoupler_VARS.h" ! contains #include "BFM_var_list.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM_FIELDS.h"
#include "EXF_FIELDS.h"
#ifdef USE_QSW
#include "FFIELDS.h"
#endif

#ifdef ALLOW_LONGSTEP
#include "LONGSTEP.h"
#endif


C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C

C     === Global variables ===
       INTEGER bi,bj,imin,imax,jmin,jmax
       INTEGER myIter
       _RL myTime
       INTEGER myThid

C     === Local variables ====
       _RL atten
       _RL parconv,parfrac
       _RL ligth     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef USE_SINK
       _RL zwork(21,Nr)
       _RL Fsink(21,Nr)
#endif /* USE_SINK */
       INTEGER i,j,k,jtr,ktmp,iv

#ifdef BFMv2
       _RL a(jptra),er(10),b(jptra),c(4),d(jptra_dia),d2(1)
C      _RS AtmosWIND0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
       logical bot,sur
#else
       _RL a(Nr,jptra)
       _RL er(Nr,11)
       _RL b(jptra,Nr)
       _RL c(4,Nr)
       _RL d(jptra_dia,Nr)
       _RL d2(jptra_dia_2d)
#endif

       INTEGER kBot(1-OLx:sNx+OLx,1-OLy:sNy+OLy), kB
C     msgBuf    :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER stdUnit

C     == routine arguments ==
       INTEGER mydate(4)
       INTEGER mydateday(4)
C     == local variables ==
       INTEGER mydate0(4)
       _RL     daysecs
       double precision POSITIVE_VALUE

C ==== LOCAL VARIABLE FOR DIAGNOSTICS ====
#ifdef ALLOW_DIAGNOSTICS
#include "BFMcoupler_VARDIAGlocal.h"
#endif /* ALLOW_DIAGNOSTICS */


CEOP
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      stdUnit = standardMessageUnit

#ifdef ALLOW_DIAGNOSTICS
c initialized variable for diagnostics
#include "BFMcoupler_VARDIAGinitializ.h"
#endif /* ALLOW_DIAGNOSTICS */



C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c initialization of kBot 

c      DO bj = myByLo(myThid), myByHi(myThid)
c       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=jmin,jmax
         DO i=imin,imax
             kBot(i,j)   = 1
         ENDDO
        ENDDO
c       ENDDO
c      ENDDO


c compute par using QSW or PAR read from files or spar already
c initialized in BFMcoupler_ini_forcing.F
       parconv= 1. _d 0/0.2174 _d 0  !conversion from W/m2 to uEin/m2/s
       parfrac=  0.4 _d 0            !fraction Qsw that is PAR

      DO j=jmin,jmax
       DO i=imin,imax
#ifdef READ_PAR
        ligth(i,j)=spar(i,j,bi,bj)*maskC(i,j,1,bi,bj)
#elif (defined USE_QSW)
#ifdef ALLOW_LONGSTEP    
        ligth(i,j)=-LS_Qsw(i,j,bi,bj)*parconv*parfrac*maskC(i,j,1,bi,bj)
#else
        ligth(i,j)=-Qsw(i,j,bi,bj)*parconv*parfrac*maskC(i,j,1,bi,bj)
#endif /* ALLOW_LONGSTEP */
#else
        ligth(i,j)=spar(i,j,bi,bj)*maskC(i,j,1,bi,bj)
c    spar already initialized in BFMcouple_ini_forcing.F at 250 uEin/m2/s
#endif
       ENDDO
      ENDDO

c For day-night cycle
c First call to extract the date array
c mydate: yyyymmdd,hhmmss,leap_year,dayOfWeek
       CALL cal_getdate(myIter,myTime,mydate,myThid)
c Then computation of the day time in seconds: daysecs
       mydate0=mydate
       mydate0(2)=0
       mydateday=mydate-mydate0
       mydateday(3)=0   
       mydateday(4)=-1 
       CALL cal_toseconds(mydateday,daysecs,myThid)
c Cosine modulation
       POSITIVE_VALUE=cos(2*PI*DBLE(daysecs)/86400.-PI)


c atten coeff read from file xESP, or computed as shelf shading of chla
      DO j=jmin,jmax
       DO i=imin,imax
             atten=0. _d 0 ! 2 steps for light attenation since half
             DO k=1,Nr     ! depth of cells is used
                atten=atten+ ( 
#ifdef USE_SHADE
#ifdef READ_xESP
     &   xESP(i,j,bi,bj)  
#else 
     &   p_espWATER
#endif
c TO DO: should be generalized!
     &   + p_espP1*pTracer(i,j,k,bi,bj,14) ! P1i diatom
     &   + p_espP2*pTracer(i,j,k,bi,bj,19) ! P2i flag
     &   + p_espP3*pTracer(i,j,k,bi,bj,23) ! P3i picopl
     &   + p_espP4*pTracer(i,j,k,bi,bj,27) ! P4i dinofl
#elif (defined (READ_xESP))
     &   xESP(i,j,bi,bj) 
#else
     &   0. _d 0 ! no attenuation  :: model is stopped in readparams.F
#endif
     &    )*.5 _d 0 *drF(k)

                if (k.gt.1) then
                     atten=atten+ ( 
#ifdef USE_SHADE
#ifdef READ_xESP
     &   xESP(i,j,bi,bj)            
#else 
     &   p_espWATER
#endif
c TO DO: should be generalized!
     &   + p_espP1*pTracer(i,j,k-1,bi,bj,14) ! P1l diatom
     &   + p_espP2*pTracer(i,j,k-1,bi,bj,19) ! P2l flag
     &   + p_espP3*pTracer(i,j,k-1,bi,bj,23) ! P3l picopl
     &   + p_espP4*pTracer(i,j,k-1,bi,bj,27) ! P4l dinofl
#elif (defined (READ_xESP))
     &     xESP(i,j,bi,bj)
#else
     &   0. _d 0 ! no attenuation :: model is stopped in readparams.F
#endif
     &    )*.5 _d 0 *drF(k-1)
                endif

#ifdef USE_CYCLE
#ifdef READ_PAR
                xpar(i,j,k,bi,bj)=(max(1.0d-3,
     & PI*POSITIVE_VALUE*ligth(i,j)))*exp(-atten)
#else 
                xpar(i,j,k,bi,bj)=ligth(i,j)*exp(-atten)
#endif
#else 
                xpar(i,j,k,bi,bj)=ligth(i,j)*exp(-atten)
#endif
             ENDDO
          ENDDO
         ENDDO

#ifdef BFMv2

c      DO j=1-OLy,sNy+OLy
c      DO i=1-OLx,sNx+OLx
      DO k=1,Nr
      DO j=jmin,jmax
      DO i=imin,imax
            do jtr=1,jptra_dia 
               d(jtr)=0.0   ! initialized BFM diagnostics to 0 
            enddo
            do jtr=1,4
               c(jtr)=0.0   ! initialized BFM sinking velocity to 0 
            enddo


            if (HFacC(i,j,k,bi,bj).gt.0. _d 0) then

             DO jtr=1, jptra
                 b(jtr)=0.0
                 if (pTracer(i,j,k,bi,bj,jtr).gt.0) then
                   a(jtr) = pTracer(i,j,k,bi,bj,jtr)
                 else
                   a(jtr)=1.0 _d -10 ! in the case of negative values
                 endif
             END DO
c set bot=1.0 if the layer below is not ocean
             ktmp=min(nR,k+1)
             if(hFacC(i,j,ktmp,bi,bj).eq.0. _d 0.or.k.eq.Nr) then
               bot =.TRUE.
               kBot(i,j)=ktmp-1
             else
               bot =.FALSE.
             endif
c set sur = 1.0 if the layer is equal 1 (the land is already checked)
             if (k.eq.1) then
                sur=.TRUE.
             else
                sur=.FALSE.
             endif
! Environmental regulating factors (er)
#ifdef ALLOW_LONGSTEP
             er(1)  = min(30.,LS_theta(i,j,k,bi,bj))  ! from LONGSTEP.h
             er(2)  = LS_salt(i,j,k,bi,bj)   ! from LONGSTEP.h
             er(3)  = rhoConst + rhoInSitu(i,j,k,bi,bj) ! rhoInSitu (density anomaly [kg/m^3] at cell center level) +
                                                        ! + RhoConst (set equal to rhoNil=999.8 (set_defaults.F and Ini_parms.F))
#else
             er(1)  = min(30.,theta(i,j,k,bi,bj))  ! from DYNVARS.h
             er(2)  = salt(i,j,k,bi,bj)   ! from DYNVARS.h
             er(3)  = rhoConst + rhoInSitu(i,j,k,bi,bj) ! rhoInSitu (density anomaly [kg/m^3] at cell center level) +
                                                        ! + RhoConst (set equal to rhoNil=999.8 (set_defaults.F and Ini_parms.F))
#endif /* ALLOW_LONGSTEP  */
             er(4)  = 0.0                 ! NO ICE is implemented
             er(5)  = AtmosPCO2(i,j,bi,bj)! PCO2 atmopheric from BFMcoupler_VAR.h
             er(6)  = max(1.0,xpar(i,j,k,bi,bj))   ! PAR [uEin/s]  ATTENTION DAILY CYCLE 
             er(7)  = 24.0                ! 24 hours
             er(8)  =  drF(k)*hFacC(i,j,k,bi,bj)! H(m) of cell considering the fract of grid box above topography
             er(9)  = wspeed(i,j,bi,bj)   ! wind speed [m/s] from EXF_FIELDS.H
              if (er(9).LE.0.0) then
                    er(9)=0.001
              endif
             er(10) = PH(i,j,k,bi,bj)     ! from BFMcoupler_VAR.h

#ifdef BFMcoupler_DEBUG
c ===== CHECK DEBUGING
            if (i.eq.8.and.j.eq.8.and.k.eq.1) then
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'BFM1:'
     &                 ,er(1),er(2),er(3),er(4),er(5)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'BFM2:'
     &                 ,er(6),er(7),er(8),er(9),er(10)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'BFM3:'
     &                 ,a(1),a(2),a(3),a(4),a(5)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'BFM4:'
     &                 ,a(6),a(7),a(8),a(9),a(10)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
              endif
#endif /* endif BFMcoupler_debug */


c  ===== CALL TO BFM ====
             CALL BFM0D_Input_EcologyDynamics(sur,bot,a,jptra,er)
             CALL BFM0D_reset()
             CALL EcologyDynamics()
             CALL BFM0D_Output_EcologyDynamics_surf(b,c,d,d2)
c            CALL BFM0D_Output_EcologyDynamics(b,jptra,c,d,jptra_dia-1)
c  ===== return of BFM tendency ====
             DO jtr=1, jptra
                 gchemTendency(i,j,k,bi,bj,jtr) = 
     &                 gchemTendency(i,j,k,bi,bj,jtr)+b(jtr) ! tendency [mmol/m3/s]
             END DO


#ifdef BFMcoupler_DEBUG
c ===== CHECK DEBUGING
            DO jtr=1,jptra
            if (isnan2(b(jtr))) then
       WRITE(msgBuf,'(A,I2,1x,I3,1x,I3,1x,I3,1x,I3,1x,I3,1x,L,1x,L)'),
     &               'NAN:',jtr,i,j,k,bi,bj,sur,bot
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'ER1:'
     &                 ,er(1),er(2),er(3),er(4),er(5)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'ER2:'
     &                 ,er(6),er(7),er(8),er(9),er(10)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A1:'
     &                 ,a(1),a(2),a(3),a(4),a(5)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A2:'
     &                 ,a(6),a(7),a(8),a(9),a(10)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A3:'
     &                 ,a(11),a(12),a(13),a(14),a(15)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A4:'
     &                 ,a(16),a(17),a(18),a(19),a(20)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A5:'
     &                 ,a(21),a(22),a(23),a(24),a(25)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A6:'
     &                 ,a(26),a(27),a(28),a(29),a(30)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A7:'
     &                 ,a(31),a(32),a(33),a(34),a(35)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A8:'
     &                 ,a(36),a(37),a(38),a(39),a(40)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A9:'
     &                 ,a(41),a(42),a(43),a(44),a(45)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'A10:'
     &                 ,a(46),a(47),a(48),a(49),a(50)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3)'),'A11:',a(51)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)

       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B1:'
     &                 ,b(1),b(2),b(3),b(4),b(5)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B2:'
     &                 ,b(6),b(7),b(8),b(9),b(10)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B3:'
     &                 ,b(11),b(12),b(13),b(14),b(15)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B4:'
     &                 ,b(16),b(17),b(18),b(19),b(20)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B5:'
     &                 ,b(21),b(22),b(23),b(24),b(25)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B6:'
     &                 ,b(26),b(27),b(28),b(29),b(30)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B7:'
     &                 ,b(31),b(32),b(33),b(34),b(35)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B8:'
     &                 ,b(36),b(37),b(38),b(39),b(40)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B9:'
     &                 ,b(41),b(42),b(43),b(44),b(45)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'B10:'
     &                 ,b(46),b(47),b(48),b(49),b(50)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
       WRITE(msgBuf,'(A,E10.3)'),'B11:',b(51)
       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)

       STOP 'NAN in BFM tendency: ABNORMAL ENDED' 
            endif
            enddo



c            if (i.eq.8.and.j.eq.8.and.k.eq.1) then
c       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'BFM5:'
c     &                 ,b(1),b(2),b(3),b(4),b(5)
c       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
c       WRITE(msgBuf,'(A,E10.3,E10.3,E10.3,E10.3,E10.3)'),'BFM6:'
c     &                 ,b(6),b(7),b(8),b(9),b(10)
c       CALL PRINT_MESSAGE(msgBuf,stdUnit,SQUEEZE_RIGHT,myThid)
c             endif
#endif  /* BFMcoupler_debug */

c  ===== sinking velocity for 4 phyto
             DO jtr=1,4
                   sediPI(i,j,k,bi,bj,jtr) = c(jtr) ! sedimentation velocities of 4 phytos
             END DO
             

             PH(i,j,k,bi,bj)=d(9) ! stored localy since needed for the Follows solver guess 

#ifdef ALLOW_DIAGNOSTICS
#include "BFMcoupler_VARDIAGcopy_fromD.h"
#endif /* ALLOW_DIAGNOSTICS */

             ENDIF ! endif HFacC>0
c NEDDED to be coded
c             CALL burial
c             CALL sediment_flux

c  ====== save BFM diagnostics in MITgcm diagnostic pkg 
      
         ENDDO   ! enddo i cycle
         ENDDO   ! enddo j cycle
         ENDDO ! enddo k cycle

#else

c      DO j=1-OLy,sNy+OLy
c      DO i=1-OLx,sNx+OLx
      DO j=jmin,jmax
      DO i=imin,imax

          if (HFacC(i,j,1,bi,bj).gt.0. _d 0) then

             kBot(i,j) = bottom(i,j,bi,bj,Nr) ! bottom layer is computed here

             a = 1.0 ! default initialization

             do k=1,kBot(i,j)
                 DO jtr=1,jptra
                     if (pTracer(i,j,k,bi,bj,jtr).gt.0) then
                         a(k,jtr) = pTracer(i,j,k,bi,bj,jtr)
                     else
                         a(k,jtr)=1.0 _d -10 ! in the case of negative values
                     endif
                 END DO
             enddo

             b = 0.0 ! default initialization
             c = 0.0 ! initialized BFM sinking velocity to 0
             d = 0.0 ! initialized BFM diagnostics to 0
             d2 = 0.0 ! initialized BFM diagnostics to 0
             er = 1.0 ! default initialization

! Environmental regulating factors (er)
#ifdef ALLOW_LONGSTEP
             do k=1,kBot(i,j)
                 er(k,1) = min(30.,LS_theta(i,j,k,bi,bj)) ! from LONGSTEP.h
             enddo
             er(1:kBot(i,j),2) = LS_salt(i,j,1:kBot(i,j),bi,bj) ! from LONGSTEP.h
             er(1:kBot(i,j),3) =
     &           rhoConst + rhoInSitu(i,j,1:kBot(i,j),bi,bj) ! rhoInSitu (density anomaly [kg/m^3] at cell center level) +
                                                             ! + RhoConst (set equal to rhoNil=999.8 (set_defaults.F and Ini_parms.F))
#else
             do k=1,kBot(i,j)
                 er(k,1)  = min(30.,theta(i,j,k,bi,bj)) ! from DYNVARS.h
             enddo
             er(1:kBot(i,j),2) = salt(i,j,1:kBot(i,j),bi,bj) ! from DYNVARS.h
             er(1:kBot(i,j),3) =
     &           rhoConst + rhoInSitu(i,j,1:kBot(i,j),bi,bj) ! rhoInSitu (density anomaly [kg/m^3] at cell center level) +
                                                             ! + RhoConst (set equal to rhoNil=999.8 (set_defaults.F and Ini_parms.F))
#endif /* ALLOW_LONGSTEP  */
             er(1,4) = 0.0 ! NO ICE is implemented
             er(1,5) = AtmosPCO2(i,j,bi,bj) ! PCO2 atmopheric from BFMcoupler_VAR.h
             do k=1,kBot(i,j)
                 er(k,6) = max(1.0,xpar(i,j,k,bi,bj)) ! PAR [uEin/s]  ATTENTION DAILY CYCLE 
             enddo
             er(1,7) = 24.0 ! 24 hours
             do k=1,kBot(i,j)
                 er(k,8) = drF(k)*hFacC(i,j,k,bi,bj) ! H(m) of cell considering the fract of grid box above topography
             enddo
             er(1,9) = wspeed(i,j,bi,bj) ! wind speed [m/s] from EXF_FIELDS.H
              if (er(1,9).LE.0.0) then
                    er(1,9)=0.001
              endif
             er(1:kBot(i,j),10) = PH(i,j,1:kBot(i,j),bi,bj) ! from BFMcoupler_VAR.h
             er(1:kBot(i,j),11) = 0.0

c  ===== CALL TO BFM ====
             CALL BFM1D_Input_EcologyDynamics(kBot(i,j),a,jptra,er)
             CALL BFM1D_reset()
             CALL EcologyDynamics()
             CALL BFM1D_Output_EcologyDynamics(b,c,d,d2)
c  ===== return of BFM tendency ====
             DO jtr=1, jptra
                 gchemTendency(i,j,1:kBot(i,j),bi,bj,jtr) = 
     &               gchemTendency(i,j,1:kBot(i,j),bi,bj,jtr) +
     &               b(jtr,1:kBot(i,j)) ! tendency [mmol/m3/s]
             END DO

c  ===== sinking velocity for 4 phyto
             DO jtr=1,4
                   sediPI(i,j,1:kBot(i,j),bi,bj,jtr) =
     &                 c(jtr,1:kBot(i,j)) ! sedimentation velocities of 4 phytos
             END DO

             PH(i,j,1:kBot(i,j),bi,bj)=d(9,1:kBot(i,j)) ! stored localy since needed for the Follows solver guess

#ifdef ALLOW_DIAGNOSTICS
#include "BFMcoupler_VARDIAGcopy_fromD.h"
#endif /* ALLOW_DIAGNOSTICS */

          endif ! endif HFacC>0

c NEDDED to be coded
c             CALL burial
c             CALL sediment_flux

c  ====== save BFM diagnostics in MITgcm diagnostic pkg 
      
      ENDDO   ! enddo i cycle
      ENDDO   ! enddo j cycle

#endif

c =====  SINKING ===================
#ifdef USE_SINK
      DO j=jmin,jmax
       DO i=imin,imax
          do k=1,Nr
            do iv=1,21
              zwork(iv,k)=0.0
            enddo
          enddo
c computing sinking at the open downward cell faces 
          DO k=1,Nr-1
            if (HFacC(i,j,k+1,bi,bj).gt.0. _d 0) then
      do iv=1,4  ! Particulates Tr45-48
       zwork(0+iv,k)=vsed*pTracer(i,j,k,bi,bj,44+iv)
      enddo
      do iv=1,5  ! Diatoms: Tr11-Tr1
      zwork(4+iv,k)=sediPI(i,j,k,bi,bj,1)*pTracer(i,j,k,bi,bj,10+iv)
      enddo
      do iv=1,4  ! Flagellates: Tr16-Tr19
       zwork(9+iv,k)=sediPI(i,j,k,bi,bj,2)*pTracer(i,j,k,bi,bj,15+iv)
      enddo
      do iv=1,4  ! Picophytoplankton: Tr20-Tr23 
      zwork(13+iv,k)=sediPI(i,j,k,bi,bj,3)*pTracer(i,j,k,bi,bj,19+iv)
      enddo
      do iv=1,4  ! Dinoflagellates: Tr24-Tr27
      zwork(17+iv,k)=sediPI(i,j,k,bi,bj,4)*pTracer(i,j,k,bi,bj,23+iv)
      enddo
            endif
          enddo ! enddo k
c compute sinking fluxes for the cells
      do iv=1,21
        do k=1,Nr
           Fsink(iv,k)=0.0
        enddo
        if  (HFacC(i,j,1,bi,bj).gt.0. _d 0) then
c        Fsink(iv,1)=-zwork(iv,1)*
c     &         1/(drF(1)*hFacC(i,j,1,bi,bj)+EtaN(i,j,bi,bj)) ! 2014-08-20 no eta for mass conservation
         Fsink(iv,1)=-zwork(iv,1)*recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        endif
        do k=2,Nr
           if (HFacC(i,j,k,bi,bj).gt.0. _d 0) then
             Fsink(iv,k)=(zwork(iv,k-1)-zwork(iv,k)) *
     &             recip_drF(k)*recip_hFacC(i,j,k,bi,bj)
           endif
        enddo
      enddo
c update tendency with sinking fluxes
         do k=1,Nr
           if (HFacC(i,j,k,bi,bj).gt.0. _d 0) then
             do iv=1,4 ! Particulates  Tr45-48
       gchemTendency(i,j,k,bi,bj,44+iv)=gchemTendency(i,j,k,bi,bj,44+iv)
     &          +  Fsink(0+iv,k)/86400. ! tendency mmol/m3/d -> mmol/m3/s
             enddo
             do iv=1,5 ! Diatoms: Tr11-Tr15
       gchemTendency(i,j,k,bi,bj,10+iv)=gchemTendency(i,j,k,bi,bj,10+iv)
     &          +  Fsink(4+iv,k)/86400.! tendency mmol/m3/d ->mmol/m3/s
             enddo
             do iv=1,4 ! Flagellates: Tr16-Tr19
       gchemTendency(i,j,k,bi,bj,15+iv)=gchemTendency(i,j,k,bi,bj,15+iv)
     &          +  Fsink(9+iv,k)/86400.! tendency mmol/m3/d ->mmol/m3/s
             enddo
             do iv=1,4 ! Picophytoplankton: Tr20-Tr23    
       gchemTendency(i,j,k,bi,bj,20+iv)=gchemTendency(i,j,k,bi,bj,20+iv)
     &          +  Fsink(13+iv,k)/86400.! tendency mmol/m3/d ->mmol/m3/s
             enddo
             do iv=1,4 ! Dinoflagellates: Tr24-Tr27
       gchemTendency(i,j,k,bi,bj,23+iv)=gchemTendency(i,j,k,bi,bj,23+iv)
     &          +  Fsink(17+iv,k)/86400. ! tendency mmol/m3/d ->mmol/m3/s
             enddo
           endif
         enddo ! enddo k
         ENDDO   ! enddo j cycle
         ENDDO   ! enddo i cycle
#endif /* USE_SINK */

C  ==== add not only atmospheric deposition but also the contribution by 
C  ==== rivers for N1p and N3n to the tendency of the surface cell. 
C  ==== N1p_surfF and N3n_surfF read in the files are in mmol/m2/s. 
C  ==== Introduction of N5s_surfF, O3c_surfF and O3h_surfF contribution
C  ==== by rivers, where N5s_surfF and O3h_read are in mmol/m2/s,
C  ==== O3c_surfF read is in mgC/m2/s.
C  ==== These values are here divided by the effective height of the
C  ==== surface cell. 
      IF ( BFMcoupler_N1pSurfForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      gchemTendency(i,j,1,bi,bj,ppN1p)= gchemTendency(i,j,1,bi,bj,ppN1p)
     &              +N1p_surfF(i,j,bi,bj)
c     & *1/(drF(1)*hFacC(i,j,1,bi,bj)+EtaN(i,j,bi,bj))  ! tendency is mmol/m3/s
     &            *recip_drF(1)*recip_hFacC(i,j,1,bi,bj) !tendency is mmol/m3/s  
      ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_N3nSurfForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
      DO j=jmin,jmax
        DO i=imin,imax
      gchemTendency(i,j,1,bi,bj,ppN3n)= gchemTendency(i,j,1,bi,bj,ppN3n)
     &              +N3n_surfF(i,j,bi,bj) 
     &            *recip_drF(1)*recip_hFacC(i,j,1,bi,bj) !tendency is mmol/m3/s  
         ENDDO   
         ENDDO  
       ENDIF
       IF ( BFMcoupler_N5sSurfForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      gchemTendency(i,j,1,bi,bj,ppN5s)= gchemTendency(i,j,1,bi,bj,ppN5s)
     &              +N5s_surfF(i,j,bi,bj) 
     &    *recip_drF(1)*recip_hFacC(i,j,1,bi,bj) !tendency is mmol/m3/s          
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_O3cSurfForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      gchemTendency(i,j,1,bi,bj,ppO3c)= gchemTendency(i,j,1,bi,bj,ppO3c)
     &              +O3c_surfF(i,j,bi,bj) 
     &  *recip_drF(1)*recip_hFacC(i,j,1,bi,bj) !tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_O3hSurfForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      gchemTendency(i,j,1,bi,bj,ppO3h)= gchemTendency(i,j,1,bi,bj,ppO3h)
     &              +O3h_surfF(i,j,bi,bj) 
     &    *recip_drF(1)*recip_hFacC(i,j,1,bi,bj) !tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_N1pBotForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      kB=kBot(i,j)
      gchemTendency(i,j,kB,bi,bj,ppN1p)=
     &    gchemTendency(i,j,kB,bi,bj,ppN1p)+N1p_botF(i,j,bi,bj)
     &    *recip_drF(kB)*recip_hFacC(i,j,kB,bi,bj)!tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_N3nBotForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      kB=kBot(i,j)
      gchemTendency(i,j,kB,bi,bj,ppN3n)=
     &    gchemTendency(i,j,kB,bi,bj,ppN3n)+N3n_botF(i,j,bi,bj)
     &    *recip_drF(kB)*recip_hFacC(i,j,kB,bi,bj)!tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_N4nBotForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      kB=kBot(i,j)
      gchemTendency(i,j,kB,bi,bj,ppN4n)=
     &    gchemTendency(i,j,kB,bi,bj,ppN4n)+N4n_botF(i,j,bi,bj)
     &    *recip_drF(kB)*recip_hFacC(i,j,kB,bi,bj)!tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_O2oBotForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      kB=kBot(i,j)
      gchemTendency(i,j,kB,bi,bj,ppO2o)=
     &    gchemTendency(i,j,kB,bi,bj,ppO2o)+O2o_botF(i,j,bi,bj)
     &    *recip_drF(kB)*recip_hFacC(i,j,kB,bi,bj)!tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_O3cBotForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      kB=kBot(i,j)
      gchemTendency(i,j,kB,bi,bj,ppO3c)=
     &    gchemTendency(i,j,kB,bi,bj,ppO3c)+O3c_botF(i,j,bi,bj)
     &    *recip_drF(kB)*recip_hFacC(i,j,kB,bi,bj)!tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF
      IF ( BFMcoupler_O3hBotForcFile.NE. ' '  ) THEN
c        DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
       DO j=jmin,jmax
        DO i=imin,imax
      kB=kBot(i,j)
      gchemTendency(i,j,kB,bi,bj,ppO3h)=
     &    gchemTendency(i,j,kB,bi,bj,ppO3h)+O3h_botF(i,j,bi,bj)
     &    *recip_drF(kB)*recip_hFacC(i,j,kB,bi,bj)!tendency is mmol/m3/s  
        ENDDO
        ENDDO
      ENDIF


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C ==== fill the diagnostic variables 
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
#include "BFMcoupler_VARDIAG_fill_diags.h"
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

c      WRITE(msgBuf,'(A)')
c     &   ' prova messaggio'
c      CALL PRINT_MESSAGE( msgBuf, stdUnit,
c     I                    SQUEEZE_RIGHT, myThid )

      RETURN

      CONTAINS
c====== FUNCTIONS ============
       LOGICAL FUNCTION isnan2(A)
        implicit none
        real*8 A
        if ( A.eq.A ) then
             isnan2 = .FALSE.
         else
             isnan2 = .TRUE.
        end if
        RETURN

        END FUNCTION isnan2
        
        integer function bottom(i, j, bi, bj, Nr)
            
            implicit none
            
            integer, intent(in) :: i, j, bi, bj, Nr
            integer :: z, z1
            
            z = 1
c           TO DO: dangerous condition
            do z1 = 1,Nr 
               if (ceiling(HfacC(i, j, z1, bi, bj)) .eq. 1) then
                  z = z + 1
               endif
            enddo
            
            bottom = z - 1

        end function bottom
c=====================================
      END

#endif  /*ALLOW_BFMCOUPLER*/
#endif  /*ALLOW_PTRACERS*/


C============================================================================
