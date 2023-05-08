!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BFMcoupler_FIELDS_LOAD
!
! DESCRIPTION
!    Reads if necessary, AtmosPCO2, AtmosP, spar, xESP, N1p_surfF,
!    N3n_surfF, N5s_surfF, O3c_surfF, O3h_surfF
!    N1p_botF, N3n_botF, N4n_botF,O2o_botF, O3c_botF, O3h_botF
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!


#include "BFMcoupler_OPTIONS.h"


C !INTERFACE: ==========================================================
      SUBROUTINE BFMcoupler_FIELDS_LOAD (
     I           myIter,myTime,myThid)


      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "BFMcoupler_VARS.h"
#include "BFMcoupler_LOAD.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
      INTEGER myIter
      _RL myTime
      INTEGER myThid

#ifdef ALLOW_BFMCOUPLER

c !LOCAL VARIABLES: ===================================================
      INTEGER bi, bj, i, j
      INTEGER intimeP, intime0, intime1
      _RL aWght,bWght

#ifdef READ_xESP  
      CHARACTER*(MAX_LEN_MBUF) msgBuf_xeps
#endif
#ifdef READ_PAR  
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#endif
CEOP

      IF (  BFMcoupler_forcingCycle.gt.0. _d 0 ) THEN

C--   Now calculate whether it is time to update the forcing arrays
       CALL GET_PERIODIC_INTERVAL(
     O               intimeP, intime0, intime1, bWght, aWght,
     I               BFMcoupler_forcingCycle, BFMcoupler_forcingPeriod,
     I               deltaTclock, myTime, myThid )

       bi = myBxLo(myThid)
       bj = myByLo(myThid)
c#ifdef ALLOW_DEBUG
c       IF ( debugLevel.GE.debLevB ) THEN
c        _BEGIN_MASTER(myThid)
c        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
c     &   ' BFMcoupler_FIELDS_LOAD,', myIter,
c     &' : iP,iLd,i0,i1=',intimeP,BFMcoupler_ldRec(bi,bj),intime0,intime1
c     &  ,' ; Wght=', bWght, aWght
c        _END_MASTER(myThid)
c       ENDIF
c#endif /* ALLOW_DEBUG */

#ifdef ALLOW_AUTODIFF_TAMC
C-    assuming that we call S/R BFMcoupler_FIELDS_LOAD at each time-step and
C     with increasing time, this will catch when we need to load new records;
C     But with Adjoint run, this is not always the case => might end-up using
C     the wrong time-records
       IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
#else /* ALLOW_AUTODIFF_TAMC */
C-    Make no assumption on sequence of calls to BFMcoupler_FIELDS_LOAD ;
C     This is the correct formulation (works in Adjoint run).
C     Unfortunatly, produces many recomputations <== not used until it is fixed
       IF ( intime1.NE.BFMcoupler_ldRec(bi,bj) ) THEN
#endif /* ALLOW_AUTODIFF_TAMC */

C--   If the above condition is met then we need to read in
C     data for the period ahead and the period behind myTime.
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' BFMcoupler_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, BFMcoupler_ldRec(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        _BARRIER

        IF ( BFMcoupler_atmosPFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( BFMcoupler_atmospFile,atmosP0,intime0,
     &        myIter,myThid )
         CALL READ_REC_XY_RS( BFMcoupler_atmospFile,atmosP1,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( BFMcoupler_PCO2File .NE. ' '  ) THEN
        CALL READ_REC_XY_RS(BFMcoupler_PCO2File,AtmosPCO20,intime0,
     &        myIter,myThid )
        CALL READ_REC_XY_RS(BFMcoupler_PCO2File,AtmosPCO21,intime1,
     &        myIter,myThid )
        ENDIF
        IF ( BFMcoupler_N1pSurfForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_N1pSurfForcFile,N1p_surfF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_N1pSurfForcFile,N1p_surfF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_N3nSurfForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_N3nSurfForcFile,N3n_surfF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_N3nSurfForcFile,N3n_surfF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_N5sSurfForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_N5sSurfForcFile,N5s_surfF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_N5sSurfForcFile,N5s_surfF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_O3cSurfForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_O3cSurfForcFile,O3c_surfF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_O3cSurfForcFile,O3c_surfF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_O3hSurfForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_O3hSurfForcFile,O3h_surfF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_O3hSurfForcFile,O3h_surfF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_N1pBotForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_N1pBotForcFile,N1p_botF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_N1pBotForcFile,N1p_botF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_N3nBotForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_N3nBotForcFile,N3n_botF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_N3nBotForcFile,N3n_botF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_N4nBotForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_N4nBotForcFile,N4n_botF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_N4nBotForcFile,N4n_botF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_O2oBotForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_O2oBotForcFile,O2o_botF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_O2oBotForcFile,O2o_botF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_O3cBotForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_O3cBotForcFile,O3c_botF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_O3cBotForcFile,O3c_botF1,
     &       intime1,myIter,myThid )
        ENDIF
        IF ( BFMcoupler_O3hBotForcFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_O3hBotForcFile,O3h_botF0,
     &       intime0,myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_O3hBotForcFile,O3h_botF1,
     &       intime1,myIter,myThid )
        ENDIF
#ifdef READ_xESP
        IF ( BFMcoupler_xESPFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_xESPFile,xESP0,intime0,
     &       myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_xESPFile,xESP1,intime1,
     &       myIter,myThid )
         ENDIF
#endif
#ifdef READ_PAR
        IF ( BFMcoupler_sparFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS(BFMcoupler_sparFile,spar0,intime0,
     &       myIter,myThid )
         CALL READ_REC_XY_RS(BFMcoupler_sparFile,spar1,intime1,
     &       myIter,myThid )
         ENDIF
#endif

C--   fill-in overlap after loading temp arrays:
        _EXCH_XY_RS(AtmosP0, myThid )
        _EXCH_XY_RS(AtmosP1, myThid )
        _EXCH_XY_RS(AtmosPCO20, myThid )
        _EXCH_XY_RS(AtmosPCO21, myThid )
        _EXCH_XY_RS(N1p_surfF0, myThid )
        _EXCH_XY_RS(N1p_surfF1, myThid )
        _EXCH_XY_RS(N3n_surfF0, myThid )
        _EXCH_XY_RS(N3n_surfF1, myThid )
        _EXCH_XY_RS(N5s_surfF0, myThid )
        _EXCH_XY_RS(N5s_surfF1, myThid )
        _EXCH_XY_RS(O3c_surfF0, myThid )
        _EXCH_XY_RS(O3c_surfF1, myThid )
        _EXCH_XY_RS(O3h_surfF0, myThid )
        _EXCH_XY_RS(O3h_surfF1, myThid )
        _EXCH_XY_RS(N1p_botF0, myThid )
        _EXCH_XY_RS(N1p_botF1, myThid )
        _EXCH_XY_RS(N3n_botF0, myThid )
        _EXCH_XY_RS(N3n_botF1, myThid )
        _EXCH_XY_RS(N4n_botF0, myThid )
        _EXCH_XY_RS(N4n_botF1, myThid )
        _EXCH_XY_RS(O2o_botF0, myThid )
        _EXCH_XY_RS(O2o_botF1, myThid )
        _EXCH_XY_RS(O3c_botF0, myThid )
        _EXCH_XY_RS(O3c_botF1, myThid )
        _EXCH_XY_RS(O3h_botF0, myThid )
        _EXCH_XY_RS(O3h_botF1, myThid )
#ifdef READ_xESP
        _EXCH_XY_RS(xESP0, myThid )
        _EXCH_XY_RS(xESP1, myThid )
#endif
#ifdef READ_PAR
        _EXCH_XY_RS(spar0, myThid )
        _EXCH_XY_RS(spar1, myThid )
#endif
C-    save newly loaded time-record
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           BFMcoupler_ldRec(bi,bj) = intime1
         ENDDO
        ENDDO
C-------------------------------------------
C-     end if-bloc (time to load new fields)
       ENDIF
C------------------------------------------
    
      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
  
         IF (BFMcoupler_AtmosPFile .NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             AtmosP(i,j,bi,bj) = bWght*atmosp0(i,j,bi,bj)
     &                         + aWght*atmosp1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF (BFMcoupler_PCO2File .NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             AtmosPCO2(i,j,bi,bj) = bWght*atmospco20(i,j,bi,bj)
     &                         + aWght*atmospco21(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
#ifdef READ_xESP
         IF ( BFMcoupler_xESPFile .NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             xESP(i,j,bi,bj) = bWght*xESP0(i,j,bi,bj)
     &                       + aWght*xESP1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
#endif
#ifdef READ_PAR
         IF ( BFMcoupler_sparFile .NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             spar(i,j,bi,bj) = bWght*spar0(i,j,bi,bj)
     &                      + aWght*spar1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
#endif
         IF ( BFMcoupler_N1pSurfForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             N1p_surfF(i,j,bi,bj) = bWght*N1p_surfF0(i,j,bi,bj)
     &                          + aWght*N1p_surfF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_N3nSurfForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             N3n_surfF(i,j,bi,bj) = bWght*N3n_surfF0(i,j,bi,bj)
     &                          + aWght*N3n_surfF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_N5sSurfForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             N5s_surfF(i,j,bi,bj) = bWght*N5s_surfF0(i,j,bi,bj)
     &                          + aWght*N5s_surfF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_O3cSurfForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             O3c_surfF(i,j,bi,bj) = bWght*O3c_surfF0(i,j,bi,bj)
     &                          + aWght*O3c_surfF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_O3hSurfForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             O3h_surfF(i,j,bi,bj) = bWght*O3h_surfF0(i,j,bi,bj)
     &                          + aWght*O3h_surfF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_N1pBotForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             N1p_botF(i,j,bi,bj) = bWght*N1p_botF0(i,j,bi,bj)
     &                          + aWght*N1p_botF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_N3nBotForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             N3n_botF(i,j,bi,bj) = bWght*N3n_botF0(i,j,bi,bj)
     &                          + aWght*N3n_botF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_N4nBotForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             N4n_botF(i,j,bi,bj) = bWght*N4n_botF0(i,j,bi,bj)
     &                          + aWght*N4n_botF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_O2oBotForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             O2o_botF(i,j,bi,bj) = bWght*O2o_botF0(i,j,bi,bj)
     &                          + aWght*O2o_botF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_O3cBotForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             O3c_botF(i,j,bi,bj) = bWght*O3c_botF0(i,j,bi,bj)
     &                          + aWght*O3c_botF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
         IF ( BFMcoupler_O3hBotForcFile.NE. ' '  ) THEN
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             O3h_botF(i,j,bi,bj) = bWght*O3h_botF0(i,j,bi,bj)
     &                          + aWght*O3h_botF1(i,j,bi,bj)
            ENDDO
           ENDDO
         ENDIF
        ENDDO
       ENDDO

C endif for BFMcoupler_forcingCycle
      ENDIF


#endif /* ALLOW_BFMCOUPLER */
      RETURN
      END
