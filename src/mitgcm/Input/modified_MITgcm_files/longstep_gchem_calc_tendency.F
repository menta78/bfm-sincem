C $Header: longstep_gchem_calc_tendency.F,v 1.5 2015/04/01 02:52:57 jmc Exp $


#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: LONGSTEP_GCHEM_CALC_TENDENCY

C !INTERFACE: ==========================================================
      SUBROUTINE LONGSTEP_GCHEM_CALC_TENDENCY(
     I     myTime, myIter, myThid )

C !DESCRIPTION:
C version of GCHEM_CALC_TENDENCY called by LONGSTEP package
C
C In the case of GCHEM_SEPARATE_FORCING not defined,
C this subroutine computes the tracer tendencies due to a
C bio-geogchemistry or ecosystem model and stores them on an array
C gchemTendency, that will be incorporated into regular timestepping in
C in ptracers_intergrate.F
C The current example uses the CFC package, but his is meant to
C be replaced by anything that the user provides.

C  This routine has been modified to use the BFM library through
C  the BFMCOUPLER package in Jan 2017 for to by G.Cossarini (gcossari@inogs.it)


C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"
#include "GCHEM_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time
C  myIter               :: current timestep
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_GCHEM
# ifndef GCHEM_SEPARATE_FORCING
# ifdef ALLOW_LONGSTEP
C !LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  iTracer              :: ptracer number
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER i,j
      INTEGER iTracer
      INTEGER bi,bj,iMin,iMax,jMin,jMax,k

C     gchemTendency is re-initialized here
      DO iTracer = 1, PTRACERS_numInUse
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            gchemTendency(i,j,k,bi,bj,iTracer) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C
C     Here is the place for code to compute bio-geochemical
C     tendency terms (sometimes referred to as source-minus-sink
C     terms). The tendencies are stored on gchemTendency, as show
C     in the CFC-example.
C
C     loop over tiles
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C     define horizontal loop ranges
        iMin = 1
        iMax = sNx
        jMin = 1
        jMax = sNy
c       DO iTracer = 1, PTRACER_numInUse
c        DO k = 1, Nr
c         DO j = 1-OLy, sNy+OLy
c          DO i = 1-OLx, sNx+OLx
c           gchemTendency(i,j,k,bi,bj,iTracer) = your specific model
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO

C------------------------
C chemical forcing      |
C------------------------
#ifdef ALLOW_CFC
        IF ( useCFC ) THEN
         CALL CFC11_FORCING(
     I                       pTracer     (1-OLx,1-OLy,1,bi,bj,1),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,1),
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       myTime, myIter, myThid )
         CALL CFC12_FORCING(
     I                       pTracer     (1-OLx,1-OLy,1,bi,bj,2),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,2),
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       myTime, myIter, myThid )
        ENDIF
#endif /* ALLOW_CFC */


C     end of tile-loops
       ENDDO
      ENDDO
C------------------------
C BFM coupler           |
C------------------------
#ifdef ALLOW_BFMCOUPLER
      IF ( useBFMcoupler ) THEN
C$taf loop = parallel
        DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)

        jMin=1
        jMax=sNy
        iMin=1
        iMax=sNx
C BFMcoupler operates on bi,bj part only, but needs to get full arrays i
C because of last index (iPtr)
         CALL BFMcoupler_calc_tendency(bi,bj,imin,imax,jmin,jmax,
     &                            myIter,myTime,myThid)
        ENDDO
        ENDDO
       ENDIF
#endif /* ALLOW_BFMCOUPLER */

#endif /* ALLOW_LONGSTEP */
#endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
