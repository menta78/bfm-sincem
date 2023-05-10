!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BFMcoupler_readparms
!
! DESCRIPTION
!    Initialize Memory of BFM and read BFM package parameters
!     from BFM namelist
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
      SUBROUTINE BFMcoupler_READPARMS( myThid )
      use api_bfm

      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "BFMcoupler_VARS.h"
       

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid    :: My Thread Id. number
      INTEGER myThid
CEOP
#ifdef ALLOW_BFMCOUPLER
C     === Local variables ===
C     msgBuf    :: Informational/error message buffer
C     iUnit     :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /BFMcoupler_FORCING/
     &          BFMcoupler_atmospFile,
     &          BFMcoupler_PCO2File
#ifdef READ_xESP
     &          ,BFMcoupler_xESPFile
#endif

     &          ,BFMcoupler_sparFile

     &          ,BFMcoupler_N1pSurfForcFile,BFMcoupler_N3nSurfForcFile,
     &          BFMcoupler_forcingPeriod, BFMcoupler_forcingCycle,
     &          BFMcoupler_N5sSurfForcFile,BFMcoupler_O3cSurfForcFile,
     &          BFMcoupler_O3hSurfForcFile,BFMcoupler_N1pBotForcFile, 
     &          BFMcoupler_N3nBotForcFile,BFMcoupler_N4nBotForcFile,
     &          BFMcoupler_O2oBotForcFile,BFMcoupler_O3cBotForcFile,
     &          BFMcoupler_O3hBotForcFile
#ifdef USE_SINK
     &          ,vsed
#endif
#ifdef USE_SHADE
     &          ,p_espP1,p_espP2,p_espP3,p_espP4,p_espWATER  
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _BEGIN_MASTER(myThid)

       BFMcoupler_atmospFile= ' '
       BFMcoupler_PCO2File  = ' '
       BFMcoupler_xESPFile  = ' '
       BFMcoupler_sparFile  = ' '
       BFMcoupler_N1pSurfForcFile= ' '
       BFMcoupler_N3nSurfForcFile= ' '
       BFMcoupler_N5sSurfForcFile= ' '
       BFMcoupler_O3cSurfForcFile= ' '
       BFMcoupler_O3hSurfForcFile= ' '
       BFMcoupler_N1pBotForcFile= ' '
       BFMcoupler_N3nBotForcFile= ' '
       BFMcoupler_N4nBotForcFile= ' '
       BFMcoupler_O2oBotForcFile= ' '
       BFMcoupler_O3cBotForcFile= ' '
       BFMcoupler_O3hBotForcFile= ' '
       BFMcoupler_forcingPeriod = externForcingPeriod
       BFMcoupler_forcingCycle  = externForcingCycle
#ifdef USE_SINK
       vsed=3.0
#endif
#ifdef USE_SHADE
       p_espP1=10.0E-3 ! [m2 mgCHL-1] use default if not present in NAMELIST
       p_espP2=8.0E-3 ! [m2 mgCHL-1] use default if not present in NAMELIST
       p_espP3=8.0E-3 ! [m2 mgCHL-1] use default if not present in NAMELIST
       p_espP4=8.0E-3 ! [m2 mgCHL-1] use default if not present in NAMELIST
       p_espWATER=0.035 ! [m-1] use default if not present in NAMELIST
#endif

      WRITE(msgBuf,'(A)')'BFMcoupler_READPARMS:check on light use: '
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )
#ifdef USE_SHADE
      WRITE(msgBuf,'(A)')'USE_SHADE defined: atten from shelfshading'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )
#ifdef READ_xESP
      WRITE(msgBuf,'(A)')'xESPfile used as background values for water'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )
#else
      WRITE(msgBuf,'(A)')'background from p_espWATER parameter'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )
#endif
#else
#ifdef READ_xESP
      WRITE(msgBuf,'(A)')'atten read from xESPfile  '
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )
#else
      WRITE(msgBuf,'(A)')'USE_SHADE and READ_xESP not defined::  error'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )
      STOP  'ABNORMAL END: BFMcoupler_readparams'
#endif 
#endif
      WRITE(msgBuf,'(A)')'BFMcoupler_READPARMS: opening data.BFMcoupler'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      CALL OPEN_COPY_DATA_FILE('data.BFMcoupler','BFMcoupler_READPARMS',
     O                          iUnit, myThid )
C-    forcing filenames and parameters
      READ(UNIT=iUnit,NML=BFMcoupler_FORCING)

      WRITE(msgBuf,'(A)')
     &   ' BFMcoupler_READPARMS: finished reading data.BFMcoupler'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

C--   Close the open data file
      CLOSE(iUnit)


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C Inizialize Libraries of BFM model 

      WRITE(msgBuf,'(A)') 'BFMcoupler: Inizializing BFM libraries model'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )


      call BFM1D_NO_BOXES(Nr,1,1,Nr,1)
      parallel_rank=MyProcId
      call Init_bfm()
      call BFM1D_INIT_IO_CHANNELS()


      call Initialize()
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|



#endif /* ALLOW_BFMCOUPLER */

      RETURN
      END




