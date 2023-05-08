!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BFMcoupler_ini_forcing
!
! DESCRIPTION
!    Initializes forcing fields,
!    taking in account surface and bottom contributions
!
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
      SUBROUTINE BFMcoupler_ini_forcing( myThid )


      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "BFMcoupler_VARS.h"
#include "BFMcoupler_LOAD.h"  

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid


#ifdef ALLOW_BFMCOUPLER

c !LOCAL VARIABLES: ===================================================
       INTEGER bi,bj,i,j,k
#if (defined (READ_PAR)  || defined (USE_QSW))
       CHARACTER*(MAX_LEN_MBUF) msgBuf   !!! reading of KPAR
#endif

      CALL LEF_ZERO( atmosP0,myThid )
      CALL LEF_ZERO( atmosP1,myThid )
      CALL LEF_ZERO( atmosPCO20,myThid )
      CALL LEF_ZERO( atmosPCO21,myThid )
#ifdef READ_xESP
      CALL LEF_ZERO( xESP0,myThid ) ! KPAR
      CALL LEF_ZERO( xESP1,myThid ) 
#endif
#ifdef READ_PAR
      CALL LEF_ZERO( spar0,myThid ) ! surface PAR
      CALL LEF_ZERO( spar1,myThid ) 
#endif
      CALL LEF_ZERO( N1p_surfF0,myThid ) ! N1p surface forcing 
      CALL LEF_ZERO( N1p_surfF1,myThid ) 
      CALL LEF_ZERO( N3n_surfF0,myThid ) ! N3n surface forcing
      CALL LEF_ZERO( N3n_surfF1,myThid ) 
      CALL LEF_ZERO( N5s_surfF0,myThid ) ! N5s surface forcing
      CALL LEF_ZERO( N5s_surfF1,myThid )
      CALL LEF_ZERO( O3c_surfF0,myThid ) ! O3c surface forcing
      CALL LEF_ZERO( O3c_surfF1,myThid )
      CALL LEF_ZERO( O3h_surfF0,myThid ) ! O3h surface forcing
      CALL LEF_ZERO( O3h_surfF1,myThid )
      CALL LEF_ZERO( N1p_botF0,myThid )  ! N1p bottom forcing 
      CALL LEF_ZERO( N1p_botF1,myThid )
      CALL LEF_ZERO( N3n_botF0,myThid )  ! N3n bottom forcing 
      CALL LEF_ZERO( N3n_botF1,myThid )
      CALL LEF_ZERO( N4n_botF0,myThid )  ! N4n bottom forcing 
      CALL LEF_ZERO( N4n_botF1,myThid )
      CALL LEF_ZERO( O2o_botF0,myThid )  ! O2o bottom forcing 
      CALL LEF_ZERO( O2o_botF1,myThid )
      CALL LEF_ZERO( O3c_botF0,myThid )  ! O3c bottom forcing 
      CALL LEF_ZERO( O3c_botF1,myThid )
      CALL LEF_ZERO( O3h_botF0,myThid )  ! O3h bottom forcing 
      CALL LEF_ZERO( O3h_botF1,myThid )
#ifdef READ_PAR
#ifdef USE_QSW
        WRITE(msgBuf,'(2A)')
     &       ' BFMcoupler_INI_FORCING: You can not use READ_PAR ',
     &       ' and USE_QSW'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BFMcoupler_INI_FORCING'
#endif
#endif

c set reasonable values where/when needed
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-Oly,sNy+Oly
          DO i=1-Olx,sNx+Olx
             AtmosP(i,j,bi,bj)   = 1. _d 0*maskC(i,j,1,bi,bj)
             AtmosPCO2(i,j,bi,bj) = 390. _d 0*maskC(i,j,1,bi,bj) ! CO2 atm conc. [ppm]
             xESP(i,j,bi,bj)     = 0.04*maskC(i,j,1,bi,bj) ! kpar [m-1]
             spar(i,j,bi,bj)    = 250. _d 0*maskC(i,j,1,bi,bj) ! surface par [uEin/m2/s]
             N1p_surfF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj) ! N1p surface forcing
             N3n_surfF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj) ! N3n surface forcing
             N5s_surfF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj) ! N5s surface forcing
             O3c_surfF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj) ! O3c surface forcing 
             O3h_surfF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj) ! O3h surface forcing
             N1p_botF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj) ! N1p bottom forcing 
             N3n_botF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj)
             N4n_botF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj)
             O2o_botF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj)
             O3c_botF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj)
             O3h_botF(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj)
             DO k = 1,Nr
                PH(i,j,k,bi,bj)=8.0
             ENDDO

          ENDDO
         ENDDO
        ENDDO
       ENDDO


#endif /* ALLOW_BFMCOUPLER */
      RETURN
      END
