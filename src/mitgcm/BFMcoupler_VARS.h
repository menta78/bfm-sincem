! file BFMcoupler_VARS.h

C     *==========================================================*
C     | BFMcoupler_VARS.h
C     *==========================================================*

C--   COMMON /BFMcoupler_FILENAMES/
C  BFMcoupler_atmospFile    :: file name of atmospheric pressure
C  BFMcoupler_PCO2File      :: file name of atmospheric pCO2
C  BFMcoupler_xESPFile      :: file name of satellite light extintion factor
C  BFMcoupler_N1pSurfForcFile    :: file name of N1p atmospheric dry deposition and contribution 
C  BFMcoupler_N3nSurfForcFile    :: file name of N3n atmospheric dry deposition and contribution
C  BFMcoupler_sparFile      :: file name of surface photosynthetically available radiation (PAR)
C  BFMcoupler_forcingPeriod :: periodic forcing parameter specific for BFMcoupler (seconds)
C  BFMcoupler_forcingCycle  :: periodic forcing parameter specific for BFMcoupler (seconds)
C
C  BFMcoupler_N5sSurfForcFile    :: file name of N5s contribution
C  BFMcoupler_O3cSurfForcFile    :: file name of O3c contribution
C  BFMcoupler_O3hSurfForcFile    :: file name of O3h contribution
C
C  BFMcoupler_N1pBotForcFile    :: file name of N1p contribution
C  BFMcoupler_N3nBotForcFile    :: file name of N3n contribution
C  BFMcoupler_N4nBotForcFile    :: file name of N4n contribution
C  BFMcoupler_O2oBotForcFile    :: file name of O2o contribution
C  BFMcoupler_O3cBotForcFile    :: file name of O3c contribution
C  BFMcoupler_O3hBotForcFile    :: file name of O3h contribution

      COMMON /BFMcoupler_FILENAMES/
     &    BFMcoupler_atmospFile, BFMcoupler_PCO2File,
     &    BFMcoupler_xESPFile, BFMcoupler_N1pSurfForcFile, 
     &    BFMcoupler_N3nSurfForcFile,BFMcoupler_sparFile,
     &    BFMcoupler_forcingPeriod, BFMcoupler_forcingCycle,
     &    BFMcoupler_N5sSurfForcFile, BFMcoupler_O3cSurfForcFile,
     &    BFMcoupler_O3hSurfForcFile, BFMcoupler_N1pBotForcFile,
     &    BFMcoupler_N3nBotForcFile, BFMcoupler_N4nBotForcFile,
     &    BFMcoupler_O2oBotForcFile, BFMcoupler_O3cBotForcFile,
     &    BFMcoupler_O3hBotForcFile
     

      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_atmospFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_PCO2File
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_xESPFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_sparFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_N1pSurfForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_N3nSurfForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_N5sSurfForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_O3cSurfForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_O3hSurfForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_N1pBotForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_N3nBotForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_N4nBotForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_O2oBotForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_O3cBotForcFile
      CHARACTER*(MAX_LEN_FNAM) BFMcoupler_O3hBotForcFile
      _RL     BFMcoupler_forcingPeriod
      _RL     BFMcoupler_forcingCycle

C ==================================================
C include the list of BFM variables 
C and list of diagnostic of BFM
C The BFM_var_list.h must be picked up from BFM code
C ==================================================
#include "BFM_var_list.h"

       COMMON /CARBON_NEEDS/
     &              AtmosPCO2, AtmosP,AtmosWIND,PH
      _RL  AtmosPCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  AtmosP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  AtmosWIND(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  PH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy) 



       COMMON /LIGHT_NEEDS/
     &          spar,xESP,xpar
      _RL  spar(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  xESP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  xpar(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy) 

       COMMON /SEDIMENTATION_RATE/
     &      sediPI,vsed
      _RL sediPI(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,4) ! 4 phytoplankton
      _RL vsed ! sinking rate for detritus (constant)

       COMMON /NUT_DEPOSITION/
     &          N1p_surfF,N3n_surfF,N5s_surfF,O3c_surfF,
     &          O3h_surfF,N1p_botF,N3n_botF,N4n_botF,
     &          O2o_botF,O3c_botF,O3h_botF
      _RL  N1p_surfF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  N3n_surfF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  N5s_surfF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  O3c_surfF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  O3h_surfF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  N1p_botF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  N3n_botF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  N4n_botF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  O2o_botF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  O3c_botF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  O3h_botF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef USE_SHADE
       COMMON /LIGHT_SHADE_VARS/
     &      p_espP1,p_espP2,p_espP3,p_espP4,p_espWATER
      _RL p_espP1
      _RL p_espP2
      _RL p_espP3
      _RL p_espP4
      _RL p_espWATER
#endif


CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***



