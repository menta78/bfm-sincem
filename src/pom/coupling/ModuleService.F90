 !
! !MODULE: Service
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
! DEFINITION AND ALLOCATION OF PARAMETERS, SCALAR AND ARRAYS USED
! FOR THE BFM-POM COUPLING.
!
! !INTERFACE
  MODULE  Service
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use POM,ONLY: KB,ilong
!
      use global_mem, ONLY:RLEN
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
! 
!     -----SURFACE NUTRIENTS AND RUNOFF (ONLY USED IF NUTSBC_MODE == 1)-----
!
      real(RLEN)                        :: PO4SURF,NO3SURF,NH4SURF,SIO4SURF,DISSURF
!
!     -----USE_KH_EXT==.TRUE. if the vertical diffusion coefficient is loaded from an external file
!
      logical                           :: USE_KH_EXT = .FALSE.
!
!     -----vertical diffusion coefficient from an external source and scaling constat
!
      real(RLEN)                        :: KH_EXT(KB-1) = 0
      real(RLEN)                        :: KH_FACTOR = 1
!
!     ----- USE_O2_TNDC==.TRUE. if an O2 profile file was provided
!
      logical                           :: USE_O2_TNDC = .FALSE.
!
!     ----- Profiles of oxygen (ONLY USED IF NUTBC_MODE == 1)----
!
      real(RLEN)                        :: O2_TNDC(KB-1) ! tendency O2 profile. Used if NUTSBC_MODE == 1
!
!     ----- USE_W_PROFILE == .TRUE. if daily profiles of W are provided ----
!
      logical                           :: USE_W_PROFILE = .FALSE.
!
!     ----- vertical profiles of W -----
!
      real(RLEN)                        :: W_MEAN_PR(KB) ! vert. prof. of w mean
      real(RLEN)                        :: W_VRNC_PR(KB) ! vert. prof. of w variance
      real(RLEN)                        :: W_PROFILE(KB)
!
!     -----SUSPENDED INORGANIC MATTER PROFILE-----
!
      real(RLEN),public,dimension(KB-1) :: ISM
!
!     ----- relaxation coefficients for: PO4, NO3, SiO4, O2, anything else
!
      real(RLEN)                        :: L_PO4, L_NO3, L_SIO4, L_O2, L_X
!
!     -----NUTRIENT SURFACE BOUNDARY CONDITIONS MODE:
!     -----    0: surface flux is computed applying the relaxation time NRT (default)
!     -----    1: surface flux is computed applying the nutrient concentration to a river runoff
!
      INTEGER                  :: NUTSBC_MODE=0
!
!     ! short wave radiation stepping: if 1 shortwave radiation is loaded hourly if 0 monthly
!
      INTEGER                  :: SWR_FILE_STEP=0

!
!     ------FREQUENCY OF OUTPUT AVERAGING (IN HOURS)-----
!
      integer(ilong)                    :: savef

      integer(ilong)                    :: DAY_OF_SIMULATION = -999999
!
!     -----THESE ARE THE PATHWAYS FOR THE IC, RESTART AND FORCING FILES (READ TROUGH NML)-----
!
     character(200)                     :: wind_input,      &
                                           ism_input,       &
                                           Sal_input,       &
                                           Temp_input,      &
                                           Sprofile_input,  &
                                           Tprofile_input,  &
                                           Oprofile_input,  & ! file with O2 profile, used if NUTSBC_MODE == 1
                                           Kprofile_input,  & ! file with vertical diffusion coeff. profile, if available
                                           Wprofile_input,  & ! file with daily profile of W, if available
                                           heat_input,      &
                                           surfNut_input,   &
                                           read_restart
!
 end module Service
!
!EOC
