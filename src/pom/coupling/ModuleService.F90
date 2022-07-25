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
!     -----SURFACE NUTRIENTS AND RUNOFF (ONLY USED IF NUTBC_MODE == 1)-----
!
      real(RLEN)                        :: PO4SURF,NO3SURF,NH4SURF,SIO4SURF,DISSURF
!
!     ----- Profiles of currens speed (ONLY USED IF NUTBC_MODE == 1)----
      real(RLEN)                        :: CURRENTS_SPEED(KB-1) = 1 ! by default using a constant profile
!
!     -----SUSPENDED INORGANIC MATTER PROFILE-----
!
      real(RLEN),public,dimension(KB-1) :: ISM
!
!     ------FREQUENCY OF OUTPUT AVERAGING (IN HOURS)-----
!
      integer(ilong)                    :: savef
!
!     -----THESE ARE THE PATHWAYS FOR THE IC, RESTART AND FORCING FILES (READ TROUGH NML)-----
!
     character(200)                     :: wind_input,      &
                                           ism_input,       &
                                           Sal_input,       &
                                           Temp_input,      &
                                           Sprofile_input,  &
                                           Tprofile_input,  &
                                           Cprofile_input,  & ! file with horizontal currents speed profile, needed if NUTSBC_MODE == 1
                                           heat_input,      &
                                           surfNut_input,   &
                                           read_restart
!
 end module Service
!
!EOC
