!
! **************************************************************
! **************************************************************
! **                                                          **
! ** ONE-DIMENSIONAL BFM-POM  MODELING SYSTEM (BFM-POM1D)     **
! **                                                          **
! ** The modeling system originate from the direct on-line    **
! ** coupling of the 1D Version of the Princeton Ocean model  **
! ** "POM" and the Biological Flux Model "BFM".               **
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
! ** Subsequent maintance by Lorenzo Mentaschi.               **
! ** Thanks are due to Prof. George L. Mellor that allowed us **
! ** to modify use and distribute the one dimensional         **
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
! !ROUTINE: get_init_TS_IC
!
! !INTERFACE
!
     subroutine get_init_TS_IC
!
!DESCRIPTION
!
! This subroutine opens and reads files containing the T&S initial conditions
! Files are read in direct access mode
! The path to the T&S I.C. file specified in namelist pom_input.
!
!***********************************************************************************
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem, ONLY: error_msg_prn, NML_OPEN, NML_READ
!
      use CPL_VARIABLES, ONLY: Sprofile_input, &
                         Tprofile_input
!    
      use pom, ONLY: KB,T,TB,S,SB
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----LOOP COUNTER-----
!
      INTEGER :: K
!
!      -----RECORD LENGTH-----
!
      INTEGER :: RLENGTH
!
!     -----NAMELIST READING UNIT-----
!
      integer,parameter  :: namlst=10
!
!    -----OPEN FILE WITH SALINITY I.C.----
!
       inquire(IOLENGTH=rlength) SB(1)
       open(29,file=Sprofile_input,form='unformatted',access='direct',recl=rlength)
!
!
!    -----OPEN FILE WITH TEMPERATURE I.C.----
!
!
       inquire(IOLENGTH=rlength) TB(1)
       open(10,file=Tprofile_input,form='unformatted',access='direct',recl=rlength)
!
!
!    -----READ T&S INITIAL CONDITIONS-----
!
     DO K = 1,KB
!
           READ (29,REC=K) SB(K)
           READ (10,REC=K) TB(K)
!
     END DO
!
!    -----COLD START: T@(t)=T@(t-dt)-----
!
     T(:)=TB(:)
     S(:)=SB(:)
!
     return
!
!    -----PRINT IF PROBLEMS WITH NML OPENING-----
!
100   call error_msg_prn(NML_OPEN,"get_init_TS_IC.F90","problem opening pom_bfm_settings.nml")
!
!    -----PRINT IF PROBLEMS WITH NML READING-----
!
102   call error_msg_prn(NML_READ,"get_init_TS_IC.F90","pom_input in pom_bfm_settings.nml")
!
      end subroutine get_init_TS_IC
