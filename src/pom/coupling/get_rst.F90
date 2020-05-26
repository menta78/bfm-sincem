!  ROUTINE: get_rst
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
! INTERFACE
!
      subroutine get_rst
!
! DESCRIPTION
!
! This subroutine opens and reads the ".fort" file containing the
! restart data for the POM component of the BFM-POM1D modeling system
!
! The input path to the restart file is specified in pom_input.nml
! The output is handled via module POM
!
!***********************************************************************
!

!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem, ONLY: error_msg_prn, NML_OPEN, NML_READ
!
      use Service, only: wind_input,     &
                         ism_input,      &
                         Sal_input,      &
                         Temp_input,     &
                         Sprofile_input, &
                         Tprofile_input, &
                         heat_input,     &
                         surfNut_input,  &
                         read_restart
!
      use POM, ONLY: time0,       &
                     u,ub,        &
                     v,vb,        &
                     t,tb,        &
                     s,sb,        &
                     q2,q2b,      &
                     q2l,q2lb,    &
                     kh,km,kq,    &
                     l,           &
                     wubot,wvbot, &
                     rho
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE 
!
!     -----NAMELIST READING UNIT-----
!
      integer,parameter  :: namlst=10
!
!     -----READ NAMELIST WITH RESTART FILE PATH-----
!
       namelist /pom_input/ wind_input,     &
                            ism_input,      &
                            Sal_input,      &
                            Temp_input,     &
                            Sprofile_input, &
                            Tprofile_input, &
                            heat_input,     &
                            surfNut_input,  &
                            read_restart
!
       open(namlst,file='pom_bfm_settings.nml',status='old',action='read',err=100)
       read(namlst,nml=pom_input, err=102)
       rewind(namlst)
       close(namlst)
!
!    -----OPEN AND READ THE RESTART FILE FOR POM-----
!
       open(70,file=read_restart,form='unformatted',status='old')
!
       read(70) time0,       &
                u,ub,        &
                v,vb,        &
                t,tb,        &
                s,sb,        &
                q2,q2b,      &
                q2l,q2lb,    &
                kh,km,kq,    &
                l,           &
                wubot,wvbot, &
                rho
!
      return
!
!     -----PRINT IF SOMETHING IS WRONG WITH NAMELIST OPENING-----
!
100   call error_msg_prn(NML_OPEN,"get_rst.F90","pom_bfm_settings.nml")
!
!     -----PRINT IF SOMETHING IS WRONG WITH NAMELIST READING-----
!
102   call error_msg_prn(NML_READ,"get_rst.F90","pom_input in pom_bfm_settings.nml")
!
      end subroutine get_rst
