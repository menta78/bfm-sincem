#include "INCLUDE.h"
#include "cppdefs.h"
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
! !ROUTINE:  restart_BFM_inPOM
!
!DESCRIPTION    
!
!  This routine writes the file needed for the BFM restart 
!
!****************************************************************
!
! !INTERFACE
!
subroutine save_restart
!
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use netcdf_bfm, only: save_rst_bfm, ncid_rst
   use constants,  only:SEC_PER_DAY
   use pom, ONLY: time 
   use global_mem, ONLY: RLEN
   use Mem
   use api_bfm
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
  IMPLICIT NONE
!
!  -----TIME IN SECONDS-----
!
   real(RLEN)               :: localtime 
!
! -----TIME ELLAPSED (IN SECONDS) SINCE THE BEGINNING OF THE SIMULATION----
!
   localtime = time*SEC_PER_DAY   
!
!  -----WRITE RESTART-----
!
   call save_rst_bfm(localtime)       
!
! -----CLOSE OUTPUT AND RESTART FILES------
!
   call close_ncdf(ncid_rst)
   call close_ncdf(ncid_bfm)

 write (6,*) 'NETCDF RESTART WRITTEN, TIME--> ', time
!
   return
end subroutine save_restart
!



subroutine load_restart
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
      use CPL_VARIABLES, only: wind_input,     &
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
end subroutine load_restart
