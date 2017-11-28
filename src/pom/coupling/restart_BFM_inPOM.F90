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
   subroutine restart_BFM_inPOM
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
   end subroutine restart_BFM_inPOM
!
