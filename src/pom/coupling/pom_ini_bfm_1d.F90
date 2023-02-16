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
! !ROUTINE: pom_ini_bfm
!
! !INTERFACE:
!
  subroutine pom_ini_bfm_1d
!
!  DESCRIPTION:
!  BFM setup and initialisation when run in coupled mode with POM
!
!********************************************************************
!
!  -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
   use mem
!
   use POM, ONLY: DTI,     &
                  KB,      &
                  H,       &
                  DZ,      &
                  ihotst,  &
                  ALAT,    &
                  ZZ,      &
                  DZZ,     &
                  NC_OUT_STARTTIME
!
   use global_mem, ONLY:ZERO, &
                        RLEN
!
   use api_bfm
!
   use init_var_bfm_local
!
   use CPL_VARIABLES, ONLY: savef
!
   use netcdf_bfm, ONLY: init_netcdf_bfm,     &
                         init_save_bfm,       &
                         init_netcdf_rst_bfm, &
                         read_rst_bfm

   use time, ONLY: timestep
!
!  -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
   IMPLICIT NONE
!
!  ---- LOOP COUNTERS-----
!
   integer    :: k,i
!  
!  ----NAMELISTS READING READING UNITS-----
!
   integer,parameter    :: namlst=10,unit=11
!  
!  -----DUMMY LONGITUDE VALUE (NOT USED IN 1D MODE)-----
!
   real(RLEN),parameter :: ALON=10.
!
!  -----MODEL TIME STEP TRANSFERED TO BFM-----
!
   timestep=dti
!
!  ----OUTPUT AVERAGING AND WRITING FREQUENCY TRANSFERED TO BFM-----
!
   out_delta = savef
!
!  ----SET THE BFM ARRAYS DIMENSION-----
!
!    *************************************************************************
!    *************************************************************************
!    **                                                                     **
!    ** SINCE THIS IS A 1D IMPLEMENTATION THE SIZE OF THE ARRAYS IS SHRUNK  **
!    ** TO A 1D VECTOR  ALONG the "VERTICAL" DIMENSION (1, 1, KB-1)         **
!    **                                                                     **
!    *************************************************************************
!    *************************************************************************
!            
     NO_BOXES=KB-1
     NO_BOXES_X=1
     NO_BOXES_Y=1
     NO_BOXES_Z=NO_BOXES
     NO_BOXES_XY=1
!
!   -----ALLOCATE AND SET LOGICAL MASKS FOR:-----
!
!           --WHOLE WATER COLUMN---
   allocate(SEAmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!           ---BOTTOM---
   allocate(BOTmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!           ---SURFACE--
   allocate(SRFmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!
   SEAmask = .TRUE.
   BOTmask = .TRUE.
   SRFmask = .TRUE.
!
!  -----ALLOCATE CPL_VARIABLES ARRAY AND INITIALISE (ZEROING)-----
!
   allocate(ZEROS(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z)) 
   ZEROS = ZERO
!
!  -----TOTAL NUMBER OF VARIABLES TO BE COMPUTED-----
!
!
!  -----PELAGIC+BENTHIC SETUP-----
!
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES_BEN 
   NO_BOXES_Z_BEN  = 1
   NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
   NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN
!
!  -----COMPRESSED COORDINATES FOR NETCDF OUTPUT-----
!
   lon_len = NO_BOXES_X 
   lat_len = NO_BOXES_Y
!
!  -----ALLOCATE LAND-SEA MASK----
!
   allocate(ocepoint(NO_BOXES))
!
!  -----INITIALISE LAND (0) SEA (1) MASK. IN 1D ONLY SEA....-----
!
   ocepoint = 1 
!
!  *********************************************************
!  *********************************************************
!  **                                                     **
!  ** ALLOCATE AND LOAD THE ARRAY CONTAINING THE INDICES  **
!  ** OF THE PELAGIC GRIDPOINTS HAVING A WATER-SEDIMENT   **
!  ** INTERFACE.                                          **
!  **                                                     **
!  ** IN 1D MODE ONLY THE GRID POINT AT K=KB-1            **
!  ** EXCHANGES WITH THE SEDIMENTS                        **
!  **                                                     **
!  *********************************************************
!  *********************************************************
!
   allocate(BOTindices(NO_BOXES_XY)); BOTindices = KB-1
!
!  *********************************************************
!  *********************************************************
!  **                                                     **
!  ** ALLOCATE AND LOAD THE ARRAY CONTAINING THE INDICES  **
!  ** OF THE PELAGIC GRIDPOINTS HAVING A WATER-ATMOSPHERE **
!  ** INTERFACE.                                          **
!  **                                                     **
!  ** IN 1D MODE ONLY THE GRID POINT AT K=1               **
!  ** EXCHANGES WITH THE ATMOSPHERE                       **
!  **                                                     **
!  *********************************************************
!  *********************************************************
!
   allocate(SRFindices(NO_BOXES_XY)); SRFindices = 1
!
!  ***********************************************************************
!  ***********************************************************************
!  **                                                                   **
!  ** READ THE BFM_GENERAL NAMELIST, ALLOCATE ARRAYS WITH ATTRIBUTES OF **
!  ** STATE VARIABLES AND WRITE SETUP INFO IN RUN LOG.                  **
!  **                                                                   **
!  ***********************************************************************
!  ***********************************************************************
!
   call init_bfm
!
!  -----INITIALISE STATE VARIABLE NAMES AND DIAGNOSTICS (FOR NETCDF OUTPUT)-----
!
   call set_var_info_bfm
!
!  ****************************************************************************
!  ****************************************************************************
!  **                                                                        **
!  ** ALLOCATE BFM ARRAYS AND PROVIDE A VERY GENERAL INITIALISATION.         **
!  **                                                                        **
!  ** N.B.: AT THIS STAGE THE INITIALISATION IS VERY GENERIC (UNIFORM IN     **
!  ** SPACE CARBON CONTENT VALUES. THE CORRESPONDING NITROGEN, PHOSPHORUS    **
!  ** AND SILICON CONTENT (IF ANY) IS DEFINED ON THE BASIS OF A "REDFIELD"   **
!  ** RATIO ASSUMPTION.                                                      **
!  ** A MORE PRECISE INITIALISATION CAN BE PROVIDED TROUGH THE USE OF THE    **
!  ** SUBROUTINE SET_INITIAL CONDITION BELOW                                 **
!  **                                                                        **
!  ****************************************************************************
!  ****************************************************************************
!!
   call init_var_bfm(bio_setup)

! Initialize internal constitutents quota of functional groups
   call ini_organic_quotas()
!
!   ******************************************************************
!   ******************************************************************
!   **                                                              **
!   ** FIX THE FLAG INDICATING THE "COLD"/"HOT" START               **
!   **                                                              **
!   ** N.B.: IN THE BFM-POM1D SETUP THE LEADING  "COLD"/'HOT"       **
!   **       START FLAG IS IHOTST (DEFINED IN THE params_POMBFM     **
!   **       NAMELIST.                                              **
!   **       THE bfm_init VALUE READ IN BFM_GENERAL IS OVERWRITTEN  **
!   **       WITH THE IHOTST VALUE                                  **
!   **                                                              **
!   ******************************************************************
!   ******************************************************************
!
    bfm_init=ihotst
!
!    -----SET THE THICKNESS OF THE WATER COLUMN LAYERS-----
!
    do k = 1 , NO_BOXES_Z
       Depth(k) = abs(dz(k)*h)
    end do
!
!  -----INITIALISE THE OUTPUT TIME AVERAGING PROCEDURE-----
!
   call calcmean_bfm(INIT)
!
!  -----INITIALISE NETCDF OUTPUT-----
!
   call init_netcdf_bfm(title=out_title,                      &
                        start_time=TRIM(NC_OUT_STARTTIME),              &
                        expinfo="BFM_POM",                    &
                        time_unit=0,                          &
                        lat=alat,                             &
                        lon=alon,                             &
                        z=zz*h,                                 &
                        dz=dzz*h,                               &
                        oceanpoint=ocepoint,                  &
                        surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
                        bottompoint=(/(i,i=1,NO_BOXES_XY)/))
!
!  *******************************************************************
!  *******************************************************************
!  **                                                               **
!  ** ALLOCATE AND INITIALISE  THE ARRAY TO BE LOADED WITH STATE    **
!  ** VARIABLES (PELAGIC AND BENTHIC COMPUTED AT t=t-dt.            **
!  ** THIS IS NEEDED FOR THE LEAPGROG NUMERICAL INTEGRATION SCHEME  **
!  ** USED BY THE BFM-POM1D SYSTEM                                  **
!  **                                                               **
!  *******************************************************************
!  *******************************************************************
!
   allocate(D3STATEB(NO_BOXES,NO_D3_BOX_STATES)) !PELAGIC
!
!  -----ZEROING-----
!
   D3STATEB = ZERO
!
   allocate(D2STATEB_BEN(NO_BOXES_XY,NO_D2_BOX_STATES_BEN)) ! BENTHIC
!
!  -----ZEROING-----
!
   D2STATEB_BEN = ZERO
!
!  -----SAVE INITIAL CONDITIONS IN OUTPUT FILE-----
!
   call init_save_bfm
!
! -----IF "COLD" START (bfm-init=0) STATE VAR'S VALUES AT t AND t-DT ARE THE SAME-----
!
!
     D3STATEB = D3STATE !PELAGIC
!
     D2STATEB_BEN = D2STATE_BEN ! BENTHIC
!
!
!  -----IF "HOT" START (Bfm_init = 1) READ BFM RESTART FILE -----
!
if (bfm_init == 1) call read_rst_bfm(in_rst_fname)
!
!      -----INITIALISE BFM RESTART FILE-----
!
       call init_netcdf_rst_bfm(out_rst_fname,                            &
                                start_time=TRIM(NC_OUT_STARTTIME),              &
                                time_unit=0,                          &
                                lat=alat,                             &
                                lon=alon,                             &
                                z=zz*h,                                 &
                                dz=dzz*h,                               &
                                oceanpoint=ocepoint,                  &
                                surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
                                bottompoint=(/(i,i=1,NO_BOXES_XY)/))
!
   return
!
 end subroutine pom_ini_bfm_1d
!EOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

