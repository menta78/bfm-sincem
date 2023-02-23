!
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
! !ROUTINE: Opendat
!
!DESCRIPTION    
!
! This subroutine opens the files containing the forcing data. 
! the path to the data file comes from the reading of the
! Namelist "pom_input".
!
! The forcing data files are unformatted and written in direct access.
! User can obviously change the the data "format" and "access" mode,
! consistently with the data they have at hand.
!
! Any change in the open arguments, requires, however, coherent reading
! procedure in subroutine FORCING_MANAGER.
!
!*****************************************************************************
!
! !INTERFACE
!
    subroutine opendat
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
    use global_mem,ONLY:RLEN,LOGUNIT,error_msg_prn,NML_OPEN,NML_READ
!
    use CPL_VARIABLES, only:  wind_input,     & ! Wind Stress
                        Sal_input,      & ! Salinity Initial Conditions
                        Temp_input,     & ! Temperature Initial Conditions
                        Sprofile_input, & ! Time varying (monthly) salinity profiles
                        Tprofile_input, & ! Time varying (monthly) temperature profiles
                        Kprofile_input, & ! Time varying (monthly) KH profiles
                        Wprofile_input, & ! Time varying (monthly) vertical velocity profiles
                        heat_input,     & ! Time varying (monthly) surface heat flux
                        surfNut_input,  & ! Time varying (monthly) surface nutrient
                        ism_input,      & ! Inorganic suspended matter Initial Conditions
                        BENT_INT_SCHEME,& ! scheme for benthic time integration: 0: leapfrog, 1: euler-forward. Default: 0
                        USE_KH_EXT,     & ! true if KH is loaded from an external source
                        KH_FACTOR,      & ! scaling factor for KH_EXT
                        SEDI_FACTOR,    & ! adjustment factor of sedimentation in vertical advection
                        USE_W_PROFILE,  &
                        NUTSBC_MODE,    & ! NUTRIENT SURFACE BOUNDARY CONDITIONS: 0 (default): concentrations. 1: fluxes
                        L_PO4,          &
                        L_NO3,          &
                        L_SIO4,         &
                        L_O2,           &
                        L_X,            &
                        ASURF_PO4,      &
                        ASURF_NO3,      &
                        ASURF_NH4,      &
                        ASURF_SIO4,     &
                        SWR_FILE_STEP,  & ! if 1 shortwave radiation is loaded hourly if 0 monthly
                        read_restart      ! Model restart file
!
     use pom, ONLY :    TB, SB, ZZ, H
!
     use forcing, ONLY: WSU1,WSV1,                 &
                        ISM1,                      &
                        SCLIM1,                    &
                        TCLIM1,                    &
                        SWRAD1,                    &
                        NO3_1,NH4_1,PO4_1, SIO4_1, & 
                        KH_1,                      &
                        WMN1, WVR1,                &
                        QCORR1                     ! NO MORE IN USE!!!!!
!
!    -----IMPLICIT TYPING IS NEVER ALLOWED----
!
     IMPLICIT NONE
!
!    ----- INTEGER SCALAR-----
!
     INTEGER :: RLENGTH
!
!    -----NAMELIST READING UNIT-----
!
     integer,parameter  :: namlst=10
!
!    -----OPEN, READ AND CLOSE NAMELIST-----
!
     namelist /pom_input/ wind_input,     &
                          ism_input,      &
                          Sal_input,      &
                          Temp_input,     &
                          Sprofile_input, &
                          Tprofile_input, &
                          Kprofile_input, &
                          Wprofile_input, &
                          heat_input,     &
                          surfNut_input,  &
                          NUTSBC_MODE,    &
                          KH_FACTOR,      &
                          SEDI_FACTOR,    &
                          L_PO4,          &
                          L_NO3,          &
                          L_SIO4,         &
                          L_O2,           &
                          L_X,            &
                          ASURF_PO4,      &
                          ASURF_NO3,      &
                          ASURF_NH4,      &
                          ASURF_SIO4,     &
                          SWR_FILE_STEP,  &
                          read_restart
!
     open(namlst,file='pom_bfm_settings.nml',status='old',action='read',err=100)
     read(namlst,nml=pom_input, err=102)
     close(namlst)
!
     write(6,*) 'Loaded namelist:'
     write(6,pom_input)
!
!
!    -----OPEN WIND STRESS FILE-----
!
     inquire(IOLENGTH=rlength) WSU1,WSV1
     open(11,file=wind_input, form='unformatted', access='direct', recl=rlength)
     write(6,*) 'open 11 done'
!
!    -----OPEN INORGANIC SUSPENDED MATTER PROFILES FILE-----
!
     inquire(IOLENGTH=rlength) ISM1(1)
     open(19, file=ism_input, form='unformatted', access='direct', recl=rlength)
     write(6,*) 'open 19 done'
!
!    -----OPEN CLIMATOLOGICAL SALINITY PROFILES FILE-----
!
     inquire(IOLENGTH=rlength) SCLIM1(1)
     open(20, file=Sal_input, form='unformatted', access='direct', recl=rlength)
     write(6,*) 'open 20 done'
!
!    -----OPEN CLIMATOLOGICAL TEMPERATURE PROFILES FILE-----
!
     inquire(IOLENGTH=rlength) TCLIM1(1)
     open(15, file=Temp_input, form='unformatted', access='direct', recl=rlength)
     write(6,*) 'open 15 done'

!
!    -----OPEN HEAT FLUX FILE-----
!
!    ****************************************
!    ****************************************
!    **                                    **
!    ** N.B. QCORR1 IS NO LONGER USED      **
!    **                                    **
!    ****************************************
!    ****************************************
!
     SELECT CASE (SWR_FILE_STEP)
        CASE (1) ! the file contains hourly data
           inquire(IOLENGTH=rlength) SWRAD1
           open(21, file=heat_input, form='unformatted', access='direct', recl=rlength)
           write(6,*) 'open 21 done'
        CASE DEFAULT ! the file contains monthly means
           inquire(IOLENGTH=rlength) SWRAD1,SWRAD1,QCORR1
           open(21, file=heat_input, form='unformatted', access='direct', recl=rlength)
           write(6,*) 'open 21 done'
     END SELECT
!
!    -----OPEN NUTRIENTS FILE-----
!
     inquire(IOLENGTH=rlength) NO3_1,NH4_1,PO4_1,SIO4_1
     open(18, file=surfNut_input, form='unformatted',access='direct',recl=rlength)
     write(6,*) 'open 18 done'
!
!    -----OPEN KH PROFILES, IF PROVIDED---
!
     inquire(FILE=Kprofile_input, EXIST=USE_KH_EXT)
     IF (USE_KH_EXT) THEN
         inquire(IOLENGTH=rlength) KH_1(1)
         write(6,*) 'KH profile file exists, opening it :',Kprofile_input
         open(34, file=Kprofile_input, form='unformatted',access='direct',recl=rlength)
         write(6,*) 'open 34 done'
     ELSE
         write(6,*) 'KH profile file ',Kprofile_input,' does not exist, setting USE_KH_EXT=.FALSE.'
     END IF
!
!    -----OPEN W PROFILE, IF PROVIDED----
!
     inquire(FILE=Wprofile_input, EXIST=USE_W_PROFILE)
     IF (USE_W_PROFILE) THEN
         inquire(IOLENGTH=rlength) WMN1(1),WVR1(1)
         write(6,*) 'W profile file exists, opening it :',Wprofile_input
         open(35, file=Wprofile_input, form='unformatted',access='direct',recl=rlength)
         write(6,*) 'open 35 done'
     ELSE
         write(6,*) 'W profile file ',Wprofile_input,' does not exist, setting USE_W_PROFILE=.FALSE.'
     END IF

     write(6,*) 'open units done'
!
#ifdef SAVEFORCING
     open(400,file="surface_forcing.txt")
     write(400,'(a8,8a18)') 'MONTH','WSU','WSV','SWRAD','NO3','NH4','PO4','SIO4'
     open(401,file="suspmat_forcing.txt")
     write(401,'(a8,40f18.8)') 'MONTH',ABS(zz*h)
     open(402,file="salvert_forcing.txt")
     write(402,'(a8,40f18.8)') 'MONTH',ABS(zz*h)
     open(403,file="temvert_forcing.txt")
     write(403,'(a8,40f18.8)') 'MONTH',ABS(zz*h)
#endif
!
     return
!
!    -----PRINT IF PROBLEM WITH NAMELIST OPENING-----
!
100   call error_msg_prn(NML_OPEN,"opendat.F90","pom_bfm_settings.nml")
!
!     -----PRINT IF PROBLEM WITH NAMELIST READING-----
!
102   call error_msg_prn(NML_READ,"opendat.F90","pom_input in pom_bfm_settings.nml")

      end subroutine opendat
!
