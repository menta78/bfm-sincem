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
    use Service, only:  wind_input,     & ! Wind Stress
                        Sal_input,      & ! Salinity Initial Conditions
                        Temp_input,     & ! Temperature Initial Conditions
                        Sprofile_input, & ! Time varying (monthly) salinity profiles
                        Tprofile_input, & ! Time varying (monthly) temperature profiles
                        heat_input,     & ! Time varying (monthly) surface heat flux
                        surfNut_input,  & ! Time varying (monthly) surface nutrient
                        ism_input,      & ! Inorganic suspended matter Initial Conditions
                        read_restart      ! Model restart file
!
     use pom, ONLY :    TB, SB
!
     use forcing, ONLY: WSU1,WSV1,                 &
                        ISM1,                      &
                        SCLIM1,                    &
                        TCLIM1,                    &
                        SWRAD1, WTSURF1,           &
                        NO3_1,NH4_1,PO4_1, SIO4_1, &
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
                          heat_input,     &
                          surfNut_input,  &
                          read_restart
!
     open(namlst,file='pom_input.nml',status='old',action='read',err=100)
     read(namlst,nml=pom_input, err=102)
     close(namlst)
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
     inquire(IOLENGTH=rlength) SWRAD1,WTSURF1,QCORR1
     open(21, file=heat_input, form='unformatted', access='direct', recl=rlength)
     write(6,*) 'open 21 done'
!
!    -----OPEN INORGANIC SUSPENDED MATTER FILE-----
!
     inquire(IOLENGTH=rlength) NO3_1,NH4_1,PO4_1,SIO4_1
     open(18, file=surfNut_input, form='unformatted',access='direct',recl=rlength)
     write(6,*) 'open 18 done'

      write(6,*) 'open units done'
!
     return
!
!    -----PRINT IF PROBLEM WITH NAMELIST OPENING-----
!
100   call error_msg_prn(NML_OPEN,"opendat.F90","pom_input.nml")
!
!     -----PRINT IF PROBLEM WITH NAMELIST READING-----
!
102   call error_msg_prn(NML_READ,"opendat.F90","pom_input")

      end subroutine opendat
!
