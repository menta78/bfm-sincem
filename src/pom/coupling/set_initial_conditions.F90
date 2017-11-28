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
!**                  bfm-community.eu                        **
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
! !ROUTINE: set_initial_conditions
!
!
! !INTERFACE
!
      subroutine set_initial_conditions
!
!DESCRIPTION
!
! This routine assigns initial conditons of biogeochemical variables in POM
! The I.C. definition for the LFG's and the organic CFF's is based on the carbon
! content. The Nitrogen, phosphorus and silicon content is obtained by applying 
! constant ratio with respect to carbon.
!
!********************************************************************************
!********************************************************************************
!**                                                                            **
!** IN THIS DEFAULT VERSION THE BFM STATE VAR'S ARE ASSIGNED THE VALUES SET IN **
!** THE NAMELISTS "bfm_General", and   BFM_init_nml_ben BUT OBVIOUSLY THIS     **
!** ROUTINE COULD/SHOULD BE TAILORED ACCCORDING TO THE IMPLEMENTATION          **
!** CHARACTERISTIC AND TO THE DATA AVAILABILITY.                               **
!** IF THE USER HAS A DATA FILE TO BE READ, THE READING STATEMENT              **
!** COULD BE PUT HERE IN SUBSTITUTION OF THE HARDWIRED I.C.'s                  **
!**                                                                            **
!** The initial conditions that are being specified using this default version **
!** of the subroutine are those used in:                                       **
!**                                                                            **
!** Mussap G., Zavatarelli M., Pinardi N., Celio M. (2016)                     **
!** A management oriented 1-D ecosystem model: Implementation in the Gulf of   **
!** Trieste. Regional Studies in Marine Science, 6, 109-123.                   **
!** doi http://dx.doi.org/10.1016/j.rsma2016.03.015                            **
!**                                                                            **
!** Mussap G., Zavatarelli M. (2017)                                           **
!** A numerical study of the benthic pelagic                                   **
!** coupling in a shallow shelf sea (Gulf of Trieste). Regional Studies in     **
!** Marine Science, 9, 24-34.                                                  **
!** doi http://dx.doi.org/10.1016/j.rsma2016.11.002                            **
!**                                                                            **
!** Mussap G., Zavatarelli M., Pinardi N. (2017)                               **
!** Linking coastal ocean modelling to environmental management: an ensemble   **
!** approach. Ocean dynamics, submitted                                        **
!**                                                                            **
!********************************************************************************
!********************************************************************************
!
!*******************************************************************************
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem, ONLY: ZERO,NML_OPEN,NML_READ,NMLUNIT,error_msg_prn
      use Mem 
      use POM
      use Service 
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----LOCAL SCALARS-----
!

      real(RLEN) :: d1cc, d1cn, d1cp, d1cs, d1ci
!
!     -----CONVERSION FACTORS-----
!
      real(RLEN),parameter ::     &
!
!                N/C RATIO
!
                 p_nRc=0.0126,    &
!
!                P/C RATIO
!
                 p_pRc=0.7862e-3, &
!
!                Si/C RATIO
!
                 p_sRc=0.0118,    &
!
!                Chl-a/C RATIO
!
                 p_iRc=1./25.
!
!
!     -----INITIAL VALUES (PELAGIC)FROM THE BFM_general NAMELIST-----
!

     real(RLEN) :: O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, B1n0, B1p0, &
      P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, P3n0, P3p0, &
      P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, Z4p0, Z5c0, &
      Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, R6c0, R6n0, &
      R6p0, R6s0, O3c0, O3h0
!
     namelist /bfm_init_nml/ O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, &
      B1n0, B1p0, P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, &
      P3n0, P3p0, P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, &
      Z4p0, Z5c0, Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, &
      R6c0, R6n0, R6p0, R6s0, O3c0, O3h0
!
!     -----INITIAL VALUES (BENTHIC) FROM THE BFM_init_nml_ben NAMELIST-----
!
#ifdef INCLUDE_BEN
!
     real(RLEN) :: G3c0, G3h0, G13c0, G13h0, G23c0, G23h0, Y1c0, Y1n0, Y1p0, &
      Y2c0, Y2n0, Y2p0, Y3c0, Y3n0, Y3p0, Y4c0, Y4n0, Y4p0, Y5c0, Y5n0, Y5p0, &
      Q6c0, Q6n0, Q6p0, Q6s0, Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, H1c0, H1n0, &
      H1p0, H2c0, H2n0, H2p0, K1p0, K11p0, K21p0, K4n0, K14n0, K24n0, K6r0, &
      K16r0, K26r0, K3n0, K5s0, G2o0, G4n0, D1m0, D2m0, D6m0, D7m0, D8m0, D9m0
!
     namelist /bfm_init_nml_ben/ G3c0, G3h0, G13c0, G13h0, G23c0, G23h0, Y1c0, &
      Y1n0, Y1p0, Y2c0, Y2n0, Y2p0, Y3c0, Y3n0, Y3p0, Y4c0, Y4n0, Y4p0, Y5c0, &
      Y5n0, Y5p0, Q6c0, Q6n0, Q6p0, Q6s0, Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, &
      H1c0, H1n0, H1p0, H2c0, H2n0, H2p0, K1p0, K11p0, K21p0, K4n0, K14n0, K24n0, &
      K6r0, K16r0, K26r0, K3n0, K5s0, G2o0, G4n0, D1m0, D2m0, D6m0, D7m0, D8m0, &
      D9m0
!
#endif
!
!    -----OPEN, READ AND CLOSE NAMELISTS-----
!
     open(NMLUNIT,file='BFM_General.nml',status='old',action='read',err=100)
     read(NMLUNIT,nml=bfm_init_nml,err=101)
     close(NMLUNIT)
!
#ifdef INCLUDE_BEN
!
     open(NMLUNIT,file='BFM_General.nml',status='old',action='read',err=102)
     read(NMLUNIT,nml=bfm_init_nml_ben,err=102)
     close(NMLUNIT)
!
#endif
!
!     ******************************************************
!     ******************************************************
!     **                                                  **
!     ** START INITIALIZATION OF THE PELAGIC STATE VAR'S  **
!     **                                                  **
!     ******************************************************
!     ******************************************************
!
!     -----INITIALIZE LIGHT PROFILE-----
!
!
        eir(:) = ZERO
!
!     -----INITIALIZE OXYGEN-----
!
        O2o(:) = O2o0
!
!     -----INITIALIZE CARBON DIOXYDE-----
!
      O3c(:)=O3c0
!
!     -----INITIALIZE DISSOLVED NUTRIENTS-----
!
        n1p(:) = N1p0
        n3n(:) = N3n0
        n4n(:) = N4n0
        n5s(:) = N5s0
        O4n(:) = O4n0
        N6r(:) = N6r0
!
!     -----INITIALIZE PARTICULATE ORGANIC MATTER-----
!
        R6c(:) = R6c0 
        R6n(:) = r6c(:)*p_nRc
        R6p(:) = r6c(:)*p_pRc
        R6s(:) = r6c(:)*p_sRc
!
!     -----INITIALIZE DISSOLVED ORGANIC MATTER-----
!
         R1c(:) = R1c0
         R1n(:) = R1c(:)*p_nRc*0.5
         R1p(:) = R1c(:)*p_pRc*0.5
         R2c(:) = R2c0
         R3c(:) = R3c0
!
!     -----INITIALIZE PHYTOPLANKTON (DIATOMS)-----
!
        d1cc=P1c0
        d1cn=d1cc*p_nRc
        d1cp=d1cc*p_pRc
        d1cs=d1cc*p_sRc
        d1ci=d1cc*p_iRc
      if (CalcPhytoPlankton(1)) then
        P1c(:) = d1cc
        P1n(:) = d1cn
        P1p(:) = d1cp
        P1s(:) = d1cs
        P1l(:) = d1ci
      else
        P1c(:) = ZERO
        P1n(:) = ZERO
        P1p(:) = ZERO
        P1s(:) = ZERO
        P1l(:) = ZERO
      end if
!
!     -----INITIALIZE PHYTOPLANKTON (NANOFLAGELLATES)-----
!
        d1cc = P2c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
      if (CalcPhytoPlankton(2)) then
        P2c(:) = d1cc
        P2n(:) = d1cn
        P2p(:) = d1cp
        P2l(:) = d1ci
      else
        P2c(:) = ZERO
        P2n(:) = ZERO
        P2p(:) = ZERO
        P2l(:) = ZERO
      end if
!
!     -----INITIALIZE PHYTOPLANKTON (PICOPHYTOPLANKTON)-----
!
        d1cc = P3c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
      if (CalcPhytoPlankton(3)) then
        P3c(:) = d1cc
        P3n(:) = d1cn
        P3p(:) = d1cp
        P3l(:) = d1ci
      else
        P3c(:) = ZERO
        P3n(:) = ZERO
        P3p(:) = ZERO
        P3l(:) = ZERO
      end if
!
!     -----INITIALIZE PHYTOPLANKTON (DINOFLAGELLATES)-----
!
        d1cc = P4c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
      if (CalcPhytoPlankton(4)) then
        P4c(:) = d1cc
        P4n(:) = d1cn
        P4p(:) = d1cp
        P4l(:) = d1ci
      else
        P4c(:) = ZERO
        P4n(:) = ZERO
        P4p(:) = ZERO 
        P4l(:) = ZERO
      end if
!
!     -----INITIALIZE MESOZOOPLANKTON (CARNIVOROUS)---_
!
      if (CalcMesoZooPlankton(1)) then
        Z3c(:) = Z3c0
        Z3n(:) = Z3c(:)*p_nRc
        Z3p(:) = Z3c(:)*p_pRc
      else
        Z3c(:) = ZERO
        Z3n(:) = ZERO
        Z3p(:) = ZERO
      end if
!
!     -----INITIALIZE MESOZOOPLANKTON (OMNIVOROUS)-----
!
      if (CalcMesoZooPlankton(2)) then
        Z4c(:) = Z4c0
        Z4n(:) = Z4c(:)*p_nRc
        Z4p(:) = Z4c(:)*p_pRc
      else
        Z4c(:) = ZERO 
        Z4n(:) = ZERO
        Z4p(:) = ZERO
       end if

!
!     -----INITIALIZE MICROZOOPLANKTON (GENERIC)-----
!
      if (CalcMicroZooPlankton(1)) then
         Z5c(:) = Z5c0
         Z5n(:) = Z5c(:)*p_nRc
         Z5p(:) = Z5c(:)*p_pRc
      else
         Z5c(:) = ZERO
         Z5n(:) = ZERO
         Z5p(:) = ZERO
      end if
!
!     -----INITIALIZE MICROZOOPLANKTON (HETEROTROPHIC FLAGELLATES)-----
!
      if (CalcMicroZooPlankton(2)) then
         Z6c(:) = Z6c0
         Z6n(:) = Z6c(:)*p_nRc
         Z6p(:) = Z6c(:)*p_pRc
      else
         Z6c(:) = ZERO
         Z6n(:) = ZERO
         Z6p(:) = ZERO
      end if
!
!     -----INITIALIZE BACTERIA-----
!
      if (CalcPelBacteria(1)) then 
        B1c(:) = B1c0
        B1n(:) = B1c(:)*p_nRc
        B1p(:) = B1c(:)*p_pRc
      else
        B1c(:) = ZERO
        B1n(:) = ZERO
        B1p(:) = ZERO
      end if
!
!     ******************************************************
!     ******************************************************
!     **                                                  **
!     ** START INITIALIZATION OF THE BENTHIC STATE VAR'S  **
!     **                                                  **
!     ******************************************************
!     ******************************************************
!
#ifdef INCLUDE_BEN
!
!     -----INITIALIZE MACROBENTHOS-----
!
      if (CalcBenOrganisms(1)) then
        Y1c(:) = Y1c0
      else
        Y1c(:) = ZERO
      end if
!
!     -----INITIALIZE DETRITIVORES-----
!

      if (CalcBenOrganisms(2)) then
        Y2c(:) = Y2c0
      else
        Y2c(:) = ZERO
      end if
!
!     -----INITIALIZE FILTER FEEDERS-----
!

      if (CalcBenOrganisms(3)) then
        Y3c(:) = Y3c0
      else
        Y3c(:) = ZERO
      end if
!
!     -----INITIALIZE MEIOBENTHOS-----
!
      if (CalcBenOrganisms(4)) then
        Y4c(:) = Y4c0
      else
        Y4c(:) = ZERO
      end if
!
!     -----INITIALIZE PREDATORS-----
!

      if (CalcBenOrganisms(5)) then
        Y5c(:) = Y5c0
      else
        Y5c(:) = ZERO 
      end if

!
!     -----INITIALIZE BACTERIA-----
!
      if (CalcBenBacteria(1)) then
        H1c(:) = H1c0
      else
        H1c(:) = ZERO
      end if

      if (CalcBenBacteria(2)) then
        H2c(:) = H2c0
      else  
        H2c(:) = ZERO 
      end if
!
!     -----INITIALIZE NUTRIENTS-----
!
       K5s(:)  = 20.75
       K6r(:)  = K6r0
       k4n(:)  = K4n0
       K14n(:) = K14n0
       k24n(:) = K24n0
       k1p(:)  = K1p0
       K11p(:) = K11p0
       K21p(:) = K21p0
       K3n(:)  = K3n0
!
!     -----INITIALIZE DISSOLVED ORGANIC MATTER-----
!
       Q1c(:)  = Q1c0
       Q11c(:) = Q11c0
!
!     -----INITIALIZE PARTICULATE ORGANIC MATTER-----
!
       Q6c(:) = Q6c0
       Q6n(:) = Q6n0
       Q6p(:) = Q6p0
       Q6s(:) = Q6s0
!
!     -----INITIALIZE OXYGEN-----
!
       G2o(:) = G2o0
!
!      -----INITIALIZE CARBON DIOXYDE-----
!
!
       G3c(:) = G3c0
!
!      -----INITIALIZE MOLECULAR NITROGEN-----
!
       G4n(:) = 37.0
!
!      -----INITIALISE BENTHIC LAYERS THICKNESS-----
!
        D1m(:) = D1m0
        D2m(:) = D2m0
        D6m(:) = D6m0
        D7m(:) = D7m0
        D8m(:) = D8m0
        D9m(:) = D9m0
!
#endif
!
      return
!
!    -----PRINT IF ERRORS WITH NAMELISTS READING-----
!
 100 call error_msg_prn(NML_OPEN,"InitParam.f90","BFM_General.nml")
 101 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters")
 102 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters_ben")
!
      end subroutine set_initial_conditions
!
