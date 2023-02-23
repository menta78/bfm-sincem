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

! !ROUTINE: vert_integration
!
! !INTERFACE
!
  SUBROUTINE vert_integration
!
!
!DESCRIPTION
!
!  This subroutine computes the forward in time BFM state variables
!  according to the Source splitting coupling technique, as described in:
! 
!  Butenschoen M., Zavatarelli M., Vichi M. (2012)
!  Sensitivity of a marine coupled physical biogeochemical model to time 
!  resolution,integration scheme and time splitting method.
!  Ocean Modeling, 52/53, 36-53.
!
!********************************************************************
!
!    -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
     use global_mem, ONLY:RLEN,ZERO
!
     use CPL_VARIABLES
!
     use api_bfm, ONLY:D3STATEB
!
     use constants, ONLY: SEC_PER_DAY
!
     use mem_Param, ONLY: p_small, AssignPelBenFluxesInBFMFlag
!
     use mem_PelSinkSet
!
     use Mem, ONLY:D3STATE,D3SOURCE,NO_D3_BOX_STATES,                     &
                   ppO2o,                                                 &
                   ppO3c,                                                 &
                   ppN1p,ppN3n,ppN4n,ppN5s,                               &
                   ppR6c,ppR6n,ppR6p,ppR6s,                               &
                   ppP1c,ppP1n,ppP1p,ppP1s,ppP1l,                         &
                   ppP2c,ppP2n,ppP2p,ppP2l,                               &
                   ppP3c,ppP3n,ppP3p,ppP3l,                               &
                   ppP4c,ppP4n,ppP4p,ppP4l,                               &
                   sediR6,sediPPY,                                        &
                   iiP1,iiP2,iiP3,iiP4,                                   &
                   n1p,n3n,n4n,n5s,                                       &
                   jsurO2o,jsurO3c,                                       &
                   Depth,                                                 &
                   iiPhytoPlankton,                                       &
                   jbotN1p,jbotN4n,jbotN3n,                               &
                   jbotO2o,jbotO3c,jbotN5s
!
     use POM, ONLY:SMOTH,KB,H,DTI,DZ,DZR,NRT,KH,NBCBFM,UMOLBFM,NTP,GRAV,RHO

!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
!
  IMPLICIT NONE
!
!-----COUNTER & FLAGS-----
!
 integer                  :: k,m,n,nbc
!
!-----BFM STATE VAR. @ time t,  t+DTI, T-DTI RESPECTIVELY-----
!
 real(RLEN)               :: fbio(KB), ffbio(KB), fbbio(KB), advtnd(KB) 
!
!-----SURFACE FLUX STORAGE-----
!
 real(RLEN)               :: surflux
!
!-----BOTTOM FLUX STORAGE-----
!
 real(RLEN)               :: botflux       
!
!-----RELAXATION VELOCITY FOR NUTRIENTS SURFACE FLUX-----
!
 real(RLEN)               :: vrelax
!
!-----SEDIMENTATION VELOCITY-----
!
 real(RLEN)               :: sink(KB), W(KB)
!
!-----TWICE THE TIME STEP-----
!
 real(RLEN)               :: dti2  
!
!-----lateral relaxation and flux
!
 real(RLEN)               :: LAMBDA, LFLUX(KB)
!
!
    dti2 = DTI*2.

!
!-----COMPUTING RELAXATION VELOCITY FOR SURFACE NUTRIENTS FLUX-----
!
    vrelax=NRT/SEC_PER_DAY
!
!-----LOOP OVER BFM STATE VAR'S-----
!

  IF (USE_KH_EXT) THEN
      KH(:KB-1) = KH_EXT*KH_FACTOR
  END IF

  do m = 1 , NO_D3_BOX_STATES
!
!         -----LOAD BFM STATE VAR.-----
!
          fbio(:KB-1)  = D3STATE(:,m)
          fbbio(:KB-1) = D3STATEB(:,m)
!
!         -----ZEROING SINKING VELOCITY-----
! 
          sink(:)=ZERO
!
!         -----COSMETIC: VALUES AT K=KB ARE NOT USED-----
!
          fbio(kb) =fbio(kb-1)
          fbbio(kb)=fbbio(kb-1)
!
!         -----SURFACE AND BOTTOM FLUXES-----
!
          select case ( m )
!
!                -----OXYGEN-----
!
                 case (ppO2o)
!
                       surflux = -(jsurO2o(1)/SEC_PER_DAY)
!
                       botflux = jbotO2o(1)/SEC_PER_DAY
!
!                -----CARBON DIOXYDE-----
!
                 case (ppO3c)
!
                       surflux= -(jsurO3c(1)/SEC_PER_DAY)
!
                       botflux = jbotO3c(1)/SEC_PER_DAY
!
!                -----PHOSPHATE----
!
                 case (ppN1p)
!
                      SELECT CASE (NUTSBC_MODE)
                           CASE (1)
                              surflux = -PO4SURF*ASURF_PO4    ! the input file provides the fluxes. units: mmol/s/m^2
                           CASE DEFAULT
                              surflux = -(PO4SURF-n1p(1))*vrelax ! units: mmol/s/m^2
                      END SELECT 

!
                      botflux = jbotN1p(1)/SEC_PER_DAY   
!
!                -----NITRATE-----
!
                 case (ppN3n)
!
                      SELECT CASE (NUTSBC_MODE)
                           CASE (1)
                              surflux = -NO3SURF*ASURF_NO3    ! the input file provides the fluxes. units: mmol/s/m^2
                           CASE DEFAULT
                              surflux = -(NO3SURF-n3n(1))*vrelax ! units: mmol/s/m^2
                      END SELECT 

                      botflux = jbotN3n(1)/SEC_PER_DAY
!
!                -----AMMONIUM-----
!
                 case (ppN4n)
!
                      SELECT CASE (NUTSBC_MODE)
                           CASE (1)
                              surflux = -NH4SURF*ASURF_NH4    ! the input file provides the fluxes. units: mmol/s/m^2
                           CASE DEFAULT
                              surflux = -(NH4SURF-n4n(1))*vrelax ! units: mmol/s/m^2
                      END SELECT 
!
                      botflux = jbotN4n(1)/SEC_PER_DAY
!
!                -----SILICATE-----
!
                 case (ppN5s)
!
                      SELECT CASE (NUTSBC_MODE)
                           CASE (1)
                              surflux = -SIO4SURF*ASURF_SIO4    ! the input file provides the fluxes. units: mmol/s/m^2
                           CASE DEFAULT
                              surflux = -(SIO4SURF-n5s(1))*vrelax ! units: mmol/s/m^2
                      END SELECT 
!
!               ***************************************************
!               ***************************************************
!               **                                               **
!               ** The simple benthic return module used in this **
!               ** implemantation has no silicate regenerated in **
!               ** the sediments                                 **
!               **                                               **
!               ***************************************************
!               ***************************************************
!
                    botflux = jbotN5s(1)/SEC_PER_DAY
!
                case default
!
!               ***************************************************
!               ***************************************************
!               **                                               **
!               ** The surface and bottom fluxes for the         **
!               ** BFM state variables not considered above      **
!               ** is nil.                                       **
!               **                                               **
!               ***************************************************
!               ***************************************************
!
                      surflux = ZERO
!
                      botflux = ZERO
!
          end select
!
!         *********************************************
!         *********************************************
!         **                                         **
!         ** Computing sedimentation fluxes          **
!         **                                         **
!         *********************************************
!         *********************************************
!
          select case ( m )
!
!                -----PARTICULATE ORGANIC MATTER-----
!
                 case (ppR6c:ppR6s)
!
                       do k = 1 , KB - 1
!
                          sink(k) = -sediR6(k)/SEC_PER_DAY
!
                       end do
!
!                      ***********************************************
!                      ***********************************************
!                      **                                           **
!                      ** Set the  Particulate organic matter       **
!                      ** sinking velocity (burial velocity) at the **
!                      ** water-sediment interface.                 **
!                      **                                           **
!                      ***********************************************
!                      ***********************************************
!
                       sink(kb)=-p_burvel_R6/SEC_PER_DAY
                       if(ABS(sink(kb-1)).lt.ABS(sink(kb))) &
                       sink(kb) = sink(kb-1)
!
                  end select
!
!         -----PHYTOPLANKTON-----
!
          if (m.GE.ppP1c .AND. m.LE.ppP4l) then
!
                     select case (m)
!
                               case (ppP1c:ppP1s)
!
                                     n = iiP1
!
                               case (ppP2c:ppP2l)
!
                                     n = iiP2
!
                               case (ppP3c:ppP3l)
!
                                     n = iiP3
!
                               case (ppP4c:ppP4l)
!
                                     n = iiP4
!
                     end select
!
                         do k = 1 , KB - 1
!
                            sink(k) = -sediPPY(k,n)/SEC_PER_DAY
!
                         enddo
!
!
!                      ***********************************************
!                      ***********************************************
!                      **                                           **
!                      ** Set the  Phytoplankton sinking velocity   **
!                      ** (burial velocity) at the water-sediment   **
!                      ** interface.                                **
!                      **                                           **
!                      ***********************************************
!                      ***********************************************
!
                       sink(kb)=-p_burvel_PI/SEC_PER_DAY
                       if(ABS(sink(kb-1)).lt.ABS(sink(kb))) &
                       sink(kb) = sink(kb-1)
!
        end if

!
!                 -----SINKING: UPSTREAM VERTICAL ADVECTION-----
!
!                 *******************************************************
!                 *******************************************************
!                 **                                                   **
!                 ** The subroutine output the vertical advection      **
!                 ** (sinking) flux stored into argument ffbio.        **
!                 ** It is called for all the BFM pelagic state        **
!                 ** variables, but for sink=0.0 it outputs ffbio=0.0  **
!                 **                                                   **
!                 *******************************************************
!                 *******************************************************
!

                  W = sink * SEDI_FACTOR(MONTH_OF_SIMULATION)
                  IF (USE_W_PROFILE) THEN
                     W = sink + W_PROFILE
                  END IF
                  IF (AssignPelBenFluxesInBFMFlag) THEN
                     ! the fluxes at the bottom are already accounted for in D3SOURCE and D2SOURCE_BEN
                     botflux = 0
                     W(KB:) = 0
                  END IF
                  CALL VERT_ADV_UPWIND_SCHEME(FBIO, W, DZ, KB, H, ADVTND) ! computing the advection tendency on the central time step

!
!                 ----- GETTING THE LATERAL RELAXATION COEFF.
!
                  SELECT CASE (m)
                     CASE (ppN1p)
                        LAMBDA = L_PO4(MONTH_OF_SIMULATION)
                     CASE (ppN3n)
                        LAMBDA = L_NO3(MONTH_OF_SIMULATION)
                     CASE (ppN5s)
                        LAMBDA = L_SIO4(MONTH_OF_SIMULATION)
                     CASE (ppO2o)
                        LAMBDA = L_O2(MONTH_OF_SIMULATION)
                     CASE DEFAULT
                        LAMBDA = L_X(MONTH_OF_SIMULATION)
                  END SELECT
                  LFLUX = LAMBDA*fbio ! estimating the lateral flux. Do it at the central time step!
                  do K=1,KB-1
!            
                   ! Integration in time: leap-frog of advection, lateral flux and reactive term
                    ffbio(k)=fbbio(k)+DTI2*(advtnd(k)+D3SOURCE(k,m)+LFLUX(K))
!            
                  end do

!
!     *******************************************************************
!     *******************************************************************
!     **                                                               **
!     ** COMPUTE VERTICAL DIFFUSION AND TERMINATE COMPUTATION WITH     **
!     ** IMPLICIT TIME INTEGRATION                                     **
!     **                                                               **
!     *******************************************************************
!     *******************************************************************
!
      CALL PROF_TRACERS(ffbio,surflux,botflux,ZERO,ZERO,NBCBFM,DTI2,NTP,UMOLBFM)
!
!     -----CLIPPING......IF NEEDED-----
!
      do k=1, KB-1
!
            ffbio(k)=max(p_small,ffbio(k))
      end do
!
!     ---MIX SOLUTIONS AND RESTORE TIME SEQUENCE-----
!
     do n = 1, KB-1
!
      D3STATEB(n,m)=fbio(n)+0.5_RLEN*smoth*(ffbio(n)+ fbbio(n)-2.0_RLEN*fbio(n))
      D3STATE(n,m)=ffbio(n)
!
     enddo
!

   enddo ! loop over NO_D3_BOX_STATES
!
end subroutine vert_integration

