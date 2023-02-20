!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PelagicBenthicCoupling
!
! DESCRIPTION
!   This routine solve the exchanges of OMT and Inorganic Nutrients 
!   between Pelagic and Benthic Systems, dealing with the following : 
!   - bottom deposition of phytoplankton (PPY) and detritus (OMT)
!   - uptake from pelagic states due to filtering organisms (and OMT release)
!   - buffering of nutrients in sediment under very-low OMT fluxes
!   - solve fluxes of OMT and Nutrients between Pelagic and Benthic Systems 
!
!  NOTES :
!  1) If Bethic system has active filter feeders group at this stage
!     bottom fluxes of R6 already contains OMT fluxes due to either 
!     excretion and/or pseudo faeces production.
!  2) Burial velocity is defined in ModuleSettling.F90, which controls 
!     the inflow rate of OMT and phytoplankton from the pelagic to benthic
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! INCLUDE
#include "DEBUG.h"
#include "INCLUDE.h"
!
! INTERFACE
   subroutine PelagicBenthicCoupling
!
! USES
   use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
   use mem
#else
   use mem, ONLY: D3STATE, PELBOTTOM, PhytoPlankton, NO_BOXES_XY, Depth,  &
           PhytoPlankton, ppPhytoPlankton, iiPhytoPlankton,               &
           MicroZooPlankton, ppMicroZooPlankton, iiMicroZooPlankton,      &
           R1c, R2c, R3c, R6c, R1n, R6n, R1p, R6p, R6s,                   &
           ppR1c, ppR6c, ppR1n, ppR6n, ppR1p, ppR6p, ppR6s, ppR2c, ppR3c, &
           ppO2o, ppN1p, ppN3n, ppN4n, ppN5s, ppN6r,                      &
           jbotR6c, jbotR6n, jbotR6p, jbotR6s, jbotR1c, jbotR1n, jbotR1p, &
           jbotO2o, jbotN1p, jbotN3n, jbotN4n, jbotN5s, jbotN6r, jbotR3c, &
           sediPPY, sediR2, sediR3, sediR6,                               &
           qpcPPY, qncPPY, qscPPY, qlcPPY, qncMIZ, qpcMIZ, qccPPY,        &
           iiC, iiN, iiP, iiL, iiS, iiLastElement, iiBen, iiPel,          &
           flux, flux_vector
#ifdef INCLUDE_PELFE
   use mem, ONLY: iiF, R6f, ppR6f, jbotR6f, ppR1f, jbotR1f, qfcPPY
#endif
#ifdef INCLUDE_PELCO2
   use mem, ONLY: sediO5, O5c, jbotO5c, ppO5c
#if defined INCLUDE_BENCO2 || ( ! defined BENTHIC_BIO && ! defined BENTHIC_FULL )
  use mem, ONLY: ppO3h, ppO3c, jbotO3c, jbotO3h
#endif
#endif
  use mem, ONLY: Q6c, Q6n, Q6p, Q6s, Q1c, Q1n, Q1p
  use mem, ONLY: ppQ6c, ppQ6n, ppQ6p, ppQ6s, ppQ1c, ppQ1n, ppQ1p, &
                 ppQ16c, ppQ16n, ppQ16p, ppQ16s
#endif
   use mem_PelSinkSet
   use constants, ONLY: MW_C
   use mem_Param, ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p , p_small,   &
                 CalcBenthicFlag, AssignPelBenFluxesInBFMFlag
   use mem_PelBac, ONLY: p_suhR1, p_sulR1, p_suR2, p_suR6,        &
                         p_qncPBA, p_qpcPBA
#if defined BENTHIC_BIO || defined BENTHIC_FULL
   use mem, ONLY: jPIY3c, jZIY3c, ZI_Fc, jRIY3c, jRIY3n, jRIY3p, jRIY3s,  &
                  D1m, D6m, D7m, D8m, D9m, ppD6m, ppD7m, ppD8m, ppD9m
#endif
 use api_bfm, ONLY: BOTindices
 use time,    ONLY: GetDelta

   IMPLICIT NONE
 
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Local Vectors used to change the scope of group vectors
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton, lcl_MicroZooPlankton
   real(RLEN), dimension(:), allocatable :: jbotR2R6, jbotR2R1
   real(RLEN), dimension(:,:), allocatable :: jbotR6PPY, jbotR1PPY
   real(RLEN)            :: newDm(NO_BOXES_XY), burfrac(NO_BOXES_XY)
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Local Variables
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer  :: i, j, kbot, Box 
   real(RLEN) :: sedi, c, p, s, uptake, ruQc, ruQn, ruQp, ruQs, ruQf, ruQl, Delta
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   allocate(jbotR1PPY(NO_BOXES_XY,iiLastElement), jbotR2R1(NO_BOXES_XY),   &
            jbotR6PPY(NO_BOXES_XY,iiLastElement), jbotR2R6(NO_BOXES_XY) )
   !
   jbotR2R1  = ZERO ; jbotR2R6  = ZERO
   jbotR1PPY = ZERO ; jbotR6PPY = ZERO
   burfrac = ZERO
   !
   ! At this point the Pelagic bottom arrays already contain :
   ! - Inorganic nutrients fluxes from benthic return/remin/nutrients schemes
   ! - Form Filter-Feeders (for BENTHIC_BIO or BENTHIC_FULL)
   !   excretion and/or pseudofaeces are contained in jbotR6x, and
   !   inorganic nutrients form fixed stoichiometry adjustment in jbotN4n,jbotN1p
   !
   ! loop over the number of bottom boxes
   BOXES_XY_LOOP : do Box = 1,NO_BOXES_XY
      kbot = BOTindices(Box)

#if defined BENTHIC_BIO || defined BENTHIC_FULL
   ! 
   ! ------------------------------------------------------------------
   ! GRAZING on pelagic states due to benthic organisms (Filter-Feeders)
   ! ------------------------------------------------------------------
   !
      if ( CalcBenthicFlag ) then   
         !
         ! Phytoplankton
         do i = 1 , iiPhytoPlankton
            lcl_PhytoPlankton => PhytoPlankton(i,iiC)
            ! Check amount of carbon, this may generate a negligible,small mass loss
            if ( lcl_PhytoPlankton(kbot) > p_small) then
               uptake = jPIY3c(Box,i)
               j = ppPhytoPlankton(i,iiC)
               PELBOTTOM(Box,j) = - uptake
               j = ppPhytoPlankton(i,iiN) 
               PELBOTTOM(Box,j) = - uptake*qncPPY(kbot,i)
               j = ppPhytoPlankton(i,iiP)
               PELBOTTOM(Box,j) = - uptake*qpcPPY(kbot,i)
               j = ppPhytoPlankton(i,iiL)
               PELBOTTOM(Box,j) = - uptake*qlcPPY(kbot,i)
               j = ppPhytoPlankton(i,iiS)
               if ( j > 0) &
                  PELBOTTOM(Box,j) = - uptake*qscPPY(kbot,i)
#ifdef INCLUDE_PELFE
               j = ppPhytoPlankton(i,iiF)
               if ( j > 0) & 
                  PELBOTTOM(Box,j) = - uptake*qfcPPY(kbot,i)
#endif
            end if
         end do
         !
         ! MicroZooplankton only
         do i = 1 , iiMicroZooPlankton
            lcl_MicroZooPlankton => MicroZooPlankton(i,iiC) ; c = lcl_MicroZooPlankton(kbot)
            uptake = jZIY3c(Box) * c / ZI_Fc(Box)
            ! Check amount of carbon, this may generate a negligible,small mass loss
            if ( c < p_small) uptake = ZERO
            j = ppMicroZooPlankton(i,iiC)
            PELBOTTOM(Box,j) = - uptake
            j = ppMicroZooPlankton(i,iiN)
            if ( j > 0 ) &
              PELBOTTOM(Box,j) = -uptake* qncMIZ(kbot,i)
            j = ppMicroZooPlankton(i,iiP)
            if ( j > 0 ) &
              PELBOTTOM(Box,j) = -uptake* qpcMIZ(kbot,i)
         end do
         !
         ! Detritus (net flux= uptake - excretion of food : flux may be negative!)
         call flux(kbot, iiPel, ppR6c, ppR6c, -(jRIY3c(Box)/Depth(kbot)) )
         call flux(kbot, iiPel, ppR6n, ppR6n, -(jRIY3n(Box)/Depth(kbot)) )
         call flux(kbot, iiPel, ppR6p, ppR6p, -(jRIY3p(Box)/Depth(kbot)) )
         call flux(kbot, iiPel, ppR6s, ppR6s, -(jRIY3s(Box)/Depth(kbot)) )
         !
      endif
#endif
   !
   ! ------------------------------------------------------------------
   ! PHYSICAL SETTLING for PPY and OMT 
   ! ------------------------------------------------------------------
   !
      ! Phytoplankton
      do i = 1 , ( iiPhytoPlankton)
         sedi  =   sediPPY(kbot,i)
         if ( sedi> ZERO ) then
            lcl_PhytoPlankton => PhytoPlankton(i,iiC)
            ! carbon
            j=ppPhytoPlankton(i,iiC)
            ruQc = sedi* lcl_PhytoPlankton(kbot)
            PELBOTTOM(Box,j)   = PELBOTTOM(Box,j) - ruQc
            jbotR1PPY(Box,iiC) = jbotR1PPY(Box,iiC) - ruQc * p_pe_R1c
            jbotR6PPY(Box,iiC) = jbotR6PPY(Box,iiC) - ruQc * (ONE - p_pe_R1c)
            ! nitrogen
            j=ppPhytoPlankton(i,iiN)
            ruQn = ruQc * qncPPY(kbot,i)
            PELBOTTOM(Box,j)   = PELBOTTOM(Box,j) - ruQn
            jbotR1PPY(Box,iiN) = jbotR1PPY(Box,iiN) - ruQn * p_pe_R1n
            jbotR6PPY(Box,iiN) = jbotR6PPY(Box,iiN) - ruQn * (ONE - p_pe_R1n)
            ! phosphorus
            j=ppPhytoPlankton(i,iiP)
            ruQp = ruQc * qpcPPY(kbot,i)
            PELBOTTOM(Box,j)   = PELBOTTOM(Box,j) - ruQp
            jbotR1PPY(Box,iiP) = jbotR1PPY(Box,iiP) - ruQp * p_pe_R1p
            jbotR6PPY(Box,iiP) = jbotR6PPY(Box,iiP) - ruQp * (ONE - p_pe_R1p)
            ! chlorophyll (stored but not used)
            j=ppPhytoPlankton(i,iiL)
            ruQl = ruQc * qlcPPY(kbot,i)
            j=ppPhytoPlankton(i,iiL)
            PELBOTTOM(Box,j)=  PELBOTTOM(Box,j) - ruQl
            ! silica
            j=ppPhytoPlankton(i,iiS)
            if ( j> 0) then
               ruQs = ruQc * qscPPY(kbot,i)
               PELBOTTOM(Box,j)   = PELBOTTOM(Box,j) - ruQs
               jbotR6PPY(Box,iiS) = jbotR6PPY(Box,iiS) - ruQs
            end if
#ifdef INCLUDE_PELFE
            j=ppPhytoPlankton(i,iiF)
            if ( j> 0) then
               ruQf = ruQc * qfcPPY(kbot,i)
               PELBOTTOM(Box,j)   = PELBOTTOM(Box,j) - ruQf
               jbotR6PPY(Box,iiF) = jbotR6PPY(Box,iiF) - ruQf ! split also to R1 ?
            end if
#endif
         end if
      enddo

      ! Particulate OM (R6)
      ruQc = sediR6(kbot) * R6c(kbot) 
      ruQn = sediR6(kbot) * R6n(kbot)
      ruQp = sediR6(kbot) * R6p(kbot)
      ruQs = sediR6(kbot) * R6s(kbot)
      jbotR6c(Box) = jbotR6c(Box) - ruQc
      jbotR6n(Box) = jbotR6n(Box) - ruQn
      jbotR6p(Box) = jbotR6p(Box) - ruQp
      jbotR6s(Box) = jbotR6s(Box) - ruQs
#ifdef INCLUDE_PELFE
      ruQf = sediR6(kbot) * R6f(kbot)
      jbotR6f(Box) = jbotR6f(Box) - ruQf
#endif
      !
      ! Intermediate degradable OM (R2) fluxes to R6 and R1 
      if ( sediR2(kbot) > ZERO) then
         ruQc = sediR2(kbot) * R2c(kbot)
         ! Actual degradability of LOC ( dependent of quotum NC,PC) and 
         ! compare with bacterial preferencial quota (first bacteria group)
         p = min(ONE, R1n(kbot)/(R1c(kbot)* p_qncPBA(1)), &
                                  R1p(kbot)/(R1c(kbot)*p_qpcPBA(1)))
         ! Calculate actual degradability of R1
         s= p_suhR1(1)*p-p_sulR1(1)* (ONE-p)
         ! Calculate distribution factor for R2 between R1 and R6
         p=(p_suR2(1)-p_suR6(1))/(s-p_suR6(1))
         jbotR2R1(Box) = - ruQc *    p   
         jbotR2R6(Box) = - ruQc * (ONE-p)
      end if

      ! Refractory OM (R3) fluxes to R6/R16
      if ( sediR3(kbot) > ZERO) then
         ruQc = sediR3(kbot) * R3c(kbot)
         jbotR3c(Box) = jbotR3c(Box) - ruQc
      end if 

      !
#ifdef INCLUDE_PELCO2
      ! burial rate of PIC (calcite)
      ruQc = sediO5(kbot) * O5c(kbot)
      jbotO5c(Box) = jbotO5c(Box) - ruQc
#endif
   !
   enddo BOXES_XY_LOOP
   !
   !
   if ( .NOT. AssignPelBenFluxesInBFMFlag) return
   !
   ! ------------------------------------------------------------------
   ! FLUXES for Pelagic-Benthic Coupling
   ! ------------------------------------------------------------------
   !
      !
      ! Pelagic Organic Matter
   OMT_XY_LOOP : do Box = 1,NO_BOXES_XY
      kbot = BOTindices(Box)
      call flux(kbot, iiPel, ppR6c, ppR6c, (jbotR6c(Box)/Depth(kbot)) )
      call flux(kbot, iiPel, ppR6n, ppR6n, (jbotR6n(Box)/Depth(kbot)) )
      call flux(kbot, iiPel, ppR6p, ppR6p, (jbotR6p(Box)/Depth(kbot)) )
      call flux(kbot, iiPel, ppR6s, ppR6s, (jbotR6s(Box)/Depth(kbot)) )
      call flux(kbot, iiPel, ppR1c, ppR1c, (jbotR1c(Box)/Depth(kbot)) )
      call flux(kbot, iiPel, ppR1n, ppR1n, (jbotR1n(Box)/Depth(kbot)) )
      call flux(kbot, iiPel, ppR1p, ppR1p, (jbotR1p(Box)/Depth(kbot)) )
#ifdef INCLUDE_PELFE
      call flux(kbot, iiPel, ppR1f, ppR1f, (jbotR1f(Box)/Depth(kbot)) )
      call flux(kbot, iiPel, ppR6f, ppR6f, (jbotR6f(Box)/Depth(kbot)) )
#endif
      c = jbotR2R6(Box) + jbotR2R1(Box)
      call flux(kbot, iiPel, ppR2c, ppR2c, ( c /Depth(kbot)) )
      call flux(kbot, iiPel, ppR3c, ppR3c, ( jbotR3c(Box)/Depth(kbot)) )
      !
      ! Detritus (R6,R1,R2,PPY) fluxes to sediment are directed via R6/R1
      ! R6 -> Q6 
      jbotR6c(Box) = jbotR6c(Box) + jbotR6PPY(Box,iiC) + jbotR2R6(Box)
      jbotR6n(Box) = jbotR6n(Box) + jbotR6PPY(Box,iiN)
      jbotR6p(Box) = jbotR6p(Box) + jbotR6PPY(Box,iiP)
      jbotR6s(Box) = jbotR6s(Box) + jbotR6PPY(Box,iiS)
#ifdef INCLUDE_PELFE
      jbotR6f(Box) = jbotR6f(Box) + jbotR6PPY(Box,iiF)
#endif
      ! R1 -> Q1
      jbotR1c(Box) = jbotR1c(Box) + jbotR1PPY(Box,iiC) + jbotR2R1(Box)
      jbotR1n(Box) = jbotR1n(Box) + jbotR1PPY(Box,iiN)
      jbotR1p(Box) = jbotR1p(Box) + jbotR1PPY(Box,iiP)
#ifdef INCLUDE_PELFE
      jbotR1f(Box) = jbotR1f(Box) + jbotR1PPY(Box,iiF)
#endif
      !
   enddo OMT_XY_LOOP

#if defined BENTHIC_BIO || defined BENTHIC_FULL
   ! ------------------------------------------------------------------
   ! Control effective flux of jbotR6x (bio+physical) to benthic system
   ! and compensate with nutrients uptake to avoid Q6 pools become zero 
   ! ------------------------------------------------------------------
      call ControlBennutBuffersDynamics
#endif
      !
   FLUXES_XY_LOOP : do Box = 1,NO_BOXES_XY
      kbot = BOTindices(Box)
      !
      ! Inorganic Nutrients
      call flux(kbot, iiPel, ppO2o, ppO2o, jbotO2o(Box)/Depth(kbot) )
      call flux(kbot, iiPel, ppN1p, ppN1p, jbotN1p(Box)/Depth(kbot) )
      call flux(kbot, iiPel, ppN3n, ppN3n, jbotN3n(Box)/Depth(kbot) )
      call flux(kbot, iiPel, ppN4n, ppN4n, jbotN4n(Box)/Depth(kbot) )
      call flux(kbot, iiPel, ppN5s, ppN5s, jbotN5s(Box)/Depth(kbot) )
      call flux(kbot, iiPel, ppN6r, ppN6r, jbotN6r(Box)/Depth(kbot) )
#if defined INCLUDE_PELCO2
      call flux(kbot, iiPel, ppO5c, ppO5c, jbotO5c(Box)/Depth(kbot) )
#if defined INCLUDE_BENCO2 || ( ! defined BENTHIC_BIO && ! defined BENTHIC_FULL )
      call flux(kbot, iiPel, ppO3c, ppO3c, jbotO3c(Box)/Depth(kbot) )
      call flux(kbot, iiPel, ppO3h, ppO3h, jbotO3h(Box)/Depth(kbot) )
#endif
#endif
      ! PhytoPlankton
      do i = 1 , ( iiPhytoPlankton)
         j = ppPhytoPlankton(i,iiC)
         call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
         j = ppPhytoPlankton(i,iiN)
         call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
         j = ppPhytoPlankton(i,iiP)
         call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
         j = ppPhytoPlankton(i,iiL)
         call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
         j = ppPhytoPlankton(i,iiS)
         if ( j > 0 ) &
           call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
#ifdef INCLUDE_PELFE
         j = ppPhytoPlankton(i,iiF)
         if ( j > 0 ) &
           call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
#endif
      enddo
      !
      ! MicroZooplankton
      do i = 1 , ( iiMicroZooPlankton)
         j = ppMicroZooPlankton(i,iiC)
         call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
         j = ppMicroZooPlankton(i,iiN)
         if ( j > 0) &
           call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
         j = ppMicroZooPlankton(i,iiP)
         if ( j > 0) &
           call flux(kbot, iiPel, j, j, PELBOTTOM(Box,j)/Depth(kbot) )
      enddo
      !
   enddo FLUXES_XY_LOOP
   !
   ! ------------------------------------------------------------------
   ! FLUXES for Benthic-Pelagic Coupling
   ! ------------------------------------------------------------------
   !
   if ( CalcBenthicFlag ) then

      if ( R6DeepBurial ) then
         !
         ! Buried fraction of POC into deep sediment (Dunne et al,2007-GBC)
         burfrac = 0.013_RLEN + (0.53_RLEN * abs(p_small + jbotR6c(:)/MW_C)**2) /  &
                        ( 7.0_RLEN + abs(p_small + jbotR6c(:)/MW_C))**2
         !
         ! Burial od Pelagic Organic Matter
         call flux_vector( iiBen, ppQ16c,ppQ16c, -jbotR6c(:) * burfrac)
         call flux_vector( iiBen, ppQ16n,ppQ16n, -jbotR6n(:) * burfrac)
         call flux_vector( iiBen, ppQ16p,ppQ16p, -jbotR6p(:) * burfrac)
         call flux_vector( iiBen, ppQ16s,ppQ16s, -jbotR6s(:) * burfrac)
         call flux_vector( iiBen, ppQ16c,ppQ16c, -jbotO5c(:) )
         call flux_vector( iiBen, ppQ16c,ppQ16c, -jbotR3c(:) )
      else
         ! Here Divert O5c flux to surface Benthic OM for mass conservation 
         ! TL: it has to be included in benthic CSYS 
         call flux_vector( iiBen, ppQ6c,ppQ6c, -jbotO5c(:) )
         call flux_vector( iiBen, ppQ6c,ppQ6c, -jbotR3c(:) )

      endif

      !
      ! Sedimentation of Pelagic Organic Matter (Q6f still missing !!!)

      call flux_vector( iiBen, ppQ6c,ppQ6c, -jbotR6c(:) * (ONE - burfrac) )
      call flux_vector( iiBen, ppQ6n,ppQ6n, -jbotR6n(:) * (ONE - burfrac) )
      call flux_vector( iiBen, ppQ6p,ppQ6p, -jbotR6p(:) * (ONE - burfrac) )
      call flux_vector( iiBen, ppQ6s,ppQ6s, -jbotR6s(:) * (ONE - burfrac) )
      call flux_vector( iiBen, ppQ1c,ppQ1c, -jbotR1c(:) )
      call flux_vector( iiBen, ppQ1n,ppQ1n, -jbotR1n(:) )
      call flux_vector( iiBen, ppQ1p,ppQ1p, -jbotR1p(:) )

#if defined BENTHIC_BIO || defined BENTHIC_FULL
      !
      ! Changes in depth distribution of states due to OMT sedimentation

      Delta=GetDelta( )
      call RecalcPenetrationDepth( D1m(:), D6m(:), &
           -jbotR6c(:)*Delta, Q6c(:),newDm(:) )
      call flux_vector(iiBen, ppD6m,ppD6m,(newDM(:)- D6m(:))/Delta)
      call RecalcPenetrationDepth( D1m(:), D7m(:), &
           -jbotR6n(:)*Delta, Q6n(:),newDm(:) )
      call flux_vector(iiBen, ppD7m,ppD7m,(newDM(:)- D7m(:))/Delta)
      call RecalcPenetrationDepth( D1m(:), D8m(:), &
           -jbotR6p(:)*Delta, Q6p(:),newDm(:) )
      call flux_vector(iiBen, ppD8m,ppD8m,(newDM(:)- D8m(:))/Delta)
      call RecalcPenetrationDepth(D1m(:), D9m(:), &
           -jbotR6s(:)*Delta, Q6s(:),newDm(:) )
      call flux_vector(iiBen, ppD9m,ppD9m,(newDM(:)- D9m(:))/Delta)
#endif

   endif 
   !
   deallocate( jbotR1PPY, jbotR2R1, jbotR6PPY, jbotR2R6 )
   ! 
 
   return

end subroutine PelagicBenthicCoupling

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
