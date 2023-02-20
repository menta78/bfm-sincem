!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! ROUTINE: PelBacDynamics4
!
! DESCRIPTION
!   This process describes the dynamics of bacterioplankton used in CMCC-ESM2
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
  subroutine PelBacDynamics4(bac)
!
! USES
  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY:  R6c, R6n, R6p, R1c, R1n, R1p, R2c, O2o, N6r, BGE, &
    N4n, N1p, N3n, R3c, iiR1, iiR6, D3STATE
  use mem, ONLY: iiPelBacteria, ppPelBacteria, iiC, iiN, iiP, ppR6c, ppR6n, ppR6p, &
    ppR1c, ppR1n, ppR1p, ppR2c, ppO2o, ppN6r, ppN4n, ppN1p, ppN3n, ppR3c, ETW, EPR, &
    qncPBA, qpcPBA, qpcOMT, qncOMT, NO_BOXES, iiBen, iiPel, flux_vector, quota_flux
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c
#endif
#endif
  use constants,  ONLY: MW_C, ONE_PER_DAY
  use mem_Param,  ONLY: p_pe_R1c, p_small
  use mem_PelBac
  use mem_globalfun,   ONLY: eTq, MM_power, insw, MM, fixratio
  use bfm_error_msg, ONLY: bfm_error

  IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: bac

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer       :: i
  integer       :: ppbacc, ppbacn, ppbacp
  integer, save :: first =0
  integer       :: AllocStatus
  integer,dimension(NO_BOXES)  :: limit
  real(RLEN)    :: p_dnc = 0.86_RLEN  ! N/C ratio used in denitrification (Paulmier et al., 2009)
  real(RLEN),allocatable,save,dimension(:) :: et,eO2,rd,rrc,rrn,rrt,  &
                                          ruR1c,ruR1n,ruR1p,ruR6c,ruR6p,ruR6n, &
                                          ruR2c, reR2c, ruR3c, reR3c, &
                                          ren,rep,rut,rugch,rum,run,rug,rdru,R3ex,R1rut, R6ch, &
                                          iN1p,iNIn,iN, eN1p,eN4n, bacc, &
                                          tfluxC, tfluxN, tfluxP, pe_N4n, pe_N1p, pe_R6c,deplim
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     ALLOCATE ( et(NO_BOXES), eO2(NO_BOXES), rd(NO_BOXES),            &
        &       rrc(NO_BOXES), rrn(NO_BOXES), rrt(NO_BOXES),          &
        &       ruR1c(NO_BOXES), ruR1n(NO_BOXES), ruR1p(NO_BOXES),    &
        &       ruR6c(NO_BOXES), ruR6p(NO_BOXES), ruR6n(NO_BOXES),    &
        &       ruR2c(NO_BOXES), reR2c(NO_BOXES), rdru(NO_BOXES),     &
        &       ruR3c(NO_BOXES), reR3c(NO_BOXES), rugch(NO_BOXES),    &
        &       ren(NO_BOXES), rep(NO_BOXES), rug(NO_BOXES),          &
        &       rut(NO_BOXES), rum(NO_BOXES), run(NO_BOXES),          &
        &       iN1p(NO_BOXES), iNIn(NO_BOXES), iN(NO_BOXES),         &
        &       eN1p(NO_BOXES), eN4n(NO_BOXES), bacc(NO_BOXES),       &
        &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES), &
        &       pe_N4n(NO_BOXES), pe_N1p(NO_BOXES), pe_R6c(NO_BOXES), &
        &       R3ex(NO_BOXES), R1rut(NO_BOXES), R6ch(NO_BOXES),   &
        &       deplim(NO_BOXES), STAT = AllocStatus )
     IF( AllocStatus /= 0 ) call bfm_error('PelBacDynamics','Error allocating arrays')
     first=1
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Copy state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppbacc = ppPelBacteria(bac,iiC)
  ppbacn = ppPelBacteria(bac,iiN)
  ppbacp = ppPelBacteria(bac,iiP)
  bacc = D3STATE(:,ppbacc)

  ! Quota collectors
  tfluxC = ZERO
  tfluxN = ZERO
  tfluxP = ZERO
  BGE  = ZERO
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et = eTq(ETW(:), p_q10(bac))

  ! Pressure limitation (Tamburini, 2009 work) use exp function
  deplim = ONE
  if (p_dep(bac)> ZERO ) &
      deplim = MIN( ONE, (EPR(:)/p_dep(bac))**p_dep_exp(bac))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Oxygen environment: bacteria are both aerobic and anaerobic
  ! To provide a faster switching between the two metabolic pathways the
  ! oxygen regulating factor eO2 is cubic (eq. 19 in Vichi et al., 2004)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2 = MM_power(max(p_small,O2o(:)),  p_chdo(bac),3)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality:
  !  1. first order mortality: p_sd
  !  2. density dependent mortality (e.g. virus infection): p_sd2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd  =  ( p_sd(bac)*et + p_sd2(bac)*bacc ) * bacc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! External nutrient limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eN4n = MM(N4n(:), p_chn(bac))
  eN1p = MM(N1p(:), p_chp(bac))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
  ! Internal Nutrient limitation (always one for fixed stoichiometry)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iNIn = min(ONE, max(ZERO, qncPBA(:,bac)/p_qncPBA(bac)))  !Nitrogen
  iN1p = min(ONE, max(ZERO, qpcPBA(:,bac)/p_qpcPBA(bac)))  !Phosphorus
  iN   = min(iN1p, iNIn)

  ! DOC aggregation into POC
  rugch = ZERO
  rugch =  MIN(R1c(:), R6c(:) * p_sulR1(bac))
  call flux_vector(iiPel, ppR1c, ppR6c, rugch) 
  call flux_vector(iiPel, ppR1n, ppR6n, rugch*qncOMT(:,iiR1)) 
  call flux_vector(iiPel, ppR1p, ppR6p, rugch*qpcOMT(:,iiR1)) 

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rum = iN * p_sum(bac) * bacc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Total substrate from detritus pools
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  rut = p_small + R1c(:) + R6c(:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Gross carbon uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  R1rut = R1c(:)/rut
  ruR1c = rum * et * MM_power( R1c(:) , p_suhR1(bac), 3 ) * R1rut

  R6ch  = p_suR6(bac) *  bacc * (ONE - R1rut)
  ruR6c = rum * deplim * MM_power( R6c(:) , R6ch , 3)  * (ONE - R1rut)

  rug = ruR1c + ruR6c + p_small

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Carbon uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call quota_flux(iiPel, ppbacc, ppR1c, ppbacc, ruR1c, tfluxC)
  call quota_flux(iiPel, ppbacc, ppR6c, ppbacc, ruR6c, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR1n = qncOMT(:,iiR1)*ruR1c
  ruR6n = qncOMT(:,iiR6)*ruR6c
  call quota_flux(iiPel, ppbacn, ppR1n, ppbacn, ruR1n, tfluxN)
  call quota_flux(iiPel, ppbacn, ppR6n, ppbacn, ruR6n, tfluxN)

  ruR1p = qpcOMT(:,iiR1)*ruR1c
  ruR6p = qpcOMT(:,iiR6)*ruR6c
  call quota_flux(iiPel, ppbacp, ppR1p, ppbacp, ruR1p, tfluxP)
  call quota_flux(iiPel, ppbacp, ppR6p, ppbacp, ruR6p, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality redistributed on substrates proportionally to realized uptake
  ! In absence of substrate, mortality goes toward dissolved matter
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rdru = rd*(ONE-ruR6c/rug)
  call quota_flux(iiPel, ppbacc, ppbacc, ppR1c, rdru              , tfluxC)
  call quota_flux(iiPel, ppbacn, ppbacn, ppR1n, rdru*qncPBA(:,bac), tfluxN)
  call quota_flux(iiPel, ppbacp, ppbacp, ppR1p, rdru*qpcPBA(:,bac), tfluxP)

  rdru = rd*ruR6c/rug
  call quota_flux(iiPel, ppbacc, ppbacc, ppR6c, rdru              , tfluxC)
  call quota_flux(iiPel, ppbacn, ppbacn, ppR6n, rdru*qncPBA(:,bac), tfluxN)
  call quota_flux(iiPel, ppbacp, ppbacp, ppR6p, rdru*qpcPBA(:,bac), tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Aerobic and anaerobic respiration
  ! Anaerobic bacteria use NO3 instead of O2, with additional carbon cost.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity + Basal respiration
  rrc = (p_pu_ra(bac) * rug * p_pu_ea_R3(bac)) + p_srs(bac)* bacc* et
  ! Aerobic consumption of Oxygen
  call flux_vector( iiPel, ppO2o, ppO2o, -eO2*rrc/MW_C )
  ! Anaerobic consumption of Nitrate (denitrification-like)
  rrn = p_pu_ra_o(bac) * rug
  call flux_vector( iiPel, ppN3n, ppN4n, (ONE-eO2)* (rrc + rrn) / MW_C* p_dnc) 
  ! Total Respiration
  rrt = rrc + (ONE-eO2)*rrn
  call quota_flux( iiPel, ppbacc, ppbacc, ppO3c, rrt, tfluxC)

  ! I suppose is R3
  R3ex = p_pu_ra(bac) * rug * (ONE - p_pu_ea_R3(bac))
  call quota_flux( iiPel, ppbacc, ppbacc, ppR3c, R3ex, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Net Production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  run = rug - rrt

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Growth efficiency
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  WHERE (run.gt.p_small) BGE = run / (run+rrt)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes from bacteria (as in Vichi et al. 2007)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved Nitrogen dynamics
  ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
  ! This rate is assumed to occur with a timescale p_ruen=1 day
  ! and controlled with a Michaelis-Menten function
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ren = (qncPBA(:,bac) - p_qncPBA(bac))*bacc*p_ruen(bac)

  ! Fix-quota: actual quota comes from detritus uptake vs. realized growth
  if ( ppbacn == 0 ) & 
     ren = min( ZERO, (ruR1n+ruR6n) - p_qncPBA(bac)*run )

  call quota_flux(iiPel, ppbacn, ppbacn, ppN4n,       ren*insw( ren), tfluxN)
  call quota_flux(iiPel, ppbacn, ppN4n, ppbacn, -eN4n*ren*insw(-ren), tfluxN)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved Phosphorus dynamics
  ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
  ! This rate is assumed to occur with a timescale of p_ruep=1 day
  ! and controlled with a Michaelis-Menten function
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rep = (qpcPBA(:,bac) - p_qpcPBA(bac))*bacc*p_ruep(bac)

  ! Fix-quota: actual quota comes from detritus uptake vs. realized growth
  if ( ppbacp == 0 ) &
     rep = min(ZERO, (ruR1p+ruR6p) - p_qpcPBA(bac)*run )

  call quota_flux(iiPel, ppbacp, ppbacp, ppN1p,       rep*insw( rep), tfluxP)
  call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -eN1p*rep*insw(-rep), tfluxP)

  ! Fix-quota control
  if ( ppbacn == 0 .or. ppbacp == 0 ) then

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Eliminate the excess of the non-limiting constituent under fixed quota
     ! Determine release due to either C, P, or N limitation
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     call fixratio(tfluxC,tfluxN,tfluxP,qncPBA(:,bac),qpcPBA(:,bac),pe_R6c, pe_N4n, pe_N1p)

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Apply the correction terms depending on the limiting constituent
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     call flux_vector(iiPel, ppbacc, ppR6c, pe_R6c*(ONE-p_pe_R1c))
     call flux_vector(iiPel, ppbacc, ppR1c, pe_R6c*(p_pe_R1c)    )
     call flux_vector(iiPel, ppbacp, ppN1p, pe_N1p)
     call flux_vector(iiPel, ppbacn, ppN4n, pe_N4n)

  endif

  end subroutine PelBacDynamics4

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
