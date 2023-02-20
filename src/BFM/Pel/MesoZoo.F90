!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: MesoZoo
!
! DESCRIPTION
!   This submodel describes the carbon dynamics and associated
!   nutrient dynamics in mesozooplankton 
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
  subroutine MesoZooDynamics(zoo)
!
! USES
  use global_mem, ONLY:RLEN, ONE, ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D3STATE, O2o, N1p, N4n, R6c, R6p, R2c, &
    R6n, PhytoPlankton, MicroZooPlankton, MesoZooPlankton
  use mem, ONLY: Depth, ppO2o, ppN1p, ppN4n, ppR6c, ppR6n, ppR6p, ppR6s, &
    ppPhytoPlankton, ppMicroZooPlankton, ppMesoZooPlankton, ETW, &
    qncPPY, qpcPPY, qlcPPY, qscPPY, qncMIZ, qpcMIZ, qncMEZ, qpcMEZ, iiPhytoPlankton, &
    iiMicroZooPlankton, iiMesoZooPlankton, iiC, iiN, iiP, iiL, iiS, NO_BOXES, &
    iiBen, iiPel, flux_vector, quota_flux
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c, ppO5c, ppO3h, qccPPY
#endif
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF, qfcPPY, ppR6f
#endif
#endif
  use mem_Param,  ONLY: p_small
  use bfm_error_msg, ONLY: bfm_error
  use constants,ONLY: MIN_VAL_EXPFUN, MW_C, C2ALK
  use mem_MesoZoo
  use mem_globalfun,   ONLY: eTq, MM, MM_power, fixratio

  IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: zoo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer  :: ppzooc, ppzoon, ppzoop
  integer, save :: first =0
  integer,dimension(NO_BOXES)  :: limit
  real(RLEN),allocatable,save,dimension(:) :: sut,temp_p,temp_n,rumc,rugc,eo,  &
                                       et,rrs_c,rrs_n,rrs_p,rut_c, &
                                       rut_n,rut_p,rd_c,rd_n,rd_p,sdo,rdo_c,  &
                                       rdo_n,rdo_p,ret_c,ret_n,ret_p,ru_c, &
                                       ru_n,ru_p,pu_e_n,pu_e_p,prI,pe_R6c

  real(RLEN),allocatable,save,dimension(:) :: pe_N1p,pe_N4n,ruPPYc,ruMIZc,ruMEZc,rq6c, &
                                       rq6n,rq6p,rrc,ren,rep,tfluxC, tfluxN, tfluxP,     &
                                       zooc,zoop,zoon
  real(RLEN),allocatable,save,dimension(:,:) :: PPYc,MIZc,MEZc
  real(RLEN),allocatable,save,dimension(:) :: net,r
  integer :: AllocStatus
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     ALLOCATE ( PPYc(NO_BOXES,iiPhytoPlankton), MIZc(NO_BOXES,iiMicroZooPlankton), &
        &       MEZc(NO_BOXES,iiMesoZooPlankton),  &
        &       zooc(NO_BOXES), zoop(NO_BOXES), zoon(NO_BOXES),        &
        &       rep(NO_BOXES), ren(NO_BOXES), rrc(NO_BOXES),           &
        &       rq6p(NO_BOXES), rq6n(NO_BOXES), rq6c(NO_BOXES),        &
        &       ruPPYc(NO_BOXES), ruMIZc(NO_BOXES), ruMEZc(NO_BOXES),  &
        &       pe_N4n(NO_BOXES), pe_N1p(NO_BOXES), pe_R6c(NO_BOXES),  &
        &       prI(NO_BOXES), pu_e_p(NO_BOXES), pu_e_n(NO_BOXES),     &
        &       ru_p(NO_BOXES), ru_n(NO_BOXES), ru_c(NO_BOXES),        &
        &       ret_p(NO_BOXES), ret_n(NO_BOXES), ret_c(NO_BOXES),     &
        &       rdo_p(NO_BOXES), rdo_n(NO_BOXES), rdo_c(NO_BOXES),     &
        &       rd_p(NO_BOXES) , rd_n(NO_BOXES) , rd_c(NO_BOXES) ,     &
        &       rut_p(NO_BOXES), rut_n(NO_BOXES), rut_c(NO_BOXES),     &
        &       eo(NO_BOXES), sdo(NO_BOXES), et(NO_BOXES), sut(NO_BOXES),   &
        &       rumc(NO_BOXES), rugc(NO_BOXES), net(NO_BOXES), r(NO_BOXES), &
        &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES),  &
        &       temp_p(NO_BOXES), temp_n(NO_BOXES),                    &
        &      STAT = AllocStatus )

     IF( AllocStatus /= 0 ) call bfm_error('MesoZooDynamics','Error allocating arrays')
     first=1
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Copy state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppzooc = ppMesoZooPlankton(zoo,iiC)
  ppzoon = ppMesoZooPlankton(zoo,iiN)
  ppzoop = ppMesoZooPlankton(zoo,iiP)

  zooc = D3STATE(:,ppzooc)
  zoon = zooc * qncMEZ(:,zoo)
  zoop = zooc * qpcMEZ(:,zoo)

  ! Quota collectors
  tfluxC = ZERO
  tfluxN = ZERO
  tfluxP = ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature and oxygen response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eo = MM_power(max(p_small,O2o(:)), p_clO2o(zoo),3)
  et = eTq(ETW(:), p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total potential food given the non-dim prey availability
  ! with loops over all LFGs.
  ! There is no parameter for capture efficiency in mesozooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rumc = ZERO
  do i = 1, iiPhytoPlankton
    PPYc(:,i) = p_paPPY(zoo,i)*PhytoPlankton(i,iiC)
    rumc = rumc + PPYc(:,i)
  end do
  do i = 1, iiMicroZooPlankton
    MIZc(:,i) = p_paMIZ(zoo,i)*MicroZooPlankton(i,iiC)
    rumc = rumc + MIZc(:,i)
  end do
  do i = 1, iiMesoZooPlankton
    MEZc(:,i) = p_paMEZ(zoo,i)*MesoZooPlankton(i,iiC)
    rumc = rumc + MEZc(:,i)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake rate (eq 38 Vichi et al. 2007) and 
  ! specific uptake rate considering potentially available food (sut)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  = et*p_sum(zoo)*MM(p_vum(zoo)*rumc, p_sum(zoo))*zooc
  sut = rugc/(p_small + rumc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptakes from every LFG
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rut_c = ZERO
  rut_n = ZERO
  rut_p = ZERO

  ! Phytoplankton
  do i = 1, iiPhytoPlankton
    ruPPYc = sut*PPYc(:,i)
    call flux_vector(iiPel, ppPhytoPlankton(i,iiC), ppzooc, ruPPYc)
    call flux_vector(iiPel, ppPhytoPlankton(i,iiN), ppzoon, ruPPYc*qncPPY(:,i))
    call flux_vector(iiPel, ppPhytoPlankton(i,iiP), ppzoop, ruPPYc*qpcPPY(:,i))
    rut_c = rut_c + ruPPYc
    rut_n = rut_n + ruPPYc*qncPPY(:,i)
    rut_p = rut_p + ruPPYc*qpcPPY(:,i)
    ! Chl is transferred to the infinite sink
    call flux_vector(iiPel, ppPhytoPlankton(i,iiL), &
               ppPhytoPlankton(i,iiL), -ruPPYc*qlcPPY(:,i))
    ! silicon constituent is transferred to biogenic silicate
    if ( ppPhytoPlankton(i,iiS) .gt. 0 ) & 
       call flux_vector(iiPel, ppPhytoPlankton(i,iiS), ppR6s, ruPPYc*qscPPY(:,i))

#ifdef INCLUDE_PELFE
    ! Fe constituent is transferred to particulate iron
    if ( ppPhytoPlankton(i,iiF) .gt. 0 ) & 
       call flux_vector(iiPel, ppPhytoPlankton(i,iiF), ppR6f, ruPPYc*qfcPPY(:,i))
#endif

#if defined INCLUDE_PELCO2
    ! PIC (calcite/aragonite) production associated to the grazed biomass
    ! The idea in PISCES is that the calcite flux exists only when associated
    ! to a carbon release from phytoplankton (there is no calcite storage in phyto)
    ! Use the realized rain ratio for each phytoplankton species and assume
    ! that only a portion is egested
    ! Calcite production is parameterized as a flux between DIC and PIC
    ! that affects alkalinity
    call flux_vector( iiPel, ppO3c,ppO5c, p_pecaco3(zoo)*ruPPYc*qccPPY(:,i))
    call flux_vector( iiPel, ppO3h,ppO3h, -C2ALK*p_pecaco3(zoo)*ruPPYc*qccPPY(:,i))
#endif

  end do

  ! Microzooplankton
  do i = 1, iiMicroZooPlankton
    ruMIZc = sut*MIZc(:,i)
    call flux_vector(iiPel, ppMicroZooPlankton(i,iiC), ppzooc, ruMIZc)
    call flux_vector(iiPel, ppMicroZooPlankton(i,iiN), ppzoon, ruMIZc*qncMIZ(:,i))
    call flux_vector(iiPel, ppMicroZooPlankton(i,iiP), ppzoop, ruMIZc*qpcMIZ(:,i))
    rut_c = rut_c + ruMIZc
    rut_n = rut_n + ruMIZc*qncMIZ(:,i)
    rut_p = rut_p + ruMIZc*qpcMIZ(:,i)
  end do

  ! Mesozooplankton
  do i = 1, iiMesoZooPlankton
    ruMEZc = sut*MEZc(:, i)
    ! Note that intra-group predation (cannibalism) is not added as a flux
    if ( i/= zoo ) then
      call flux_vector(iiPel, ppMesoZooPlankton(i,iiC), ppzooc, ruMEZc)
      call flux_vector(iiPel, ppMesoZooPlankton(i,iiN), ppzoon, ruMEZc*qncMEZ(:,i))
      call flux_vector(iiPel, ppMesoZooPlankton(i,iiP), ppzoop, ruMEZc*qpcMEZ(:,i))
    end if
    rut_c = rut_c + ruMEZc
    rut_n = rut_n + ruMEZc*qncMEZ(:,i)
    rut_p = rut_p + ruMEZc*qpcMEZ(:,i)
  end do

  ! Note that tfluxC include also intra-group predation
  tfluxC = tfluxC + rut_c 
  tfluxN = tfluxN + rut_n
  tfluxP = tfluxP + rut_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity respiration and basal metabolism
  ! First compute the the energy cost of ingestion
  ! 1 - assimilation - egestion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  prI = ONE - p_puI(zoo) - p_peI(zoo)
  rrc = prI*rut_c + p_srs(zoo)*et*zooc
  call flux_vector(iiPel, ppO2o, ppO2o, -rrc/MW_C)
  call quota_flux(iiPel, ppzooc, ppzooc, ppO3c, rrc, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Specific rates of low oxygen mortality
  ! and Density dependent mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rdo_c = p_sdo(zoo)*(ONE-eo)*et*zooc
  rd_c  = p_sd(zoo)*zooc**p_sds(zoo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total egestion including pellet production 
  ! Eq. 40 and 44 Vichi et al. 2007
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rq6c = p_peI(zoo)*rut_c + rdo_c + rd_c
  rq6n = p_peI(zoo)*rut_n + qncMEZ(:,zoo)*(rdo_c + rd_c)
  rq6p = p_peI(zoo)*rut_p + qpcMEZ(:,zoo)*(rdo_c + rd_c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient remineralization 
  ! basal metabolism + excess of non-limiting nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ren = p_srs(zoo)*et*eo*zoon 
  rep = p_srs(zoo)*et*eo*zoop 
  call quota_flux(iiPel, ppzoop, ppzoop, ppN1p, rep, tfluxP)
  call quota_flux(iiPel, ppzoon, ppzoon, ppN4n, ren, tfluxN)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes to particulate organic matter
  ! Add the correction term for organic carbon release in case of
  ! nutrient limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call quota_flux(iiPel, ppzooc, ppzooc,ppR6c, rq6c, tfluxC)
  call quota_flux(iiPel, ppzoon, ppzoon,ppR6n, rq6n, tfluxN)
  call quota_flux(iiPel, ppzoop, ppzoop,ppR6p, rq6p, tfluxP)

  if ( ppzoon > 0 .and. ppzoop > 0 ) then
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Check the assimilation rate for Carbon, Nitrogen and Phosphorus
     ! Note that activity respiration does not involve nutrient utilization
     ! so more nutrients than carbon are taken up.
     ! Then compute P:C and N:C ratios in the assimilation rate
     ! Eq 41 in Vichi et al. 2007 (there is an error in the denominator,
     ! the \Iota_c should be \Iota_i, with i=n,p)
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ru_c = p_puI(zoo)*rut_c
     ru_n = (p_puI(zoo) + prI)* rut_n
     ru_p = (p_puI(zoo) + prI)* rut_p
     pu_e_n  =   ru_n/( p_small+ ru_c)
     pu_e_p  =   ru_p/( p_small+ ru_c)

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Eliminate the excess of the non-limiting constituent
     ! Determine whether C, P or N is the limiting element and assign the
     ! value to variable limit
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     limit = iiC
     temp_p  = pu_e_p/qpcMEZ(:,zoo)
     temp_n  = pu_e_n/qncMEZ(:,zoo)

     WHERE ( temp_p<temp_n .OR. abs(temp_p-temp_n)<p_small ) 
         WHERE ( pu_e_p< qpcMEZ(:,zoo) )
           limit = iiP
         END WHERE
     ELSEWHERE
         WHERE ( pu_e_n<qncMEZ(:,zoo) )
           limit = iiN
         END WHERE
     END WHERE

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Compute the correction terms depending on the limiting constituent
     ! Eq. 42 Vichi et al 2007 for a combination of N and P limitation
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     WHERE     ( limit == iiC )
         pe_R6c = ZERO
         pe_N1p = max(ZERO, (ONE - p_peI(zoo))*rut_p - p_qpcMEZ(zoo)*ru_c)
         pe_N4n = max(ZERO, (ONE - p_peI(zoo))*rut_n - p_qncMEZ(zoo)*ru_c)
     ELSEWHERE ( limit == iiP )
         pe_R6c = max(ZERO, ru_c - (ONE - p_peI(zoo))*rut_p/p_qpcMEZ(zoo))
         pe_N1p = ZERO
         pe_N4n = max( ZERO, (ONE - p_peI(zoo))*rut_n - p_qncMEZ(zoo)*(ru_c - pe_R6c))
     ELSEWHERE ( limit == iiN )
         pe_R6c = max(ZERO, ru_c - (ONE - p_peI(zoo))*rut_n/p_qncMEZ(zoo))
         pe_N1p = max(ZERO, (ONE - p_peI(zoo))*rut_p - p_qpcMEZ(zoo)*(ru_c - pe_R6c))
         pe_N4n = ZERO
     END WHERE
  else
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Eliminate the excess of the non-limiting constituent under fixed quota
     ! Determine release due to either C, P, or N limitation
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     call fixratio(tfluxC,tfluxN,tfluxP,qncMEZ(:,zoo),qpcMEZ(:,zoo),pe_R6c, pe_N4n, pe_N1p)
   
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correction term for excess of non-limiting nutrients as organic carbon 
  ! release (POC) and nutrient remineralization (PO4 and NH)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector(iiPel, ppzooc, ppR6c, pe_R6c)
  call flux_vector(iiPel, ppzoop, ppN1p, pe_N1p)
  call flux_vector(iiPel, ppzoon, ppN4n, pe_N4n)

  end subroutine MesoZooDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
