!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: MicroZoo
!
! DESCRIPTION
!  Microzooplankton dynamics
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
  subroutine MicroZooDynamics(zoo)
!
! USES
  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants,  ONLY:MW_C,C2ALK
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D3STATE, O2o, R1c, R6c, R1n, R6n, &
    R2c,R1p, R6p, N4n, N1p, PhytoPlankton, MicroZooPlankton, PelBacteria
  use mem, ONLY: ppPelBacteria, ppO2o, ppR1c, ppR6c, ppR6s, Depth,&
    ppR1n, ppR6n, ppR1p, ppR6p, ppN4n, ppN1p, ppPhytoPlankton, ppMicroZooPlankton, &
    ETW, eO2mO2, qncPBA, qpcPBA, qncPPY, qpcPPY, qncMIZ, qpcMIZ, iiPelBacteria, &
    qlcPPY, qscPPY, iiPhytoPlankton, iiMicroZooPlankton, iiC, iiN, iiP, iiL, iiS, &
    NO_BOXES, iiBen, iiPel, flux_vector, quota_flux
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c, ppO5c, ppO3h, qccPPY
#endif
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF, qfcPPY, ppR6f
#endif
#endif
  use mem_Param,     ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_small
  use bfm_error_msg, ONLY: bfm_error
  use mem_MicroZoo
  use mem_globalfun,   ONLY: eTq, MM, fixratio

  IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: zoo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer       :: i
  integer       :: ppzooc, ppzoon, ppzoop
  integer, save :: first =0
  integer       :: AllocStatus
  integer,dimension(NO_BOXES)  :: limit
  real(RLEN),allocatable,save,dimension(:) :: sut,et,eO2,rumc,  &
                                         rugc,rugn,rugp,runc,runn,runp, &
                                         rrsc,rrac,reac,rdc,rrtc,ruPBAc,ruPPYc,  &
                                         ruMIZc,rric,rr1c,rr6c,rr1p,rr1n, &
                                         rrip,rr6p,rep,rrin,zooc, tfluxC, tfluxN, tfluxP
  real(RLEN),allocatable,save,dimension(:)    :: rr6n,ren,pu_ra,r
  real(RLEN),allocatable,save,dimension(:)    :: pe_N1p, pe_N4n, pe_R6c
  real(RLEN),allocatable,save,dimension(:,:)  :: PBAc,PPYc,MIZc
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     ALLOCATE ( PBAc(NO_BOXES,iiPelBacteria),   PPYc(NO_BOXES,iiPhytoPlankton),  &
        &       MIZc(NO_BOXES,iiMicroZooPlankton),               &
        &       sut(NO_BOXES), et(NO_BOXES), eO2(NO_BOXES),      &
        &       rumc(NO_BOXES),  &
        &       rugc(NO_BOXES), rugn(NO_BOXES), rugp(NO_BOXES),  &
        &       runc(NO_BOXES), runn(NO_BOXES), runp(NO_BOXES),  &
        &       rrsc(NO_BOXES), rrac(NO_BOXES), reac(NO_BOXES),  &
        &       rdc(NO_BOXES) , rrtc(NO_BOXES), ruPBAc(NO_BOXES), ruPPYc(NO_BOXES), &
        &       ruMIZc(NO_BOXES), rric(NO_BOXES), rr1c(NO_BOXES), rr6c(NO_BOXES), &
        &       rr1p(NO_BOXES), rr1n(NO_BOXES), zooc(NO_BOXES), rrip(NO_BOXES), &
        &       rr6p(NO_BOXES), rep(NO_BOXES), rrin(NO_BOXES), rr6n(NO_BOXES), &
        &       ren(NO_BOXES), pu_ra(NO_BOXES), r(NO_BOXES), &
        &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES), &
        &       pe_N4n(NO_BOXES), pe_N1p(NO_BOXES), pe_R6c(NO_BOXES), &
        &      STAT = AllocStatus )
    
     IF( AllocStatus /= 0 ) call bfm_error('MicroZooDynamics','Error allocating arrays')
     first=1
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppzooc = ppMicroZooPlankton(zoo,iiC)
  ppzoon = ppMicroZooPlankton(zoo,iiN)
  ppzoop = ppMicroZooPlankton(zoo,iiP)
  zooc = D3STATE(:,ppzooc)

  ! Quota collectors
  tfluxC = ZERO
  tfluxN = ZERO
  tfluxP = ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et = eTq(ETW(:), p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2 = min(ONE, MM(O2o(:), p_chro(zoo)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total potential food given the non-dim prey availability
  ! and capture efficiency with loops over all LFGs.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rumc   = ZERO
  do i = 1 ,iiPelBacteria
     PBAc(:,i) = p_paPBA(zoo,i)*PelBacteria(i,iiC)* &
                 MM(PelBacteria(i,iiC), p_minfood(zoo))
     rumc = rumc + PBAc(:,i)
  end do

  do i = 1 ,iiPhytoPlankton
     PPYc(:,i) = p_paPPY(zoo,i)*PhytoPlankton(i,iiC)* &
                   MM(PhytoPlankton(i,iiC), p_minfood(zoo))
     rumc = rumc + PPYc(:,i)
  end do

  do i = 1, iiMicroZooPlankton
     MIZc(:,i) = p_paMIZ(zoo,i)*MicroZooPlankton(i,iiC)* &
                   MM(MicroZooPlankton(i,iiC), p_minfood(zoo))
     rumc = rumc + MIZc(:,i)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake rate (eq 38 Vichi et al. 2007) and 
  ! specific uptake rate considering potentially available food (sut)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  = et*p_sum(zoo)*MM(rumc, p_chuc(zoo))*zooc
  sut = rugc / (p_small + rumc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptakes from every LFG
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Bacterioplankton
  rugn = ZERO
  rugp = ZERO

  do i = 1, iiPelBacteria
    ruPBAc = sut*PBAc(:,i)
    call quota_flux(iiPel, ppzooc, ppPelBacteria(i,iiC), ppzooc, ruPBAc            , tfluxC)
    call quota_flux(iiPel, ppzoon, ppPelBacteria(i,iiN), ppzoon, ruPBAc*qncPBA(:,i), tfluxN)
    call quota_flux(iiPel, ppzoop, ppPelBacteria(i,iiP), ppzoop, ruPBAc*qpcPBA(:,i), tfluxP)
    rugn = rugn + ruPBAc*qncPBA(:,i)
    rugp = rugp + ruPBAc*qpcPBA(:,i)
  end do
  ! Phytoplankton
  do i = 1, iiPhytoPlankton
    ruPPYc = sut*PPYc(:,i)
    call quota_flux(iiPel, ppzooc, ppPhytoPlankton(i,iiC), ppzooc, ruPPYc            , tfluxC)
    call quota_flux(iiPel, ppzoon, ppPhytoPlankton(i,iiN), ppzoon, ruPPYc*qncPPY(:,i), tfluxN)
    call quota_flux(iiPel, ppzoop, ppPhytoPlankton(i,iiP), ppzoop, ruPPYc*qpcPPY(:,i), tfluxP)
    rugn = rugn + ruPPYc*qncPPY(:,i)
    rugp = rugp + ruPPYc*qpcPPY(:,i)
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
    ! Note that intra-group predation (cannibalism) is not added as a flux
    if ( i/= zoo) then
       call quota_flux(iiPel, ppzooc, ppMicroZooPlankton(i,iiC), ppzooc, ruMIZc            , tfluxC)
       call quota_flux(iiPel, ppzoon, ppMicroZooPlankton(i,iiN), ppzoon, ruMIZc*qncMIZ(:,i), tfluxN)
       call quota_flux(iiPel, ppzoop, ppMicroZooPlankton(i,iiP), ppzoop, ruMIZc*qpcMIZ(:,i), tfluxP)
    end if
    rugn = rugn + ruMIZc*qncMIZ(:,i)
    rugp = rugp + ruMIZc*qpcMIZ(:,i)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes from microzooplankton
  ! The metabolic balance is the following:
  ! Ingestion = Growth + Excretion + Respiration
  ! Assimilation efficiency p_pu = G/I
  ! Excretion E = I*p_pu_ea
  ! therefore R = (1-p_pu-p_pu_ea)*I
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Rest, activity and total respiration fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrsc = p_srs(zoo)*et*zooc
  ! the activity respiration is derived from the other constant parameters
  rrac = rugc*(ONE - p_pu(zoo) - p_pu_ea(zoo))
  rrtc = rrsc + rrac
  call quota_flux(iiPel, ppzooc, ppzooc, ppO3c, rrtc, tfluxC)
  call flux_vector(iiPel, ppO2o, ppO2o, -rrtc/MW_C)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality (rdc) + Activity Excretion (reac)
  ! and partitioning between particulate and dissolved
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rdc  = ((ONE - eO2)*p_sdo(zoo) + p_sd(zoo))*zooc
  reac = rugc* p_pu_ea(zoo)
  rric = reac + rdc
  rr1c = rric*p_pe_R1c
  rr6c = rric*(ONE - p_pe_R1c)
  call quota_flux(iiPel, ppzooc, ppzooc, ppR1c, rr1c, tfluxC)
  call quota_flux(iiPel, ppzooc, ppzooc, ppR6c, rr6c, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !     Nutrient dynamics in microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrin = rugn*p_pu_ea(zoo) + rdc*qncMIZ(:,zoo)
  rr1n = rrin*p_pe_R1n
  rr6n = rrin - rr1n
  call quota_flux(iiPel, ppzoon, ppzoon, ppR1n, rr1n, tfluxN)
  call quota_flux(iiPel, ppzoon, ppzoon, ppR6n, rr6n, tfluxN)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Phosphorus dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrip = rugp*p_pu_ea(zoo) + rdc*qpcMIZ(:,zoo)
  rr1p = rrip*p_pe_R1p
  rr6p = rrip - rr1p
  call quota_flux(iiPel, ppzoop, ppzoop, ppR1p, rr1p, tfluxP)
  call quota_flux(iiPel, ppzoop, ppzoop, ppR6p, rr6p, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved nutrient dynamics
  ! Compare the quota of the net growth rates with the optimal quota
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  runc = max(ZERO, rugc*(ONE - p_pu_ea(zoo)) - rrac)
  runn = max(ZERO, rugn*(ONE - p_pu_ea(zoo)) + rrsc*qncMIZ(:,zoo))
  runp = max(ZERO, rugp*(ONE - p_pu_ea(zoo)) + rrsc*qpcMIZ(:,zoo))
  ren  = max(ZERO, runn/(p_small + runc) - p_qncMIZ(zoo))* runc
  rep  = max(ZERO, runp/(p_small + runc) - p_qpcMIZ(zoo))* runc
  call quota_flux(iiPel, ppzoon, ppzoon, ppN4n, ren, tfluxN)
  call quota_flux(iiPel, ppzoop, ppzoop, ppN1p, rep, tfluxP)

  if ( ppzoon == 0 .or. ppzoop == 0 ) then

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Eliminate the excess of the non-limiting constituent under fixed quota
     ! Determine release due to either C, P, or N limitation
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     call fixratio(tfluxC,tfluxN,tfluxP,qncMIZ(:,zoo),qpcMIZ(:,zoo),pe_R6c, pe_N4n, pe_N1p)

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Apply the correction terms depending on the limiting constituent
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     call flux_vector(iiPel, ppzooc, ppR6c, pe_R6c*(ONE-p_pe_R1c) )
     call flux_vector(iiPel, ppzooc, ppR1c, pe_R6c*(p_pe_R1c))
     call flux_vector(iiPel, ppzoop, ppN1p, pe_N1p)
     call flux_vector(iiPel, ppzoon, ppN4n, pe_N4n)

  endif

  end subroutine MicroZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
