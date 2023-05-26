!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PhytoDynamics
!
! DESCRIPTION
!   This process describes the dynamics of all phytoplankton
!   groups. The differences in behaviour are expressed
!   by differences in parameter-values only.
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
  subroutine PhytoDynamics(phyto)
!
! USES
  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants, ONLY: MW_C,C2ALK
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: iiC,iiN,iiP,iiS,iiL
  use mem, ONLY: D3STATE, R1c, R6c, O2o, R2c, &
                 N3n, N4n, N1p, R1n, R6n, R1p, R6p, N5s
  use mem, ONLY: ppR1c, ppR6c, ppO2o, ppR2c, ppN3n, ppN4n, ppN1p, ppR1n, &
    ppR6n, ppR1p, ppR6p, ppN5s, ppR6s, SUNQ, ThereIsLight, ETW, EIR, &
    xEPS, Depth, eiPPY, sediPPY, sunPPY, qpcPPY, qncPPY, qscPPY, qlcPPY, NO_BOXES, &
    iiBen, iiPel, flux_vector, quota_flux, ppPhytoPlankton
#ifdef BFM_OGS
  use mem, ONLY: exPPYR2ac
#endif
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c, ppO5c, ppO3h, qccPPY
#endif
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF,N7f,qfcPPY,ppN7f,ppR6f,ppR1f
#endif
#endif
  use constants,  ONLY: SEC_PER_DAY, E2W, HOURS_PER_DAY
  use mem_Param,  ONLY: p_small, ChlDynamicsFlag
  use mem_PAR,    ONLY: LightPeriodFlag, LightLocationFlag
  use mem_Phyto
  use mem_globalfun,   ONLY: eTq, MM, insw

  IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: phyto

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, save :: first=0
  integer :: ppphytoc, ppphyton, ppphytop, ppphytos, ppphytol 
  integer,dimension(NO_BOXES)  :: limit
  real(RLEN),allocatable,save,dimension(:) :: phytoc,phyton,phytop,phytos,phytol
                                                                                                                                                             
  real(RLEN),allocatable,save,dimension(:) :: r,tmp,et,sum,sadap,sea,sdo,rugc,  &
        & srt,slc,run,pe_R6,rupp,rump,misp,rupn,rumn3,rumn4,rumn,netgrowth,     &
        & misn,cqun3,sra,srs,rums,rups,miss,tN,fpplim,iN,iN1p,rr1n,rr1p,rr6c,   &
        & rr6n,rr6p,rr6s,runn,runn3,runn4,runp,runs,Irr,rho_Chl,rate_Chl,seo,   &
        & flPIR2c,iNIn,eN5s,rrc,rr1c,iN5s,chl_opt

  real(RLEN),allocatable,save,dimension(:) :: tfluxC, tfluxN, tfluxP, pe_R1n, pe_R1p, pe_R1c
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
#ifdef INCLUDE_PELFE
  integer :: ppphytof
  real(RLEN),allocatable,save,dimension(:) :: phytof
  real(RLEN),allocatable,save,dimension(:) :: iN7f,misf,rr1f,rr6f,rupf,rumf,runf
#endif
  integer :: AllocStatus, DeallocStatus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     ALLOCATE ( phytoc(NO_BOXES), phyton(NO_BOXES), phytop(NO_BOXES),           &
        &       phytos(NO_BOXES), phytol(NO_BOXES), r(NO_BOXES), tmp(NO_BOXES), &
        &       et(NO_BOXES), sum(NO_BOXES), sadap(NO_BOXES),                   &
        &       sdo(NO_BOXES), rugc(NO_BOXES), sra(NO_BOXES), srs(NO_BOXES),    &
        &       srt(NO_BOXES), slc(NO_BOXES), run(NO_BOXES), pe_R6(NO_BOXES),   &
        &       rupp(NO_BOXES), rump(NO_BOXES), misp(NO_BOXES), rupn(NO_BOXES), & 
        &       rumn3(NO_BOXES), rumn4(NO_BOXES), rumn(NO_BOXES),               & 
        &       netgrowth(NO_BOXES), misn(NO_BOXES), cqun3(NO_BOXES),           &
        &       rums(NO_BOXES), rups(NO_BOXES), miss(NO_BOXES), tN(NO_BOXES),   &
        &       fpplim(NO_BOXES), iN(NO_BOXES), iN1p(NO_BOXES), iNIn(NO_BOXES), &
        &       eN5s(NO_BOXES), iN5s(NO_BOXES), rrc(NO_BOXES),                  &
        &       rr1c(NO_BOXES), rr1n(NO_BOXES), rr1p(NO_BOXES),                 &
        &       rr6c(NO_BOXES), rr6n(NO_BOXES), rr6p(NO_BOXES), rr6s(NO_BOXES), &
        &       runn(NO_BOXES), runn3(NO_BOXES), runn4(NO_BOXES),               &
        &       runp(NO_BOXES), runs(NO_BOXES), Irr(NO_BOXES),                  &
        &       flPIR2c(NO_BOXES), seo(NO_BOXES), sea(NO_BOXES),                &
        &       rate_Chl(NO_BOXES), rho_Chl(NO_BOXES), chl_opt(NO_BOXES),       & 
        &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES),           &
        &       pe_R1n(NO_BOXES), pe_R1p(NO_BOXES), pe_R1c(NO_BOXES),           &
#ifdef INCLUDE_PELFE
        &       phytof(NO_BOXES), misf(NO_BOXES), iN7f(NO_BOXES),               &
        &       rr1f(NO_BOXES), rr6f(NO_BOXES),                                 &
        &       rupf(NO_BOXES), rumf(NO_BOXES), runf(NO_BOXES),                 &
#endif
        &      STAT = AllocStatus )
     IF( AllocStatus /= 0 ) call bfm_error('PhytoDynamics','Error allocating arrays')
     first=1
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppphytoc = ppPhytoPlankton(phyto,iiC)
  ppphyton = ppPhytoPlankton(phyto,iiN)
  ppphytop = ppPhytoPlankton(phyto,iiP)
  ppphytos = ppPhytoPlankton(phyto,iiS)
  ppphytol = ppPhytoPlankton(phyto,iiL)
  phytoc(:) = D3STATE(:,ppphytoc)
  phyton(:) = phytoc(:) * qncPPY(:,phyto)
  phytop(:) = phytoc(:) * qpcPPY(:,phyto)
  phytol(:) = phytoc(:) * qlcPPY(:,phyto)
  if ( ppphytos > 0 )  phytos(:) = phytoc(:) * qscPPY(:,phyto)
#ifdef INCLUDE_PELFE
  ppphytof = ppPhytoPlankton(phyto,iiF)
  if ( ppphytof > 0 ) phytof(:) = phytoc(:) * qfcPPY(:,phyto)
#endif
  ! Quota collectors
  tfluxC = ZERO
  tfluxN = ZERO
  tfluxP = ZERO
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitations (intracellular and extracellular)
  ! fpplim is the combined non-dimensional factor limiting photosynthesis
  ! Note for silicate limitation:
  !  p_switchSi =1 : external regulation of silica limitation 
  !  p_switchSi =2 : internal regulation of silica limitation 
  ! The standard Michaelis-Menten formulation contains the Contois parameter
  ! p_Contois=0: standard Michaelis Menten Formulation
  ! 0<p_Contois<=1: The Contois formulation is active. 
  !                 The limiting role of the population size (intraspecific 
  !                 competition) can be tuned by increasing p_Contois 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iN1p = min( ONE, max( p_small, ( qpcPPY(:,phyto) &
         - p_qplc(phyto))/( p_qpcPPY(phyto)- p_qplc(phyto))))
  iNIn = min( ONE, max( p_small, ( qncPPY(:,phyto) &
         - p_qnlc(phyto))/( p_qncPPY(phyto)- p_qnlc(phyto))))
  if (ppphytos > 0) then
     select case (p_switchSi(phyto)) 
       case (1)  ! external control
         eN5s = min( ONE, N5s(:)/(N5s(:) + p_chPs(phyto)+(p_Contois(phyto)*phytos(:))))
         fpplim = eN5s
         iN5s   = ONE  
       case (2) ! internal control
         iN5s = min(ONE, max( p_small, ( qscPPY(:,phyto) &
                - p_qslc(phyto))/( p_qscPPY(phyto)- p_qslc(phyto))))
         fpplim = iN5s
         eN5s   = ONE  
     end select
  else 
     iN5s   = ONE  
     eN5s   = ONE  
     fpplim = ONE
  end if
#ifdef INCLUDE_PELFE
  if (ppphytof > 0) then
     iN7f = min( ONE, max( p_small, ( qfcPPY(:,phyto) &
            - p_qflc(phyto))/( p_qfcPPY(phyto)- p_qflc(phyto))))
  end if
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Multiple nutrient limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_limnut(phyto))
    case ( 0 )
      iN  =   (iN1p* iNIn)**(0.5_RLEN)  ! geometric mean

    case ( 1 )
#ifdef INCLUDE_PELFE
      if (ppphytof > 0) then
         iN  =   min(  iN1p,  iNIn, iN7f)  ! Liebig rule
      else
         iN  =   min(  iN1p,  iNIn)  ! Liebig rule
      endif
#else
      iN  =   min(  iN1p,  iNIn)  ! Liebig rule
#endif

    case ( 2 )
      iN  =   2.0_RLEN/( ONE/ iN1p+ ONE/ iNIn)  ! combined

  end select

  ! tN only controls sedimentation of phytoplankton (Liebig)
  tN= min(iN,fpplim)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Phytoplankton
  ! Include cut-off at low temperature if p_temp>0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq(  ETW(:),  p_q10(phyto))
  et  =   max(ZERO,et-p_temp(phyto))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Light limitation with Chl dynamics
  ! If Chl is a diagnostic variable the limiting factor has been 
  ! computed in Light/PhotoAvailableRadiation.F90
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( ChlDynamicsFlag== 2) then
    ! Irradiance EIR is in uE m-2 s-1, 
    ! Irr is top, middle or average irradiance in uE m-2 day-1
    select case ( LightLocationFlag)
       case ( 1 )
          ! Light at the top of the cell
          Irr =  max(p_small,EIR(:))*SEC_PER_DAY;
       case ( 2 )
          ! Light in the middle of the cell
          Irr = max(p_small,EIR(:))*exp(-xEPS(:)*0.5_RLEN*Depth(:)) &
                *SEC_PER_DAY
       case ( 3 ) ! default
          ! Average Light in the cell
          r = xEPS(:)* Depth(:)
          r = EIR(:)/xEPS(:)/Depth(:)*(ONE-exp(-r))
          Irr = max(p_small,r*SEC_PER_DAY)
    end select
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Compute exponent E_PAR/E_K = alpha0/PBmax
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    r(:) = qlcPPY(:,phyto)*p_alpha_chl(phyto)/p_sum(phyto)* Irr
    select case ( LightPeriodFlag)
      case ( 1 ) ! instantaneous light
        ! no other factors needed
      case ( 2 ) ! daylight average is used
        ! recompute r and photsynthesis limitation using daylight scaling
        fpplim  =   fpplim*SUNQ(:)/HOURS_PER_DAY
        r(:) = r(:)*HOURS_PER_DAY/SUNQ(:)
      case ( 3 ) ! on-off
        fpplim  =   fpplim*ThereIsLight
    end select
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Light limitation factor according to Platt
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    eiPPY(:,phyto) = ( ONE- exp( - r))
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total photosynthesis
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sum  =   p_sum(phyto)*et*eiPPY(:,phyto)*fpplim

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration rate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sra  =   p_pu_ra(phyto)* sum ! activity
  srs  =   et* p_srs(phyto)                   ! basal
  srt  =   sra+ srs                           ! total
  rrc  =   srt* phytoc                        ! total actual respiration

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Stress Lysis
  sdo = ( p_thdo(phyto)/( iN+ p_thdo(phyto)))* p_sdmo(phyto)
  ! extra lysis for high-density
  sdo = sdo+ p_seo(phyto)* MM(phytoc, p_sheo(phyto))
  !
  ! Activity Excretion
  sea = sum * p_pu_ea(phyto)
  !
  ! Nutrient Stress Excretion
  if (p_netgrowth(phyto)) then
     seo = ZERO
  else 
     seo = sum*(ONE-p_pu_ea(phyto)-p_pu_ra(phyto))*(ONE- iN) 
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apportioning over R1 and R6:
  ! Cell lysis generates both DOM and POM.
  ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
  ! Assuming that this structural part is not easily degradable,
  ! at least a fraction equal to the minimum quota is released as POM.
  ! Therefore, nutrients (and C) in the structural part go to R6.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pe_R6 = min( p_qplc(phyto)/( qpcPPY(:,phyto)+ p_small), p_qnlc(phyto)/ &
          ( qncPPY(:,phyto)+ p_small))
  pe_R6 = min(  ONE,  pe_R6)
  rr6c  =     pe_R6     * sdo * phytoc
  rr1c  = (ONE - pe_R6) * sdo * phytoc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Production, productivity and C flows
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Gross Production
  rugc = sum* phytoc

  ! Specific loss terms
  slc = sea + seo + srt + sdo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Exudation of carbohydrate toward DOC
  ! 1 - to R1 (Vichi2007b)
  ! 2 - to R2 (Vichi2004). If p_netgrowth then seo = ZERO
  ! 3 - Activity to R1 and Nut-Stress to R2 (BFM Manual)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case (p_switchDOC(phyto))
    case (1)
       ! All activity excretions are assigned to R1
       rr1c = rr1c + sea*phytoc + seo*phytoc
       flPIR2c = ZERO
    case (2)
       ! Activity and Nutrient-stress excretions are assigned to R2
       flPIR2c = seo*phytoc + sea*phytoc
    case (3)
       ! Activity excretion is only assigned to R2
       rr1c = rr1c + sea*phytoc
       flPIR2c = seo*phytoc
  end select

  call quota_flux( iiPel, ppphytoc ,ppO3c,ppphytoc, rugc, tfluxC )  
  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR1c, rr1c, tfluxC )
  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR6c, rr6c, tfluxC )

  call quota_flux( iiPel, ppphytoc, ppphytoc,ppO3c, rrc, tfluxC )
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrc/ MW_C) )
  call flux_vector( iiPel, ppO2o,ppO2o, rugc/ MW_C ) 

#ifdef BFM_OGS
  exPPYR2ac(:,phyto) = sea* phytoc
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential-Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (p_netgrowth(phyto)) then
     sadap  =   max(  0.05_RLEN,  sum- slc)
  else
     sadap  =   et*p_sum(phyto)
  end if
  run  =   max(  ZERO, ( sum- (sea + seo + srt))* phytoc)  ! net production

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Uptake: calculate maximum uptake of N, P
  ! based on affinity
  ! Ammonium preference is considered if p_lN4 /= 0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cqun3  =   p_lN4(phyto)/( p_small + p_lN4(phyto)+ N4n(:))
  rumn3  =   p_qun(phyto)* N3n(:)* phytoc* cqun3  ! max pot. uptake of N3
  rumn4  =   p_qun(phyto)* N4n(:)* phytoc  ! max pot. uptake of N4
  rumn  =   rumn3+ rumn4  ! max pot. uptake of DIN
  rump  =   p_qup(phyto)* N1p(:)* phytoc  ! max pot. uptake of PO4

  if (p_netgrowth(phyto)) then
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Check which fraction of fixed C can be used for new biomass
   ! given the internal storage.
   ! N and P uptake are compared sequentially
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth = min( run, ( rumn+ max( ZERO, 0.05_RLEN* &
      rugc*( qncPPY(:,phyto)- p_qnlc(phyto))))/ p_qnlc(phyto))
      netgrowth = min( netgrowth, ( rump+ max( ZERO, &
       0.05_RLEN* rugc*( qpcPPY(:,phyto)- p_qplc(phyto))))/ p_qplc(phyto))
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Excrete C that cannot be used for growth as carbo-hydrates:
   ! Correct the net C uptake
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth  =   max(  netgrowth,  ZERO)
      flPIR2c  =   flPIR2c+ run- netgrowth
      run  =   netgrowth
  end if

  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR2c, flPIR2c, tfluxC )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Specific net growth rate (d-1)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sunPPY(:,phyto)  =   run/( p_small+ phytoc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misn  =   sadap*( p_xqn(phyto)* p_qncPPY(phyto)* phytoc- phyton)  ! Intracellular missing amount of N
  rupn  =   p_xqn(phyto)* p_qncPPY(phyto)* run  ! N uptake based on net assimilat. C
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of NI

  r  =   insw(runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of Nn
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of Nn
  call quota_flux( iiPel, ppphyton, ppN3n,ppphyton, runn3, tfluxN )  ! source/sink.n
  call quota_flux( iiPel, ppphyton, ppN4n,ppphyton, runn4, tfluxN )  ! source/sink.n
  tmp = - runn*( ONE- r)
  call quota_flux( iiPel, ppphyton, ppphyton,ppR1n,tmp, tfluxN)  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misp  =   sadap*( p_xqp(phyto)* p_qpcPPY(phyto)* phytoc- phytop)  ! intracellular missing amount of P
  rupp  =   p_xqp(phyto)* run* p_qpcPPY(phyto)  ! P uptake based on C uptake
#ifdef EXTRACOST
  rupp  =   p_xqp(phyto)* run* p_qpcPPY(phyto)-( sdo+ srs)* phytop  ! P uptake based on C uptake
#endif
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

  r  =   insw(runp)
  tmp = runp*r
  call quota_flux(iiPel, ppphytop, ppN1p,ppphytop, tmp, tfluxP)  ! source/sink.p
  tmp = - runp*( ONE- r)
  call quota_flux(iiPel, ppphytop, ppphytop,ppR1p, tmp, tfluxP)  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr6n =     pe_R6     * sdo * phyton
  rr1n = (ONE - pe_R6) * sdo * phyton

  rr6p =     pe_R6     * sdo * phytop
  rr1p = (ONE - pe_R6) * sdo * phytop

  call quota_flux( iiPel, ppphyton, ppphyton,ppR6n, rr6n, tfluxN )  ! source/sink.n
  call quota_flux( iiPel, ppphyton, ppphyton,ppR1n, rr1n, tfluxN )  ! source/sink.n

  call quota_flux( iiPel, ppphytop, ppphytop,ppR6p, rr6p, tfluxP )  ! source/sink.p
  call quota_flux( iiPel, ppphytop, ppphytop,ppR1p, rr1p, tfluxP )  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( ppphytos > 0 )  then
    select case (p_switchSi(phyto))
    case (1)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Gross uptake of silicate excluding respiratory costs
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      miss = max(ZERO, phytos - p_qscPPY(phyto)*phytoc) ! intracellular missing Si
      runs = p_qscPPY(phyto) * (sum-sra-sea-seo) * phytoc
      rr6s = (srs+sdo) * phytos + miss
    case (2)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Silicate uptake based on intracellular needs (note, no luxury)
      !  There can be efflux of dissolved silicate (M-J et al., 2000)
      !  however this generates fake remineralization and it is not implemented
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rums  =   p_qus(phyto)* N5s(:)* phytoc  ! max pot uptake based on affinity
      miss  =   sadap*(p_qscPPY(phyto)*phytoc - phytos) ! intracellular missing Si
      rups  =   run* p_qscPPY(phyto)  ! Si uptake based on net C uptake
      runs  =   min(  rums,  rups+ miss)  ! actual uptake
      rr6s  =   sdo*phytos
    end select
              
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Uptake and Losses of Si (only lysis)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    r  =   insw(runs)
    call flux_vector( iiPel, ppN5s, ppphytos, runs*r)
    call flux_vector( iiPel, ppphytos, ppR6s, rr6s - runs*(ONE-r))
  endif

#ifdef INCLUDE_PELFE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: IRON
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (ppphytof > 0) then
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Net uptake
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     rumf  =   p_quf(phyto)* N7f(:)* phytoc  ! max potential uptake
     ! intracellular missing amount of Fe
     misf  =   sadap*(p_xqf(phyto)*p_qfcPPY(phyto)*phytoc - phytof)  
     rupf  =   p_xqf(phyto)* run* p_qfcPPY(phyto)  ! Fe uptake based on C uptake
     runf  =   min(  rumf,  rupf+ misf)  ! actual uptake
     r  =   insw(runf)
     ! uptake from inorganic if shortage
     call flux_vector( iiPel, ppN7f,ppphytof, runf* r )
     ! release to dissolved organic to keep the balance if excess
     call flux_vector(iiPel, ppphytof,ppR1f,- runf*( ONE- r))
   
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Losses of Fe
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     rr6f  =   rr6c* qfcPPY(:,phyto)
     rr1f  =   sdo* phytof- rr6f
     call flux_vector( iiPel, ppphytof,ppR1f, rr1f )
     call flux_vector( iiPel, ppphytof,ppR6f, rr6f )
  end if
#endif

  if ( ChlDynamicsFlag== 2) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Chl-a synthesis and photoacclimation
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case (p_switchChl(phyto))
      case (1) ! PELAGOS
           rho_Chl = p_qlcPPY( phyto)* min(ONE, p_sum(phyto)* eiPPY(:,phyto)* phytoc/( &
                     p_alpha_chl(phyto)*( phytol+ p_small)* Irr))
           rate_Chl = rho_Chl*(sum - seo - sea - sra) * phytoc - sdo*phytol
      case (2) ! OPATM-BFM
           rho_Chl  =   p_qlcPPY(phyto)* sum/( p_alpha_chl(phyto)* qlcPPY(:,phyto)* Irr)
           rate_Chl = iN* rho_Chl* run- max( p_sdchl(phyto)*( ONE - iN), sdo)* &
               phytol+ min( ZERO, sum- slc+ sdo)* max( ZERO, phytol- p_qlcPPY(phyto)* phytoc)
      case (3) ! UNIBO
           rho_Chl = p_qlcPPY(phyto)*min(ONE,          &
                     (sum-seo-sea-sra) *phytoc /          &
                     (p_alpha_chl(phyto)*(phytol+p_small) *Irr))
           ! The "optimal" chl concentration corresponds to the chl that
           ! (given the actual C biomass) would give (Epar/Ek)=p_EpEk
           chl_opt = p_EpEk_or(phyto)*p_sum(phyto)*phytoc/  &
                     (p_alpha_chl(phyto)*Irr+p_small)
           !  Actual chlorophyll concentration exceeding the "optimal" value is 
           !  discarded with a p_tochl_relt relaxation.
           rate_Chl = rho_Chl*(sum-seo-sea-sra)*phytoc-(sdo+srs)*phytol - &
                      max(ZERO,(phytol-chl_opt))*p_tochl_relt(phyto)
      case (4) ! NIOZ
          ! total synthesis, only when there is net production (run > 0)
          ! The fixed loss rate due to basal respiration is introduced to have 
          ! chl loss in the absence of light (< 1 uE/m2/s)
           rho_Chl = p_qlcPPY( phyto)* min(ONE, p_sum(phyto)* eiPPY(:,phyto)* phytoc/( &
                     p_alpha_chl(phyto)*( phytol+ p_small)* Irr))
           rate_Chl = rho_Chl*run - p_sdchl(phyto)*phytol*max( ZERO, ( p_thdo(phyto)-tN)) &
                     -srs * phytol * ONE/(Irr+ONE)
    end select
    call flux_vector( iiPel, ppphytol,ppphytol, rate_Chl )
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Sedimentation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( p_res(phyto)> ZERO) then
    sediPPY(:,phyto) = sediPPY(:,phyto) &
                   + p_res(phyto)* max( ZERO, ( p_esNI(phyto)-tN))
  end if

#if defined INCLUDE_PELCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! PIC (calcite/aragonite) production
  ! The idea in PISCES is that the calcite flux exists only when associated
  ! to a carbon release from phytoplankton (there is no calcite storage in phyto)
  ! First compute the realized rain ratio for each phytoplankton species
  ! The presence of PIC in phytoplankton group is controlled by p_caco3r
  ! with the following regulating factors:
  !  - nutrient limitation
  !  - temperature enhancement
  !  - density enhancement
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( p_caco3r(phyto) > ZERO ) then
     qccPPY(:,phyto) = p_caco3r(phyto)
     ! Calcite production represented as a flux between DIC and PIC, impacting ALK
     call flux_vector( iiPel, ppO3c,ppO5c, qccPPY(:,phyto)*rr6c )
     call flux_vector( iiPel, ppO3h,ppO3h, -C2ALK*qccPPY(:,phyto)*rr6c )
  endif
#endif

  end subroutine PhytoDynamics
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
