!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PelBac
!
! DESCRIPTION
!   This process describes the dynamics of bacterioplankton
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
  subroutine PelBacDynamics(bac)
!
! USES
  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY:  R6c, R6n, R6p, R1c, R1n, R1p, R2c, O2o, N6r, &
    N4n, N1p, N3n, R3c, iiR1, iiR6, D3STATE
  use mem, ONLY: iiPelBacteria, ppPelBacteria, iiC, iiN, iiP, ppR6c, &
    ppR6n, ppR6p, ppR1c, ppR1n, ppR1p, &
    ppR2c, ppO2o, ppN6r, ppN4n, ppN1p, ppN3n, ppR3c, flPTN6r, Depth, ETW, &
    qncPBA, qpcPBA, qpcOMT, qncOMT, NO_BOXES, iiBen, iiPel, flux_vector,quota_flux
#ifdef BFM_OGS
    use mem, ONLY: BAC_ACT_FACT
#endif
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c
#endif
#endif
  use constants,  ONLY: MW_C, ONE_PER_DAY
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro, p_small
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
  real(RLEN),allocatable,save,dimension(:) :: et,eO2,r,flN6rPBA,rrc,  &
                                          rd,ruR1c,ruR1n,ruR1p,ruR2c,ruR3c,  &
                                          ruR6c,ruR6p,ruR6n,rump,  &
                                          rumn,rumn3,rumn4,ren,rep,reR2c, &
                                          reR3c,rut,rum,run,sun,rug, &
                                          suR2,cuR6,cuR1,iN1p,iNIn,iN, &
                                          eN1p,eN4n,huln, hulp, bacc, &
                                          tfluxC, tfluxN, tfluxP, pe_N4n, &
                                          pe_N1p, pe_R6c
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     ALLOCATE ( et(NO_BOXES), eO2(NO_BOXES), r(NO_BOXES),             &
        &       flN6rPBA(NO_BOXES), rrc(NO_BOXES), rd(NO_BOXES),      &
        &       ruR1c(NO_BOXES), ruR1n(NO_BOXES), ruR1p(NO_BOXES),    &
        &       ruR6c(NO_BOXES), ruR6p(NO_BOXES), ruR6n(NO_BOXES),    &
        &       ruR2c(NO_BOXES), ruR3c(NO_BOXES),                     &
        &       rump(NO_BOXES), rumn(NO_BOXES),                       &
        &       rumn3(NO_BOXES), rumn4(NO_BOXES), ren(NO_BOXES),      &
        &       rep(NO_BOXES), reR2c(NO_BOXES), reR3c(NO_BOXES),      &
        &       rut(NO_BOXES), rum(NO_BOXES), run(NO_BOXES),          &
        &       sun(NO_BOXES), rug(NO_BOXES), suR2(NO_BOXES),         &
        &       cuR6(NO_BOXES), cuR1(NO_BOXES),                       &
        &       iN1p(NO_BOXES), iNIn(NO_BOXES), iN(NO_BOXES),         &
        &       eN1p(NO_BOXES), eN4n(NO_BOXES), huln(NO_BOXES),       &
        &       hulp(NO_BOXES), bacc(NO_BOXES),                       &
        &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES), &
        &       pe_N4n(NO_BOXES), pe_N1p(NO_BOXES), pe_R6c(NO_BOXES), &
        &      STAT = AllocStatus )
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
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et = eTq(ETW(:), p_q10(bac))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Oxygen environment: bacteria are both aerobic and anaerobic
  ! To provide a faster switching between the two metabolic pathways the
  ! oxygen regulating factor eO2 is cubic
  ! (eq. 19 in Vichi et al., 2004)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2 = MM_power(max(p_small,O2o(:)),  p_chdo(bac),3)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! External nutrient limitation (used by some parametrizations)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eN4n = MM(N4n(:), p_chn(bac))
  eN1p = MM(N1p(:), p_chp(bac))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Mortality:
  !   1. first order mortality: p_sd 
  !   2. density dependent mortality due to virus infection: p_sd2
  !
  !   It is assumed that mortality is distributed in the same way over
  !   DOC (R1) and detritus (R6) as for phytoplankton and microzooplankton
  !   using the p_pe_R1x parameters defined in Param
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd  =  ( p_sd(bac)*et + p_sd2(bac)*bacc ) * bacc

  call quota_flux(iiPel, ppbacc, ppbacc, ppR6c, rd*(ONE-p_pe_R1c)              , tfluxC)
  call quota_flux(iiPel, ppbacn, ppbacn, ppR6n, rd*qncPBA(:,bac)*(ONE-p_pe_R1n), tfluxN)
  call quota_flux(iiPel, ppbacp, ppbacp, ppR6p, rd*qpcPBA(:,bac)*(ONE-p_pe_R1p), tfluxP)

  call quota_flux(iiPel, ppbacc, ppbacc, ppR1c, rd*p_pe_R1c              , tfluxC)
  call quota_flux(iiPel, ppbacn, ppbacn, ppR1n, rd*qncPBA(:,bac)*p_pe_R1n, tfluxN) 
  call quota_flux(iiPel, ppbacp, ppbacp, ppR1p, rd*qpcPBA(:,bac)*p_pe_R1p, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Substrate availability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_version(bac) )

    case ( BACT3 )  ! Polimene et al. (2006)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      ! Note: oxygen control in eq 5 (Polimene et al. 2006) is not included
      !        as bacteria are both aerobic and anaerobic
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  =  p_sum(bac)*et*bacc

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! No correction of organic material quality
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      cuR1 = ONE
      cuR6 = ONE

    case ( BACT1,BACT2 )  ! Vichi et al. (2004,2007), Lazzari et al. (2012) 

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
      ! Nutrient limitation (intracellular, eq. 51 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      iNIn = min(ONE, max(ZERO, qncPBA(:,bac)/p_qncPBA(bac)))  !Nitrogen
      iN1p = min(ONE, max(ZERO, qpcPBA(:,bac)/p_qpcPBA(bac)))  !Phosphorus
      iN   = min(iN1p, iNIn)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria (eq. 50 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  = iN*et*p_sum(bac)*bacc

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! correction of substrate quality depending on nutrient content
      ! (eq. 52 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      cuR1 = min(ONE, qpcOMT(:,iiR1)/p_qpcPBA(bac), qncOMT(:,iiR1)/ p_qncPBA(bac))
      cuR6 = min(ONE, qpcOMT(:,iiR6)/p_qpcPBA(bac), qncOMT(:,iiR6)/ p_qncPBA(bac))

  end select
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate the realized substrate uptake rate depending on the
  ! type of detritus and quality (cuRx)
  ! See eq 27 in Vichi et al., 2004 for R2
  ! and eq 6 in Polimene et al., 2006 for R3 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ruR1c = (p_suhR1(bac)*cuR1 + p_sulR1(bac)*(ONE-cuR1))*R1c(:)
  ruR2c = p_suR2(bac)*R2c(:)
  ruR3c = p_suR3(bac)*R3c(:)
  ruR6c = p_suR6(bac)*cuR6*R6c(:)
  rut   = p_small + ruR1c + ruR2c + ruR6c + ruR3c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria (eq. 50 Vichi et al. 2007)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rug = min( rum, rut )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR1c = rug*ruR1c/rut
  ruR2c = rug*ruR2c/rut
  ruR3c = rug*ruR3c/rut
  ruR6c = rug*ruR6c/rut

  call quota_flux(iiPel, ppbacc, ppR1c, ppbacc, ruR1c, tfluxC)
  call quota_flux(iiPel, ppbacc, ppR2c, ppbacc, ruR2c, tfluxC)
  call quota_flux(iiPel, ppbacc, ppR3c, ppbacc, ruR3c, tfluxC)
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
  ! Aerobic and anaerobic respiration
  ! Anaerobic bacteria use N6r as electron acceptor, with additional carbon cost
  ! If nitrate is present, the rate of consumption of N6r is converted to N3n
  ! consumption (eq 19 Vichi et al., 2004 and PelChem.F90)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity + Basal respiration
  rrc = (p_pu_ra(bac) * rug) + p_srs(bac)* bacc* et
  ! Aerobic consumption of Oxygen
  call flux_vector( iiPel, ppO2o, ppO2o, -eO2*rrc/MW_C )
  ! Anaerobic consumption of N6r
  flN6rPBA =( ((ONE-eO2)*rrc) + (p_pu_ra_o(bac)*(ONE-eO2)*rug) )/ MW_C* p_qro
  call flux_vector( iiPel, ppN6r, ppN6r, flN6rPBA )
  ! Bookkeeping of reduction equivalent formation rate
  flPTN6r(:) = flPTN6r(:) + flN6rPBA
  ! Total Respiration
  rrc = rrc + (p_pu_ra_o(bac)*(ONE-eO2)*rug)
  call quota_flux( iiPel, ppbacc, ppbacc, ppO3c, rrc, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Net Production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  run = rug - rrc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !             Fluxes from bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_version(bac))

    case ( BACT1 ) ! Vichi et al. 2007
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! There is no Carbon excretion in Vichi et al. 2007
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
      ! This rate is assumed to occur with a timescale p_ruen=1 day
      ! and controlled with a Michaelis-Menten function
      ! Fix-ratio: actual quota comes from detritus uptake vs. realized growth
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ren  =  (qncPBA(:,bac) - p_qncPBA(bac))*bacc*p_ruen(bac)
      if ( ppbacn == 0 ) & 
         ren = min( ZERO, (ruR1n+ruR6n) - p_qncPBA(bac)*run )
      call quota_flux(iiPel, ppbacn, ppbacn, ppN4n,       ren*insw( ren), tfluxN)
      call quota_flux(iiPel, ppbacn, ppN4n, ppbacn, -eN4n*ren*insw(-ren), tfluxN)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
      ! This rate is assumed to occur with a timescale of p_ruep=1 day
      ! and controlled with a Michaelis-Menten function
      ! Fix-ratio: actual quota comes from detritus uptake vs. realized growth
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rep  =  (qpcPBA(:,bac) - p_qpcPBA(bac))*bacc*p_ruep(bac)
      if ( ppbacp == 0 ) &
         rep = min(ZERO, (ruR1p+ruR6p) - p_qpcPBA(bac)*run )
      call quota_flux(iiPel, ppbacp, ppbacp, ppN1p,       rep*insw( rep), tfluxP)
      call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -eN1p*rep*insw(-rep), tfluxP)

    case ( BACT2 ) ! Vichi et al. 2004
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium remineralization  (Eq. 28 Vichi et al. 2004, note that there 
      ! is a bug in the paper as there should be no division by B1c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      huln = (ruR6n + ruR1n) - p_qncPBA(bac)*run
      ren  = huln*insw(huln)
      call quota_flux(iiPel, ppbacn, ppbacn, ppN4n, ren, tfluxN)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Nitrogen uptake  (Eq. 29 Vichi et al. 2004, there is a bug
      ! here as well, the min should be a max as the numbers are negative!)
      ! from Ammonium and Nitrate when N is not balanced (huln<0)
      ! (nitrate uptake with ammonium inhibition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rumn3 = p_qun(bac)*N3n(:)*bacc*(ONE-eN4n)
      rumn4 = p_qun(bac)*N4n(:)*bacc
      rumn  = rumn3 + rumn4
      ren   = max(-rumn,huln)*insw(-huln)
      call quota_flux(iiPel, ppbacn, ppN4n, ppbacn, -ren*rumn4/rumn, tfluxN)
      call quota_flux(iiPel, ppbacn, ppN3n, ppbacn, -ren*rumn3/rumn, tfluxN)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Phosphate remineralization  (Eq. 28 Vichi et al. 2004, note that there 
      ! is an error in the paper as there should be no division by B1c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      hulp = (ruR6p + ruR1p) - p_qpcPBA(bac)*run
      rep  = hulp*insw(hulp)
      call quota_flux(iiPel, ppbacp, ppbacp, ppN1p, rep, tfluxP)
  
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Phosphorus uptake  (Eq. 29 Vichi et al. 2004, there is a bug
      ! here as well, the min should be a max as the numbers are negative!)
      ! from dissolved phosphate whene P is not balanced (hulp<0)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rump = p_qup(bac)*N1p(:)*bacc
      rep  = max(-rump,hulp)*insw(-hulp)
      call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -rep, tfluxP)
  
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Excess carbon (also considering dissolved nutrient uptake ren and rep) 
      ! is released as R3c, no other excretion (reR2c=0)
      ! (eq. 30 Vichi et al. 2004, unfortunately there is another error in 
      ! the paper, the flux of dissolved nutrient is not written)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      r     = min(run, (ruR6n+ruR1n-ren)/p_qlnc(bac))
      reR3c = run - min(r, (ruR6p+ruR1p-rep)/p_qlpc(bac))
      reR3c = max(ZERO, reR3c)
      call quota_flux( iiPel, ppbacc, ppbacc,ppR3c, reR3c ,tfluxC)

    case ( BACT3 ) ! Polimene et al. (2006)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Carbon excretion as Semi-Labile (R2) and Semi-Refractory (R3) DOC
      ! The R2 rate is assumed to occur with a timescale of 1 day
      ! (eq 8 Polimene et al., 2006)
      ! The renewal of capsular material is a constant rate, equivalent
      ! to about 1/4 of the respiration rate, ~5% of uptake
      ! (Stoderegger and Herndl, 1998)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      reR2c = max((ONE-(qpcPBA(:,bac)/p_qpcPBA(bac))), &
              (ONE-(qncPBA(:,bac)/p_qncPBA(bac))))*p_rec(bac)
      reR2c = max(ZERO,reR2c)*bacc
      reR3c = rug*(ONE-p_pu_ra(bac))*(p_pu_ra(bac)*p_pu_ea_R3(bac))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
      ! This rate is assumed to occur with a timescale p_ruen=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ren = (qncPBA(:,bac) - p_qncPBA(bac))*bacc*p_ruen(bac)
      call quota_flux(iiPel, ppbacn, ppbacn, ppN4n,       ren*insw( ren), tfluxN)
      call quota_flux(iiPel, ppbacn, ppN4n, ppbacn, -eN4n*ren*insw(-ren), tfluxN)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
      ! This rate is assumed to occur with a timescale of p_ruep=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rep  =  (qpcPBA(:,bac) - p_qpcPBA(bac))*bacc*p_ruep(bac)
      call quota_flux(iiPel, ppbacp, ppbacp, ppN1p,       rep*insw( rep), tfluxP)
      call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -eN1p*rep*insw(-rep), tfluxP)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Excretion fluxes (only losses to R2 and R3)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call quota_flux( iiPel, ppbacc, ppbacc, ppR2c, reR2c, tfluxC)
      call quota_flux( iiPel, ppbacc, ppbacc, ppR3c, reR3c, tfluxC)
  
  end select
#ifdef BFM_OGS
  ! Degradation of R3c ( time scale 100 years)
  rrc = 1/(200.0D0*365.D0) * BAC_ACT_FACT * R3c
  call flux_vector(iiPel, ppR3c,    ppO3c, rrc)
  call flux_vector(iiPel, ppO2o,    0    , rrc/ MW_C)
#endif

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

  end subroutine PelBacDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
