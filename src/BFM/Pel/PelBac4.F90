#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelBacDynamics4
!
! DESCRIPTION
!   This process describes the dynamics of bacterioplankton
!    
!
! !INTERFACE
  subroutine PelBacDynamics4(bac)
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY:  R6c, R6n, R6p, R1c, R1n, R1p, R2c, O2o, N6r, BGE, &
    N4n, N1p, N3n, R3c, iiR1, iiR6, D3STATE, BRUM, BRUT, TEMP1, TEMP2, TEMP3, TEMP4
  use mem, ONLY: iiPelBacteria, ppPelBacteria, iiC, iiN, iiP, ppR6c, &
    ppR6n, ppR6p, ppR1c, ppR1n, ppR1p, &
    ppR2c, ppO2o, ppN6r, ppN4n, ppN1p, ppN3n, ppR3c, flPTN6r, Depth, ETW, &
    qncPBA, qpcPBA, eO2mO2, qpcOMT, qncOMT, NO_BOXES, iiBen, iiPel, flux_vector,quota_flux
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c
#endif
#endif
  use constants,  ONLY: MW_C, ONE_PER_DAY
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro, p_small
  use mem_PelBac
  use mem_globalfun,   ONLY: eTq, MM_power, insw, MM, fixratio
  use bfm_error_msg, ONLY: bfm_error
!  
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
! !INPUT:
  integer,intent(IN)  :: bac

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer       :: i
  integer       :: ppbacc, ppbacn, ppbacp
  integer, save :: first =0
  integer       :: AllocStatus
  integer,dimension(NO_BOXES)  :: limit
  real(RLEN)    :: p_dnc = 0.86_RLEN  ! N/C ratio of denitrification (Paulmier et al., 2009)
  real(RLEN),allocatable,save,dimension(:) :: et,eO2,r,flN6rPBA,rrc,rro,  &
                                          rd,ruR1c,ruR1n,ruR1p,ruR2c,ruR3c,  &
                                          ruR6c,ruR6p,ruR6n,rump, rtot,rnorm,  &
                                          rumn,rumn3,rumn4,ren,rep,reR2c, &
                                          reR3c,rut,rum,run,sun,rug, &
                                          suR2,cuR6,cuR1,iN1p,iNIn,iN, &
                                          eN1p,eN4n,huln, hulp, bacc, sut, &
                                          tfluxC, tfluxN, tfluxP, pe_N4n, pe_N1p, pe_R6c
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
        &       ruR2c(NO_BOXES), ruR3c(NO_BOXES), rro(NO_BOXES),      &
        &       rump(NO_BOXES), rumn(NO_BOXES), rtot(NO_BOXES),       &
        &       rumn3(NO_BOXES), rumn4(NO_BOXES), ren(NO_BOXES),      &
        &       rep(NO_BOXES), reR2c(NO_BOXES), reR3c(NO_BOXES),      &
        &       rut(NO_BOXES), rum(NO_BOXES), run(NO_BOXES),          &
        &       sun(NO_BOXES), rug(NO_BOXES), suR2(NO_BOXES),         &
        &       cuR6(NO_BOXES), cuR1(NO_BOXES), rnorm(NO_BOXES),      &
        &       iN1p(NO_BOXES), iNIn(NO_BOXES), iN(NO_BOXES),         &
        &       eN1p(NO_BOXES), eN4n(NO_BOXES), huln(NO_BOXES),       &
        &       hulp(NO_BOXES), bacc(NO_BOXES), sut(NO_BOXES),        &
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
  bacc = D3STATE(ppbacc,:)

  ! Quota collectors
  tfluxC = ZERO
  tfluxN = ZERO
  tfluxP = ZERO
  BRUM = ZERO 
  BGE  = ZERO
  
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
  ! Mortality:
  !  1. first order mortality: p_sd 
  !  2. density dependent mortality (e.g. virus infection): p_sd2
  !
  !  It is assumed that mortality is distributed in the same way over
  !  DOC (R1) and detritus (R6) using cytoplasm composition (p_pe_R1x in Param)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd  =  ( p_sd(bac)*et + p_sd2(bac)*bacc ) * bacc

  !call quota_flux(iiPel, ppbacc, ppbacc, ppR6c, rd*(ONE-p_pe_R1c)              , tfluxC)
  !call quota_flux(iiPel, ppbacn, ppbacn, ppR6n, rd*qncPBA(bac,:)*(ONE-p_pe_R1n), tfluxN)
  !call quota_flux(iiPel, ppbacp, ppbacp, ppR6p, rd*qpcPBA(bac,:)*(ONE-p_pe_R1p), tfluxP)

  !call quota_flux(iiPel, ppbacc, ppbacc, ppR1c, rd*p_pe_R1c              , tfluxC)
  !call quota_flux(iiPel, ppbacn, ppbacn, ppR1n, rd*qncPBA(bac,:)*p_pe_R1n, tfluxN) 
  !call quota_flux(iiPel, ppbacp, ppbacp, ppR1p, rd*qpcPBA(bac,:)*p_pe_R1p, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular, always one for fixed stoichiometry)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iNIn = min(ONE, max(ZERO, qncPBA(bac,:)/p_qncPBA(bac)))  !Nitrogen
  iN1p = min(ONE, max(ZERO, qpcPBA(bac,:)/p_qpcPBA(bac)))  !Phosphorus
  iN   = min(iN1p, iNIn)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential uptake by bacteria (eq. 50 Vichi et al. 2007)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rum = iN*et*p_sum(bac)*bacc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! No correction of organic material quality (Polimene et al. 2006)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  cuR1 = ONE
  cuR6 = ONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Realized Substrate Uptake rate for each detritus type and quality (cuRx)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ruR1c = p_suhR1(bac) * cuR1 * R1c(:)
  ruR6c = p_suR6(bac) * cuR6 * R6c(:)
  rut   = p_small + ruR1c + ruR6c
  !
  ! Weighted mean of substrate uptake
  rnorm   = ONE / (p_suhR1(bac)*cuR1 + p_suR6(bac)*cuR6)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria (eq. 50 Vichi et al. 2007)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !BRUT = rum
  !TEMP1 = MM_power( ruR1c*rnorm + ruR6c*rnorm, p_chuc(bac), p_chucMM(bac) )
  !TEMP1 = MM_power( ruR1c + ruR6c, p_chuc(bac)*ruR1c/rut, p_chucMM(bac) ) ! 10s
  TEMP1 = MM_power( rut , p_chuc(bac)*max(p_qsum(bac), ruR1c/rut), p_chucMM(bac) ) ! 10h
  TEMP2 = rum
  TEMP3 = ruR1c/rut
  TEMP4 = rut

  !rug = rum * MM_power( ruR1c*rnorm + ruR6c*rnorm, p_chuc(bac), p_chucMM(bac) ) 
  !rug = rum * MM_power( ruR1c + ruR6c, p_chuc(bac)*ruR1c/rut, p_chucMM(bac) ) ! 10s
  rug = rum * MM_power( rut , p_chuc(bac)*max(p_qsum(bac), ruR1c/rut), p_chucMM(bac) ) 
  ruR1c = rug*ruR1c/rut
  ruR6c = rug*ruR6c/rut

  BRUM = rug
  BRUT = ruR1c + ruR6c

  ! Mortality toward substrates
  call quota_flux(iiPel, ppbacc, ppbacc, ppR6c, rd*ruR6c/rut              , tfluxC)
  call quota_flux(iiPel, ppbacn, ppbacn, ppR6n, rd*qncPBA(bac,:)*ruR6c/rut, tfluxN)
  call quota_flux(iiPel, ppbacp, ppbacp, ppR6p, rd*qpcPBA(bac,:)*ruR6c/rut, tfluxP)

  call quota_flux(iiPel, ppbacc, ppbacc, ppR1c, rd*ruR1c/rut              , tfluxC)
  call quota_flux(iiPel, ppbacn, ppbacn, ppR1n, rd*qncPBA(bac,:)*ruR1c/rut, tfluxN)
  call quota_flux(iiPel, ppbacp, ppbacp, ppR1p, rd*qpcPBA(bac,:)*ruR1c/rut, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Carbon uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call quota_flux(iiPel, ppbacc, ppR1c, ppbacc, ruR1c, tfluxC)
  call quota_flux(iiPel, ppbacc, ppR6c, ppbacc, ruR6c, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR1n = qncOMT(iiR1,:)*ruR1c
  ruR6n = qncOMT(iiR6,:)*ruR6c
  call quota_flux(iiPel, ppbacn, ppR1n, ppbacn, ruR1n, tfluxN)
  call quota_flux(iiPel, ppbacn, ppR6n, ppbacn, ruR6n, tfluxN)

  ruR1p = qpcOMT(iiR1,:)*ruR1c
  ruR6p = qpcOMT(iiR6,:)*ruR6c
  call quota_flux(iiPel, ppbacp, ppR1p, ppbacp, ruR1p, tfluxP)
  call quota_flux(iiPel, ppbacp, ppR6p, ppbacp, ruR6p, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Aerobic respiration
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity + Basal respiration
  rrc = (p_pu_ra(bac) * rug) + p_srs(bac)* bacc* et
  ! Aerobic consumption of Oxygen
  call flux_vector( iiPel, ppO2o, ppO2o, -eO2*rrc/MW_C )
  ! Anaerobic consumption of Nitrate (implicit denitrification)
  rro = p_pu_ra_o(bac) * rug
  call flux_vector( iiPel, ppN3n, ppN4n, (ONE-eO2)* (rrc + rro) / MW_C* p_dnc) 
  ! Total Respiration
  rtot = rrc + (ONE-eO2)*rro
  call quota_flux( iiPel, ppbacc, ppbacc, ppO3c, rtot, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Net Production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  run = rug - rtot

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Growth efficiency
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  WHERE (run.gt.p_small) BGE = run / (run+rtot)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !             Fluxes from bacteria
  !
  ! Vichi et al. 2007 (There is no Carbon excretion) 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved Nitrogen dynamics
  ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
  ! This rate is assumed to occur with a timescale p_ruen=1 day
  ! and controlled with a Michaelis-Menten function
  ! Fix-ratio: actual quota comes from detritus uptake vs. realized growth
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ren  =  (qncPBA(bac,:) - p_qncPBA(bac))*bacc*p_ruen(bac)
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
  rep  =  (qpcPBA(bac,:) - p_qpcPBA(bac))*bacc*p_ruep(bac)
  if ( ppbacp == 0 ) &
     rep = min(ZERO, (ruR1p+ruR6p) - p_qpcPBA(bac)*run )
  call quota_flux(iiPel, ppbacp, ppbacp, ppN1p,       rep*insw( rep), tfluxP)
  call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -eN1p*rep*insw(-rep), tfluxP)

  ! Fix quota control
  if ( ppbacn == 0 .or. ppbacp == 0 ) then

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Eliminate the excess of the non-limiting constituent under fixed quota
     ! Determine release due to either C, P, or N limitation
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     call fixratio(tfluxC,tfluxN,tfluxP,qncPBA(bac,:),qpcPBA(bac,:),pe_R6c, pe_N4n, pe_N1p)

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Apply the correction terms depending on the limiting constituent
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     call flux_vector(iiPel, ppbacc, ppR6c, pe_R6c*(ONE-p_pe_R1c))
     call flux_vector(iiPel, ppbacc, ppR1c, pe_R6c*(p_pe_R1c)    )
     call flux_vector(iiPel, ppbacp, ppN1p, pe_N1p)
     call flux_vector(iiPel, ppbacn, ppN4n, pe_N4n)

  endif

  end subroutine PelBacDynamics4

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
