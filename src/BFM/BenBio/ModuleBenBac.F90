!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BenBac
!
! DESCRIPTION
!   The parameter value file for benthic bacteria (H1,H2)   
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
!
! INTERFACE
  module mem_BenBac
!
! USES
  use global_mem
  use mem,  ONLY: iiBenBacteria, iiH1, iiH2

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenBac PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_cdm(iiBenBacteria)  ! Half mortality layer        (m)
  real(RLEN)  :: p_qncBBA(iiBenBacteria)  ! Optimal Internal N quota
  real(RLEN)  :: p_qpcBBA(iiBenBacteria)  ! Optimal Internal P quota
  real(RLEN)  :: p_qlnc(iiBenBacteria)  ! Optimal Internal N quota
  real(RLEN)  :: p_qlpc(iiBenBacteria)  ! Optimal Internal P quota
  real(RLEN)  :: p_q10(iiBenBacteria)  ! Q10
  real(RLEN)  :: p_suhQ6(iiBenBacteria)  ! Specific (high) uptake rate (1/d)
  real(RLEN)  :: p_sulQ6(iiBenBacteria)  ! Specific (slow) uptake rate (1/d)
  real(RLEN)  :: p_sum(iiBenBacteria)  ! Potential uptake rate       (1/d)
  real(RLEN)  :: p_pue(iiBenBacteria)  ! Fraction of Q6 degradated as Q1
  real(RLEN)  :: p_suQ1(iiBenBacteria)  ! Specific uptakte rate of Q1 (1/d)
  real(RLEN)  :: p_pur(iiBenBacteria)  ! Fraction of uptake respired
  real(RLEN)  :: p_srr(iiBenBacteria)  ! Specific respiration        (1/d)
  real(RLEN)  :: p_sumKIn(iiBenBacteria)  ! max. uptake of KIn (mmN/m2)
  real(RLEN)  :: p_sumKIp(iiBenBacteria)  ! max. uptake of KIp (mmP/m2)
  real(RLEN)  :: p_sd(iiBenBacteria)  ! Specific mortality          (1/d)
  integer  :: p_iK4(iiBenBacteria)  ! BenthicAmmonium(p_IK4) =1,2 --> K4n K14n
  integer  :: p_iK1(iiBenBacteria)  ! BenthicPhosphate(p_IK1) =1,2 --> K1p K11p
  integer  :: p_iQ1(iiBenBacteria)  ! BenDetritus(p_IQ1) =1,2 --> Q1 Q11

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitBenBac

  contains

  subroutine InitBenBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenBacteria_parameters/ p_cdm, p_qpcBBA, p_qlpc, p_qncBBA, p_qlnc, p_q10, &
    p_suhQ6, p_sulQ6, p_sum, p_pue, p_suQ1, p_pur, p_srr, p_sumKIn, p_sumKIp, &
    p_sd, p_iK4, p_iK1, p_iQ1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading BenBac parameters.."
    open(NMLUNIT,file='Benthic_Ecology.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=BenBacteria_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=BenBacteria_parameters)

  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenBac.f90","Benthic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitBenBac.f90","BenBacteria_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  end subroutine InitBenBac

  end module mem_BenBac

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
