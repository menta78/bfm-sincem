!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BenOrganism
!
! DESCRIPTION
!   The parameter value file for benthic groups (Y1,Y2,Y4,Y5)
!   The values for suspension feeders (Y3) are in a separate file
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
  module mem_BenOrganisms
!
! USES
  use global_mem
  use mem,  ONLY: iiBenOrganisms, iiBenBacteria ,  &
                  iiY1, iiY2, iiY3, iiY4, iiY5

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenOrganism PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_q10(iiBenOrganisms)  !
  real(RLEN)  :: p_su(iiBenOrganisms)  !
  real(RLEN)  :: p_chu(iiBenOrganisms)  !
  real(RLEN)  :: p_clu(iiBenOrganisms)  !
  real(RLEN)  :: p_pue(iiBenOrganisms)  !
  real(RLEN)  :: p_pur(iiBenOrganisms)  !
  real(RLEN)  :: p_pudil(iiBenOrganisms)  !
  real(RLEN)  :: p_sr(iiBenOrganisms)  !
  real(RLEN)  :: p_sd(iiBenOrganisms)  !
  real(RLEN)  :: p_sdm(iiBenOrganisms)
  real(RLEN)  :: p_paBOS(iiBenOrganisms,iiBenOrganisms)  ! Y1 > YI
  real(RLEN)  :: p_paBBA(iiBenOrganisms,iiBenBacteria)  ! H1 > YI
  real(RLEN)  :: p_puQ6(iiBenOrganisms)  ! H2 > YI
  real(RLEN)  :: p_pueQ6(iiBenOrganisms)
  real(RLEN)  :: p_cm(iiBenOrganisms)
  real(RLEN)  :: p_clm(iiBenOrganisms)
  real(RLEN)  :: p_qncBOS(iiBenOrganisms)
  real(RLEN)  :: p_qpcBOS(iiBenOrganisms)
#ifdef BFM_POM
  real(RLEN)  :: p_clO2o(iiBenOrganisms)
  real(RLEN)  :: p_sdo(iiBenOrganisms)
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitBenOrganisms

  contains

  subroutine InitBenOrganisms()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenOrganisms_parameters/ p_q10, p_su, p_chu, p_clu, p_pue, p_pur, &
    p_pudil, p_sr, p_puQ6, p_pueQ6, p_cm, p_clm, p_sd, p_sdm, p_qncBOS, p_qpcBOS, &
#ifdef BFM_POM
    p_clO2o, p_sdo, &
#endif
    p_paBOS, p_paBBA
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading BenOrganisms parameters.."
    open(NMLUNIT,file='Benthic_Ecology.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=BenOrganisms_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=BenOrganisms_parameters)

  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenOrganisms.f90","Benthic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitBenOrganisms.f90","BenOrganisms_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  end subroutine InitBenOrganisms

  end module mem_BenOrganisms

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
