!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! ROUTINE: MesoZoo
!
! DESCRIPTION
!   This file contains the parameter values for the mesozooplankton submodel.
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
#include "cppdefs.h"
!
! INTERFACE
  module mem_MesoZoo
!
! USES
  use global_mem
  use mem,  ONLY: iiMesoZooPlankton, &
                  iiPhytoPlankton, iiMicroZooplankton

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  !-------------------------------------------------------------------------!
  !NAMELIST MesoZoo_parameters
  !-------------------------------------------------------------------------!
  !  MESO-ZOOPLANKTON
  ! NAME         [UNIT]/KIND            DESCRIPTION
  ! p_q10        [-]             Q10 value for physiological rates
  ! p_srs        [1/d]           Respiration rate at 10 degrees C
  ! p_sum        [1/d]           Maximal productivity at 10 degrees C
  ! p_sd         [1/d]           Background natural mortality
  ! p_vum        [m3/mgC/d]      Specific search volume
  ! p_puI        [-]             Assimilation efficiency
  ! p_peI        [-]             Fraction of Faeces production
  ! p_pecaco3    [-]             Portion of egested calcified shells during grazing
  ! p_sdo        [m3/mgC/d]      Specific density-dependent mortality
  ! p_sds        [-]             Exponent of density-dependent mortality
  ! p_qpcMEZ     [mmolP/mgC]     Maximum quotum P:C
  ! p_qncMEZ     [mmolN/mgC]     Maximum quotum N:C
  ! p_paPPY(z,p) [-]             Availability of PhytoPlankton group p
  !                              to Zooplankton group z
  ! p_paMIZ(z,m) [-]             Availability of MicroZooplankton group m
  !                              to Zooplankton group z
  ! p_paMEZ(z,m) [-]             Availability of MesoZooplankton group m
  !                              to Zooplankton group z
  !-------------------------------------------------------------------------!
  real(RLEN)  :: p_q10(iiMesoZooPlankton)
  real(RLEN)  :: p_srs(iiMesoZooPlankton)
  real(RLEN)  :: p_sum(iiMesoZooPlankton)
  real(RLEN)  :: p_sd(iiMesoZooPlankton)
  real(RLEN)  :: p_vum(iiMesoZooPlankton)
  real(RLEN)  :: p_puI(iiMesoZooPlankton)
  real(RLEN)  :: p_peI(iiMesoZooPlankton)
  real(RLEN)  :: p_pecaco3(iiMesoZooPlankton)
  real(RLEN)  :: p_sdo(iiMesoZooPlankton)
  real(RLEN)  :: p_sds(iiMesoZooPlankton)
  real(RLEN)  :: p_qpcMEZ(iiMesoZooPlankton)
  real(RLEN)  :: p_qncMEZ(iiMesoZooPlankton)
  real(RLEN)  :: p_clO2o(iiMesoZooPlankton)
  real(RLEN)  :: p_paPPY(iiMesoZooPlankton,iiPhytoPlankton)
  real(RLEN)  :: p_paMIZ(iiMesoZooPlankton,iiMicroZooPlankton)
  real(RLEN)  :: p_paMEZ(iiMesoZooPlankton,iiMesoZooPlankton)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitMesoZoo

  contains

  subroutine InitMesoZoo()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /MesoZoo_parameters/ p_q10, p_srs, p_paPPY, p_paMIZ, p_paMEZ, p_sd, &
    p_sum, p_vum, p_puI, p_peI, p_sdo, p_sds, p_qpcMEZ, p_qncMEZ, p_clO2o,     &
    p_pecaco3
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  LEVEL1 "#  Reading MesoZoo parameters.."
  open(NMLUNIT,file='Pelagic_Ecology.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=MesoZoo_parameters,err=101)
  close(NMLUNIT)
  LEVEL1 "#  Namelist is:"
  if (bfm_lwp) write(LOGUNIT,nml=MesoZoo_parameters)

  return

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitMesoZoo.f90","Pelagic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitMesoZoo.f90","MesoZoo_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end subroutine InitMesoZoo

  end module mem_MesoZoo

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
