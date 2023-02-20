!-----------------------------------------------------------------------
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: Light and other environmental forcings used in the BFM
!
! DESCRIPTION
! This routine sets the environmental forcings according to user
! choice. The currently implemented methods are:
!    1) analytical sinusoidal forcings
!    2) external values from files
!    3) interactive slab ocean from meteorological data
! Any combination of the method above must be provided by the user
! in an additional subroutine.
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
   subroutine envforcing_bfm(step)
!
! USES
   use global_mem, only: RLEN
   use envforcing

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer        :: step
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    select case (forcing_method)
    case (1) ! analytical forcings
      call analytical_forcing
    case (2) ! input data
      call analytical_forcing
      call external_forcing
#if defined BENTHIC_BIO || defined BENTHIC_FULL
      if ( use_benthic_data ) call external_benthic
#endif
    case (3) ! interactive air-sea fluxes
!      call do_air_sea(timesec,startime)
    end select
    ! Assign external data
    call external_data
    ! Assign external event data
    call event_data
#ifdef INCLUDE_SEAICE
    if ( use_seaice_data ) call external_seaice
#endif
    if (init_forcing_vars) init_forcing_vars=.false.

  end subroutine envforcing_bfm

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
