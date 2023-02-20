!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: Close forcing files in the standalone BFM
!
! DESCRIPTION
!   Read the namelist parameters for the analytical forcings
!   Also initialize additional components for other forcing methods
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
   subroutine end_envforcing_bfm
!
! USES
   use envforcing
   use global_mem, only: RLEN, bfm_lwp, LOGUNIT

   IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   LEVEL1 'end_envforcing_bfm'
    select case (forcing_method)
    case (1) ! analytical forcings
    case (2) ! input data
       LEVEL2 'Closing forcing data file:'
       LEVEL2 trim(forcing_file)
       close(unit_forcing)
       if (use_external_data) then
          LEVEL2 'Closing external data file:'
          LEVEL2 trim(data_file)
          close(unit_data)
       end if
    case (3) ! interactive air-sea fluxes
    end select
    if ( use_benthic_data ) then
       LEVEL2 'Closing benthic forcing data file:'
       LEVEL2 trim(benthic_file)
       close(unit_benthic)
    endif
#ifdef INCLUDE_SEAICE
    if ( use_seaice_data ) then
       LEVEL2 'Closing sea-ice forcing data file:'
       LEVEL2 trim(seaice_file)
       close(unit_seaice)
    endif
#endif

   return

   end subroutine end_envforcing_bfm

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
