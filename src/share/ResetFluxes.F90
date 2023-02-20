!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE ResetFluxes.F90
!
! DESCRIPTION
!   Reset the arrays for the next integration
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
#include "INCLUDE.h"
!
! INTERFACE
  subroutine ResetFluxes
!
! USES
    use global_mem, only:ZERO,LOGUNIT
#ifdef NOPOINTERS
    use mem
#else

    use mem, ONLY: NO_D3_BOX_STATES,D3SOURCE, &
         PELBOTTOM, PELSURFACE, D3FLUX_FUNC

#if defined INCLUDE_SEAICE
    use mem, ONLY: NO_D2_BOX_STATES_ICE, D2SOURCE_ICE, &
         D2FLUX_FUNC_ICE
#endif
    use mem, ONLY: NO_D2_BOX_STATES_BEN, D2SOURCE_BEN, &

         D2FLUX_FUNC_BEN

#endif

    implicit none

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    integer :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    ! Reset the source and sink arrays
    D3SOURCE(:,:) = ZERO
#if defined INCLUDE_SEAICE
    D2SOURCE_ICE(:,:) = ZERO
#endif
    D2SOURCE_BEN(:,:) = ZERO

    if (allocated(D3FLUX_FUNC)) D3FLUX_FUNC(:,:) = ZERO
#if defined INCLUDE_SEAICE
    if (allocated(D2FLUX_FUNC_ICE)) D2FLUX_FUNC_ICE(:,:) = ZERO
#endif
    if (allocated(D2FLUX_FUNC_BEN)) D2FLUX_FUNC_BEN(:,:) = ZERO

#ifdef BFM_POM
    jbotO2o(:) = ZERO
    jbotO3c(:) = ZERO
    jbotN1p(:) = ZERO
    jbotN3n(:) = ZERO
    jbotN4n(:) = ZERO
#endif

#ifdef BFM_POM
    jsurO2o(:) = ZERO
    jsurO3c(:) = ZERO
#endif

    ! reset surface and bottom fluxes
    do i=1,NO_D3_BOX_STATES
       PELSURFACE(:,i) = ZERO 
       PELBOTTOM(:,i) = ZERO
    end do

  end subroutine ResetFluxes

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
