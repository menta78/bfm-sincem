!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: ClearMem
!
! DESCRIPTION
!   Deallocate all the arrays
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
subroutine ClearMem
!
! USES
   use mem
   use api_bfm

   implicit none

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer :: i,j,origin,destination
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    ! free willy, free the fluxes
    if (allocated(D3FLUX_MATRIX)) then
    origin=0
    do i=stPelStateS,stPelStateE
       origin=origin+1
       destination=0
       do j=stPelStateS,stPelStateE
          destination=destination+1
          if( allocated(D3FLUX_MATRIX(origin,destination)%p) ) deallocate(D3FLUX_MATRIX(origin,destination)%p)
       end do
    end do
    deallocate(D3FLUX_MATRIX)
    if( allocated(D3FLUX_FUNC) ) deallocate(D3FLUX_FUNC)
    end if
#if defined INCLUDE_SEAICE
    if (allocated(D2FLUX_MATRIX_ICE)) then
    origin=0
    do i=stIceStateS,stIceStateE
       origin=origin+1
       destination=0
       do j=stIceStateS,stIceStateE
          destination=destination+1
          if( allocated(D2FLUX_MATRIX_ICE(origin,destination)%p) ) deallocate(D2FLUX_MATRIX_ICE(origin,destination)%p)
       end do
    end do
    deallocate(D2FLUX_MATRIX_ICE)
    if( allocated(D2FLUX_FUNC_ICE) ) deallocate(D2FLUX_FUNC_ICE)
    end if
#endif
    if (allocated(D2FLUX_MATRIX_BEN)) then
    origin=0
    do i=stBenStateS,stBenStateE
       origin=origin+1
       destination=0
       do j=stBenStateS,stBenStateE
          destination=destination+1
          if( allocated(D2FLUX_MATRIX_BEN(origin,destination)%p) ) deallocate(D2FLUX_MATRIX_BEN(origin,destination)%p)
       end do
    end do
    deallocate(D2FLUX_MATRIX_BEN)
    if( allocated(D2FLUX_FUNC_BEN) ) deallocate(D2FLUX_FUNC_BEN)
    end if

    ! from api_bfm 
    deallocate(var_ids)
    deallocate(var_ave)
    deallocate(var_names)
    deallocate(var_units)
    deallocate(var_long)
    deallocate(c1dim)
    if (allocated(D3ave)) deallocate(D3ave)
    if (allocated(D2ave)) deallocate(D2ave)

#ifndef NOT_STANDALONE
     deallocate(D3STATE)
     deallocate(D3SOURCE)
     deallocate(D3STATETYPE)
     deallocate(D3DIAGNOS)
     deallocate(D2DIAGNOS)

#if defined INCLUDE_SEAICE
     if ( allocated(D2ave_ice) ) deallocate(D2ave_ice)
     deallocate(D2DIAGNOS_ICE)
     deallocate(D2STATE_ICE)
     deallocate(D2SOURCE_ICE)
     deallocate(D2STATETYPE_ICE)
#endif

     if ( allocated(D2ave_ben) ) deallocate(D2ave_ben)
     deallocate(D2DIAGNOS_BEN)
     deallocate(D2STATE_BEN)
     deallocate(D2SOURCE_BEN)
     deallocate(D2STATETYPE_BEN)

#endif

end subroutine ClearMem

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
