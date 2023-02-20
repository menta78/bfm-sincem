!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: compute the average field
!
! DESCRIPTION
!   This is a sophisticated subroutine to accumulate instantaneous
!   values and compute averages over each output interval.
!   It stores only the variables defined in the {/tt namelist} variable
!   {/tt ave_save}.
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
   subroutine calcmean_bfm(mode)
!
! USES
   use api_bfm
   use mem, only: D3STATE,D3DIAGNOS,D2DIAGNOS

   use mem, only: NO_BOXES,NO_BOXES_XY

#if defined INCLUDE_SEAICE
   use mem, only: D2STATE_ICE,D2DIAGNOS_ICE,D2FLUX_FUNC_ICE
#endif
   use mem, only: D2STATE_BEN,D2DIAGNOS_BEN,D2FLUX_FUNC_BEN

   implicit none

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,intent(IN)                  :: mode

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    integer                     ::i
    integer                     ::j
    integer                     ::k
    integer                     ::rc
    logical,save                ::do_3ave,do_2ave
#ifdef INCLUDE_SEAICE
    logical,save                ::do_2ave_ice
#endif
    logical,save                ::do_2ave_ben
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#ifdef DEBUG
   LEVEL1 'calcmean_bfm, mode=',mode
#endif

   select case (mode)
      case(INIT)   ! initialization
         i=count(var_ave(stPelStart:stPelFluxE))
         if ( i > 0 ) then
            allocate(D3ave(1:NO_BOXES,1:i),stat=rc)
            if (rc /= 0) stop 'init_bio(): Error allocating D3ave'
            D3ave=0.0
            do_3ave = .true.
         else
            do_3ave = .false.
         endif
         i=count(var_ave(stPelDiag2dS:stPelEnd))
         if ( i > 0 ) then
            allocate(D2ave(1:NO_BOXES_XY,1:i),stat=rc)
            if (rc /= 0) stop 'init_bio(): Error allocating D2ave'
            D2ave=0.0
            do_2ave = .true.
         else
            do_2ave = .false.
         end if
#if defined INCLUDE_SEAICE
         i=count(var_ave(stIceStart:stIceEnd))
         if ( i > 0 ) then
            allocate(D2ave_ice(1:NO_BOXES_XY,1:i),stat=rc)
            if (rc /= 0) stop 'init_bio(): Error allocating D2ave_ice'
            D2ave_ice=0.0
            do_2ave_ice = .true.
         else
            do_2ave_ice = .false.
         end if
#endif
         i=count(var_ave(stBenStart:stBenEnd))
         if ( i > 0 ) then
            allocate(D2ave_ben(1:NO_BOXES_XY,1:i),stat=rc)
            if (rc /= 0) stop 'init_bio(): Error allocating D2ave_ben'
            D2ave_ben=0.0
            do_2ave_ben = .true.
         else
            do_2ave_ben = .false.
         end if
         ave_count=0.0

      case(RESET)
         ave_count=0.0

      case(MEAN)    ! prepare for printing
         if (do_3ave) D3ave=D3ave/ave_count
         if (do_2ave) D2ave=D2ave/ave_count
#ifdef INCLUDE_SEAICE
         if (do_2ave_ice) D2ave_ice=D2ave_ice/ave_count
#endif
         if (do_2ave_ben) D2ave_ben=D2ave_ben/ave_count
         ave_count=0.0

      case(ACCUMULATE)   ! Start of new time-step
         ave_count=ave_count+1.0
         !---------------------------------------------
         ! Compute 3D pelagic means
         !---------------------------------------------
         if (stPelEnd /= 0 .and. do_3ave ) then
            k=0

            j=0
            do i=stPelStateS,stPelStateE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D3ave(:,k)=D3STATE(:,j)
                  else
                     D3ave(:,k)=D3ave(:,k)+D3STATE(:,j)
                  end if
               end if
            end do

            j=0
            do i=stPelDiagS,stPelDiagE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D3ave(:,k)=D3DIAGNOS(:,j)
                  else
                     D3ave(:,k)=D3ave(:,k)+D3DIAGNOS(:,j)
                  end if
               end if
            end do

            j=0
            do i=stPelFluxS,stPelFluxE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  call correct_flux_output(1,j,1,NO_BOXES,c1dim)
                  if ( ave_count < 1.5 ) then
                     D3ave(:,k)=c1dim
                  else
                     D3ave(:,k)=D3ave(:,k)+c1dim
                  end if
               end if
            end do
         endif
         if (stPelDiag2dE /= 0 .and. do_2ave) then
            !---------------------------------------------
            ! Compute Pelagic 2D means
            !---------------------------------------------
            k=0
            j=0
            do i=stPelDiag2dS,stPelBotE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave(:,k)=D2DIAGNOS(:,j)
                  else
                     D2ave(:,k)=D2ave(:,k)+D2DIAGNOS(:,j)
                  end if
               end if
            end do
         end if
#if defined INCLUDE_SEAICE
         if (stIceEnd /= 0 .and. do_2ave_ice) then
            !---------------------------------------------
            ! Compute Seaice 2D means
            !---------------------------------------------
            k=0

            j=0
            do i=stIceStateS,stIceStateE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave_ice(:,k)=D2STATE_ICE(:,j)
                  else
                     D2ave_ice(:,k)=D2ave_ice(:,k)+D2STATE_ICE(:,j)
                  end if
               end if
            end do

            j=0
            do i=stIceDiag2dS,stIceDiag2dE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave_ice(:,k)=D2DIAGNOS_ICE(:,j)
                  else
                     D2ave_ice(:,k)=D2ave_ice(:,k)+D2DIAGNOS_ICE(:,j)
                  end if
               end if
            end do

            j=0
            do i=stIceFlux2dS,stIceFlux2dE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave_ice(:,k)=D2FLUX_FUNC_ICE(:,j)
                  else
                     D2ave_ice(:,k)=D2ave_ice(:,k)+D2FLUX_FUNC_ICE(:,j)
                  end if
               end if
            end do
         end if
#endif
         if (stBenEnd /= 0 .and. do_2ave_ben) then
            !---------------------------------------------
            ! Compute Benthic 2D means
            !---------------------------------------------
            k=0

            j=0
            do i=stBenStateS,stBenStateE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave_ben(:,k)=D2STATE_BEN(:,j)
                  else
                     D2ave_ben(:,k)=D2ave_ben(:,k)+D2STATE_BEN(:,j)
                  end if
               end if
            end do

            j=0
            do i=stBenDiag2dS,stBenDiag2dE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave_ben(:,k)=D2DIAGNOS_BEN(:,j)
                  else
                     D2ave_ben(:,k)=D2ave_ben(:,k)+D2DIAGNOS_BEN(:,j)
                  end if
               end if
            end do

            j=0
            do i=stBenFlux2dS,stBenFlux2dE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave_ben(:,k)=D2FLUX_FUNC_BEN(:,j)
                  else
                     D2ave_ben(:,k)=D2ave_ben(:,k)+D2FLUX_FUNC_BEN(:,j)
                  end if
               end if
            end do
         end if
   end select

end subroutine calcmean_bfm

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
