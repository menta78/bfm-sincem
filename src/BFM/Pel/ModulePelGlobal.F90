#include "cppdefs.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelGlobal
!
! DESCRIPTION
!   !
!   Define parameters and control for sinking of 3D pelagic variables
!
! !INTERFACE
  module mem_PelGlobal
!
! !USES:

  use global_mem
  use mem,       only: NO_D3_BOX_STATES, ppPelDetritus, ppPhytoPlankton,  &
      & iiPhytoPlankton, iiLastElement, iiR6, sediR2, sediR6,sediPPY
  use api_bfm,   only: var_names, BOTindices
  use mem_Phyto, only: p_res, p_rPIm
#if defined INCLUDE_PELCO2
  use mem,       only: sediO5, ppO5c
#endif

!  
!
! !AUTHORS
!   Piet Ruardij
!
! !REVISION_HISTORY
!   Created at Tue Apr 20 09:11:59 AM CEST 2004
!   Updated at Sep 2016 (T. Lovato)
!
! COPYING
!   
!   Copyright (C) 2016 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  ! Define type to control sinking of state variables
  type :: SinkControl
     logical                           :: dosink
     integer                           :: group
     real(RLEN), pointer, dimension(:) :: sedi
  end type SinkControl
  
  type(SinkControl), dimension(NO_D3_BOX_STATES) :: SINKD3STATE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! PelGlobal PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! NAME           UNIT      DESCRIPTION
  ! p_rR6m         [m/d]   detritus sinking rate
  ! p_rO5m         [m/d]   calcite sinking rate
  ! KSINK_rPPY      [m]    prescribe sinking rate for phytoplankton below this 
  !                        depth threshold to p_rR6m value. Use 0.0 to disable. 
  ! AggregateSink  logic   use aggregation = true to enhance the sink rate
  !                        and bypass the prescribed sinking
  ! depth_factor    [m]    depth factor for aggregation method
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)    :: p_rR6m = 0.0_RLEN
  real(RLEN)    :: p_rO5m = 0.0_RLEN
  real(RLEN)    :: KSINK_rPPY = 0.0_RLEN
  logical       :: AggregateSink = .FALSE.
  real(RLEN)    :: depth_factor = 2000.0_RLEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitPelGlobal
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPelGlobal()

  implicit none 

  integer :: n, m, ppstate

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#if ! defined INCLUDE_PELCO2
  namelist /PelGlobal_parameters/ p_rR6m, KSINK_rPPY, AggregateSink, depth_factor
#else   
  namelist /PelGlobal_parameters/ p_rR6m, p_rO5m, KSINK_rPPY, AggregateSink, &
                                  depth_factor
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  LEVEL1 "#  Reading PelGlobal parameters.."
  open(NMLUNIT,file='Pelagic_Environment.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=PelGlobal_parameters,err=101)
  close(NMLUNIT)
  LEVEL1 "#  Namelist is:"
  if (bfm_lwp) write(LOGUNIT,nml=PelGlobal_parameters)

  ! Initialize Pelagic variables settling dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  ! Fill SINKD3STATE array to control sinking of pelagic 3D state variables
  SINKD3STATE(:)%dosink = .FALSE.
  SINKD3STATE(:)%group = 0

  ! Particulate Organic Matter (only R6)
  if ( p_rR6m > 0.0_RLEN) then
     do n = 1, iiLastElement
         ppstate = ppPelDetritus(iiR6,n)
         if (ppstate > 0) then 
            SINKD3STATE(ppstate)%dosink = .TRUE.
            SINKD3STATE(ppstate)%sedi => sediR6
         endif
     enddo
  endif

#if defined INCLUDE_PELCO2
  ! Calcite (O5)
  if ( p_rO5m > 0.0_RLEN) then
     SINKD3STATE(ppO5c)%dosink = .TRUE.      
     SINKD3STATE(ppO5c)%sedi => sediO5   
  endif
#endif

  ! Phytoplankton (if sinking parameters are defined)
  do m = 1 , iiPhytoPlankton
      if ( p_res(m)>0_RLEN .or. p_rPIm(m)>0_RLEN ) then
         do n = 1 , iiLastElement
            ppstate = ppPhytoPlankton(m,n)
            if (ppstate > 0) then 
               SINKD3STATE(ppstate)%dosink = .TRUE.
               SINKD3STATE(ppstate)%group  = 1
               SINKD3STATE(ppstate)%sedi  => sediPPY(m,:)
            endif
         enddo
      endif
  enddo

  ! write log summary of pelagic states sinking setting
  LEVEL1 'SINK setting of pelagic 3D STATES variables'
  LEVEL1 '  ID   Variable   Group'
  do n = 1 , NO_D3_BOX_STATES
     if ( bfm_lwp .and. SINKD3STATE(n)%dosink ) then
     if (allocated(var_names)) then
       write(LOGUNIT,'(i8,a8,i9)') n,trim(var_names(n)),SINKD3STATE(n)%group
     endif
     endif
  enddo
  LEVEL1 ''

  return

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPelGlobal.f90","Pelagic_Environment.nml")
101 call error_msg_prn(NML_READ,"InitPelGlobal.f90","PelGlobal_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  end  subroutine InitPelGlobal

  end module mem_PelGlobal
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
