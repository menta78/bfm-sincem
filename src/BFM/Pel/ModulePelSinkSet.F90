!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PelGlobal
!
! DESCRIPTION
!   Define parameters and control for sinking of 3D pelagic variables
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
  module mem_PelSinkSet
!
! USES
  use global_mem
  use mem,       only: NO_D3_BOX_STATES, ppPelDetritus, ppPhytoPlankton,  &
      & iiPhytoPlankton, iiLastElement, iiR6, iiR3, sediR2, sediR3, sediR6,sediPPY
  use api_bfm,   only: var_names, BOTindices
  use mem_Phyto, only: p_res, p_rPIm
#if defined INCLUDE_PELCO2
  use mem,       only: sediO5, ppO5c
#endif

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
  integer(RLEN), allocatable, dimension(:) :: sink_var_map
  integer :: sink_vars = 0

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Pelagic sinking PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! NAME           UNIT      DESCRIPTION
  ! p_rR6m         [m/d]   detritus sinking rate
  ! p_rO5m         [m/d]   calcite sinking rate
  ! KSINK_rPPY      [m]    prescribe sinking rate for phytoplankton below this 
  !                        depth threshold to p_rR6m value. Use 0.0 to disable. 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)    :: p_rR6m = 0.0_RLEN
  real(RLEN)    :: p_rR3m = 0.0_RLEN
  real(RLEN)    :: p_rO5m = 0.0_RLEN
  real(RLEN)    :: KSINK_rPPY = 0.0_RLEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Pelagic Settling PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BURIAL VELOCITIES into the sediment
  ! NAME         [UNIT]/KIND            DESCRIPTION
  ! p_burvel_R6     [m/d]              Bottom Burial Velocity for detritus
  ! p_burvel_R2     [m/d]              Bottom Burial Velocity for dissolved
  ! p_burvel_PI     [m/d]              Bottom Burial Velocity for plankton
  ! p_burvel_O5     [m/d]              Bottom Burial Velocity for calcite
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_burvel_R6=0.0_RLEN
  real(RLEN)  :: p_burvel_R2=0.0_RLEN
  real(RLEN)  :: p_burvel_R3=0.0_RLEN
  real(RLEN)  :: p_burvel_PI=0.0_RLEN
  real(RLEN)  :: p_burvel_O5=0.0_RLEN
  logical     :: R6DeepBurial = .FALSE.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitPelSinkSet

  contains

  subroutine InitPelSinkSet()

  implicit none 

  integer :: n, m, ppstate

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#if ! defined INCLUDE_PELCO2
  namelist /PelGlobal_parameters/ p_rR6m, p_rR3m,  KSINK_rPPY
#else   
  namelist /PelGlobal_parameters/ p_rR6m, p_rR3m, p_rO5m, KSINK_rPPY
#endif
  namelist /Settling_parameters/ p_burvel_R6, p_burvel_R2, p_burvel_R3, p_burvel_PI, &
                                 p_burvel_O5, R6DeepBurial

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  LEVEL1 "#  Reading Pelagic Sinking parameters.."
  open(NMLUNIT,file='Pelagic_Environment.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=PelGlobal_parameters,err=101)
  close(NMLUNIT)
  LEVEL1 "#  Namelist is:"
  if (bfm_lwp) write(LOGUNIT,nml=PelGlobal_parameters)

  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  LEVEL1 "#  Reading Pelagic Settling parameters.."
  open(NMLUNIT,file='Benthic_Environment.nml',status='old',action='read',err=102)
  read(NMLUNIT,nml=Settling_parameters,err=103)
  close(NMLUNIT)
  LEVEL1 "#  Namelist is:"
  if (bfm_lwp) write(LOGUNIT,nml=Settling_parameters)

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
            SINKD3STATE(ppstate)%sedi => sediR6(:)
         endif
     enddo
  endif

  ! refractory DOC (R3)
  if ( p_rR3m > 0.0_RLEN) then
     do n = 1, iiLastElement
         ppstate = ppPelDetritus(iiR3,n)
         if (ppstate > 0) then
            SINKD3STATE(ppstate)%dosink = .TRUE.
            SINKD3STATE(ppstate)%sedi => sediR3(:)
         endif
     enddo
  endif

#if defined INCLUDE_PELCO2
  ! Calcite (O5)
  if ( p_rO5m > 0.0_RLEN) then
     SINKD3STATE(ppO5c)%dosink = .TRUE.      
     SINKD3STATE(ppO5c)%sedi => sediO5(:)
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
               SINKD3STATE(ppstate)%sedi  => sediPPY(:,m)
            endif
         enddo
      endif
  enddo

  ! write log summary of pelagic states sinking setting
  LEVEL1 'SINK setting of pelagic 3D STATES variables'
  LEVEL1 '  ID   Variable   Group'
  sink_vars = COUNT(SINKD3STATE(:)%dosink)
  ALLOCATE(sink_var_map(sink_vars))
  m = 0
  do n = 1 , NO_D3_BOX_STATES
     if ( SINKD3STATE(n)%dosink ) then
        m = m + 1
        sink_var_map(m) = n  
        if (allocated(var_names)) then
          if (bfm_lwp) write(LOGUNIT,'(i8,a8,i9)')  sink_var_map(m),trim(var_names(n)),SINKD3STATE(n)%group
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
102 call error_msg_prn(NML_OPEN,"InitPelGlobal.f90","Benthic_Environment.nml")
103 call error_msg_prn(NML_READ,"InitPelGlobal.f90","Settling_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end subroutine InitPelSinkSet

  end module mem_PelSinkSet

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
