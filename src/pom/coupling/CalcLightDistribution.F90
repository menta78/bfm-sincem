#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model verion 2.2
!
! SUBROUTINE
!   CalcLightDistribution
!
! FILE
!   CalcLightDistribution
!
! DESCRIPTION
!   This function #
!	%scalar3%
!
!
!  
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij
!   structure of the code based on ideas of M. Vichi.
!
! AUTHORS
!   !
!
! CHANGE_LOG
!   ! REVISION INFO
! Local  array r (NO_BOXES) added. To track changes search for "zav"
!                                              Marco.Zavatarelli@unibo.it
!                                              May 2008
!
!
!
! COPYING
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
  SUBROUTINE CalcLightDistribution()
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global (scalar) variables are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberY, NO_BOXES_Y, BoxNumberX, &
  ! NO_BOXES_X, BoxNumber, BoxNumberXY, DailyIrrAtSurface
  ! The following 1-d global Pelagic boxvars are modified : EIR
  ! The following 1-d global Pelagic boxvars are used: SUNQ, xEPS, Depth
  ! The following global parameter vars are used: p_PAR
  ! The following global constants are used: RLEN
  ! The following constants are used: HOURS_PER_DAY, E2W

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#ifdef NOPOINTERS
   use mem
#else
  use mem, ONLY: BoxNumberZ, NO_BOXES_Z, BoxNumberY, NO_BOXES_Y, BoxNumberX, &
    NO_BOXES_X, BoxNumber, BoxNumberXY, EIR, SUNQ, xEPS, Depth ,NO_BOXES_XY, &
    NO_BOXES

!    NO_BOXES_X, BoxNumber, BoxNumberXY, DailyIrrAtSurface, EIR, SUNQ, xEPS, Depth ,NO_BOXES_XY
#endif
  use global_mem, ONLY:RLEN
  use constants,  ONLY: HOURS_PER_DAY, E2W

  use mem_Param,  ONLY: p_PAR  ! DailyIrrAtSurface

!  real(RLEN) :: p_PAR

!  use TimeModule, ONLY:cycle,time
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Specifications
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: box_no
!------zav-----
real(RLEN), dimension(NO_BOXES) :: r
!--------------

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: iy
  

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            !	This is a surface box, calculate P-synthetically
            !   available radiance (Aksnes & Lie 1990)
            !   DailyIrrAtSurface is the average irradition per 24 hour
            !   However for EIR the average light in light period is used.
            !   Therefore division by the light_period!
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

             EIR(1) = EIR(1)* &
                p_PAR/ E2W



  DO BoxNumberZ=2,NO_BOXES_Z

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            !	The calculations below are for the lower boxes only.
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            box_no = BoxNumberZ
!------zav-----
             EIR(box_no) = &
               EIR(box_no-1)*  &
               exp( - 1.0D+00* xEPS( box_no-1)* &
               Depth( box_no-1))
!--------------

    end DO
      return
  end

