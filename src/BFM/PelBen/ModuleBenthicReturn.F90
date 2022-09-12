!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: mem_BenthicReturn
!
! DESCRIPTION
!   Setup benthic return interface
!
! COPYING
!
!   Copyright (C) 2022 BFM System Team (bfm_st@cmcc.it)
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
  module mem_BenthicReturn
!
! USES
  use global_mem
  use mem_Param, ONLY: CalcBenthicFlag
  use mem,       ONLY: NO_BOXES_XY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenthicReturn PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_rmnQ1c, &  ! Spec. remin. rate of Dissolved OM (d-1)
                 p_rmnQ1p, p_rmnQ1n 
  real(RLEN)  :: p_rmnQ6c, &  ! Spec. remin. rate of Particulate OM (d-1)
                 p_rmnQ6p,p_rmnQ6n, p_rmnQ6s
  real(RLEN)  :: p_pQIN3      ! Partitioning coeff. between NO3 and NH4
  real(RLEN)  :: p_depscale   ! Depth level to scale remineralization rates (m)
  real(RLEN)  :: p_q10        ! Q10 value for remineralization rate
  real(RLEN)  :: p_qBT        ! Base Temperature for Q10 remin
  real(RLEN)  :: p_chdo       ! Half-saturation Oxygen concentration for MM^2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Scaling factor for remineralization rates
  real(RLEN),ALLOCATABLE, DIMENSION(:) :: RETFAC
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitBenthicReturn

  contains

  subroutine InitBenthicReturn()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenthicReturn_parameters/ p_rmnQ1c, p_rmnQ1n, p_rmnQ1p,            &
                                      p_rmnQ6c, p_rmnQ6n, p_rmnQ6p, p_rmnQ6s,  &
                                      p_pQIN3, p_depscale, p_q10, p_qBT, p_chdo
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if ( CalcBenthicFlag ) then
     write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
     write(LOGUNIT,*) "#  Reading BenthicReturn parameters.."
     open(NMLUNIT,file='Benthic_Environment.nml',status='old',action='read',err=100)
     read(NMLUNIT,nml=BenthicReturn_parameters,err=101)
     close(NMLUNIT)
     write(LOGUNIT,*) "#  Namelist is:"
     write(LOGUNIT,nml=BenthicReturn_parameters)
  endif
  !
  ! Scaling factor of remineralzation rates
  allocate(RETFAC(NO_BOXES_XY))
  RETFAC = ONE
  ! 
  if ( p_depscale > ZERO) &
     write(LOGUNIT,'(a,f12.6)') "Set depth scaling factor for Benthic Return rates using ", p_depscale

  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenthicReturn1.f90","Benthic_Environment.nml")
101 call error_msg_prn(NML_READ,"InitBenthicReturn1.f90","BenthicReturn_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  end  subroutine InitBenthicReturn

  end module mem_BenthicReturn

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
