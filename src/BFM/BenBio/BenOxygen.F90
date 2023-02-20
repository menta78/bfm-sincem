!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: BenOxygen
!
! DESCRIPTION
!   Description of first order oxic processes in the sediment and computation 
!   of oxygen penetration depth (D1m). 
!   This is based on the stationary state solution of the following: 
!      dGO/dt = Do * d2GO/dz^2 - Mo(bt) - Mo(nit) - Mo(rox)
!   where Do is oxygen molecular diffusion, Mo() are benthic respiration terms, 
!   and GO is pore water oxygen concentration [mmolO2/m2].
!   At equilibrium, the following boundary conditions are defined 
!      Go(z=0)eq = O2o(z=-H)  ;  dGo(z=D1m)eq/dz = 0
!   with O2o as pelagic oxygen and z sediment depth.      
!   By imposing the above conditions to stationary state solution at equilibrium 
!   (dGO/dt=0) the pore water oxygen at equilibrium is 
!      GOeq = O2o - p * D1m *z + p/2 * z^2
!   where p=( Mo(bt) + Mo(nit) + Mo(rox) ) / Do
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
#include "DEBUG.h"
#include "INCLUDE.h"
#include "cppdefs.h"
!
! INTERFACE
  subroutine BenOxygenDynamics
!
! USES
  use global_mem, ONLY:RLEN,LOGUNIT
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: D1m, G2o, D2m
  use mem, ONLY: ppD1m, ppG2o, ppD2m, InitializeModel, shiftD1m, &
    jbotO2o, ETW_Ben, irrenh, rrBTo, jG2K3o, jG2K7o, O2o_Ben, NO_BOXES_XY, iiBen, &
    iiPel, flux_vector
#if defined BENTHIC_FULL
  use mem, ONLY: KNO3, M3n, BoxNumberXY_ben
  use bennut_interface, ONLY: CalculateFromSet
#endif
#endif

  use constants,  ONLY: SEC_PER_DAY, ONE_PER_DAY, STANDARD,EQUATION, ZERO_KELVIN
  use mem_Param,  ONLY: p_poro, p_small, p_d_tot, CalcBenthicFlag
  use mem_BenOxygen
  use mem_Param,  ONLY: p_d_tot, p_clD1D2m
  use time,       ONLY: bfmtime

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: Do, MG2o, p, r
  real(RLEN),dimension(NO_BOXES_XY)  :: D1mNew, G2oNew
  real(RLEN),dimension(NO_BOXES_XY)  :: GO2flux
  real(RLEN)                         :: dummy, delta

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Current timestep for transient G2o computation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  delta = bfmtime%timestep / SEC_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen diffusion from pelagic to benthic pore water [m2/d]
  ! Empirical equation from Broecker and Peng (1974, Table 2 note a)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  Do = SEC_PER_DAY* 1.0e-9_RLEN* (10.0_RLEN)**((- 984.26_RLEN/( 273.0_RLEN+ &
    ETW_Ben(:))+ 3.672_RLEN))* p_poro* irrenh(:)* p_exsaf

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total respiration in pore water from /m2 to /m3:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  MG2o = ( rrBTo(:)+ jG2K3o(:)+ jG2K7o(:) ) / ( D1m(:)* p_poro )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Determine new thickness by imposing the following boundary condition
  !    Go(z)eq = 0 for z > D1m
  ! to pore water oxygen equilibrium equation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  D1mNew = sqrt( 2.0_RLEN* Do* O2o_Ben(:) / (p_small+ mG2o) )
  D1mNew = min(D1mnew ,p_d_tot-2.0_RLEN*p_clD1D2m);

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate rate of change of thickness of the aerobic layer:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  shiftD1m(:)  =  ( max(  p_mD1m,  D1mNew)- D1m(:))/ ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Damping the change of D1m in case of large changes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( InitializeModel== 0) then
     shiftD1m(:) = shiftD1m(:)* (D1m(:)/( D1m(:)+ &
      abs(shiftD1m(:))))**(p_xdampingD1m)*( p_chD1m/( p_chD1m+ D1m(:)))
#if defined BENTHIC_FULL
     do BoxNumberXY_ben = 1,NO_BOXES_XY
        r(1) = CalculateFromSet( KNO3(BoxNumberXY_ben), EQUATION, &
               STANDARD, D1m(BoxNumberXY_ben), dummy)/M3n(BoxNumberXY_ben)
        if ( r(1) .lt.ZERO) then
          write(LOGUNIT,*) "BFM Warning: BenOxygen proportion M3n(D1m)/M3n(0..D2m)=",r(1)
        endif
     end do
     shiftD1m(:)= shiftD1m(:) * max(ZERO,min(ONE,r*2.0_RLEN));
#endif
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Damping the change of D1m in case of too thick D1m
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  r = min( shiftD1m(:),max(ZERO,p_d_tot-p_chD1m-D1m(:)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! recalculate the new D1mNew at the actual time step:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  D1mNew = D1m(:) + r * delta

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! calculate the new pore water oxygen concentration at D1mNew
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  p =  ( 2.0_RLEN* Do* O2o_Ben(:))/( D1mNew* D1mNew)

  G2oNew = D1mNew*( O2o_Ben(:)- 0.66667_RLEN* p * D1mNew* D1mNew/( 2.0_RLEN* &
    Do))* p_poro

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! flux to pelagic: correct flux for rate of change of G2o
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  GO2flux = -( rrBTo(:)+ jG2K3o(:)+ jG2K7o(:))- (G2oNew- G2o(:))/ ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Assign fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( InitializeModel== 0) then
    shiftD1m(:) = r
    call flux_vector( iiBen, ppD1m,ppD1m, r )
    call flux_vector( iiBen, ppG2o,ppG2o,-GO2flux )
    jbotO2o(:) = jbotO2o(:) + GO2flux
#if defined BENTHIC_BIO
    ! Compute shifting of the denitrification layer here in case of running only
    ! the benthic submodel and NOT the benthic nutrient model.
    r  =   p_d_tot- D2m(:)
    call flux_vector( iiBen, ppD2m,ppD2m, shiftD1m(:)* r/( r+ 0.01_RLEN) )
#endif
  else
      G2o(:)=G2oNew
  endif

  end subroutine BenOxygenDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
