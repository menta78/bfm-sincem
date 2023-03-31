#include "fabm_driver.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   This process describes the calcite dissolution 
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
! !INTERFACE

 module bfm_CalciteDissolution

   use fabm_types

   use ogs_bfm_shared
   use ogs_bfm_pelagic_base

   implicit none

   private

   type,extends(type_ogs_bfm_pelagic_base), public :: type_ogs_bfm_CalciteDissolution
!     Variable identifiers
      type (type_state_variable_id)      :: id_O5c
      type (type_diagnostic_variable_id) :: id_fCaCO3_2_O3c,id_fCaCO3_2_O3h
      type (type_state_variable_id)      :: id_O3c,id_O3h
      type (type_dependency_id)          :: id_OCalc

      real(rk) :: p_kdca, p_nomega, p_wsinkPIC
   contains
      procedure :: initialize
      procedure :: do
   end type type_ogs_bfm_CalciteDissolution

contains

   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_CalciteDissolution), intent(inout), target :: self
      integer,                              intent(in)            :: configunit
      real(rk) :: c0
!EOP
!-----------------------------------------------------------------------
!BOC
      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are
      ! given in d-1.
      self%dt = 86400._rk

      call self%get_parameter(self%p_wsinkPIC, 'p_wsinkPIC', 'm d-1', 'sinking velocity of PIC')
      call self%register_state_variable(self%id_O5c,'c','mg C/m^3','PIC calcite',vertical_movement=-self%p_wsinkPIC/self%dt,initial_value=1._rk)

      call self%get_parameter(self%p_kdca   ,'p_kdca'   ,'1/d','maximum specific dissolution rate')
      call self%get_parameter(self%p_nomega,'p_nomega','-'  ,'power of the dissolution law ')
 
      call self%get_parameter(c0,'c0','mg C/m^3','background concentration',default=0.0_rk)

      call self%register_diagnostic_variable(self%id_fCaCO3_2_O3c,'fCaCO2_2_O3c','mg C/m^3/d','calcite dissolution flux')
      call self%register_diagnostic_variable(self%id_fCaCO3_2_O3h,'fCaCO2_2_O3h','mmol/m3/d','production of O3h for calcite dissolution')
      call self%register_dependency(self%id_OCalc,'OCalc','-','calcite saturation') !!! viene da pelagicCO2.F90
      call self%register_state_dependency(self%id_O3c,'O3c','mg C/m^3','total dissolved inorganic carbon')
      call self%register_state_dependency(self%id_O3h,'O3h','mmol/m^3','alkalinity')
   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ogs_bfm_CalciteDissolution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
     real(rk) :: O5c
     real(rk) :: OCalc
     real(rk) :: rdiss
     real(rk) :: excess
      _LOOP_BEGIN_
         _GET_(self%id_O5c,O5c)
         _GET_(self%id_OCalc,OCalc)

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Dissolution of Particulate Inorganic Carbon (calcite/aragonite) in seawater
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute undersaturation
  excess = max(ZERO,ONE - OCalc)
  ! Dissolution rate of C in CaCO3 (mg C/m3/d) from Morse and Berner (1972)
  rdiss = self%p_kdca * excess**self%p_nomega * O5c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Inorganic carbon and alkalinity flux due to PIC changes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  call flux_vector( iiPel, ppO5c, ppO3c, rdiss(:) )
         _SET_ODE_(self%id_O5c, -rdiss)
         _SET_ODE_(self%id_O3c, rdiss)
!  call flux_vector( iiPel, ppO3h, ppO3h, -C2ALK*rdiss(:) )
         _SET_ODE_(self%id_O3h,C2ALK*rdiss)  ! Dissolution of CaCO3 increases alkalinity by 2 units
!!! ATTENZIONE IN BFM origin e' -C2ALK ma e' un errore => comunicare a Tomas e Paolo
         _SET_DIAGNOSTIC_(self%id_fCaCO3_2_O3c,rdiss)
         _SET_DIAGNOSTIC_(self%id_fCaCO3_2_O3h,C2ALK*rdiss)

      _LOOP_END_

   end subroutine do

end module

