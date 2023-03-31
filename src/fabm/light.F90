#include "fabm_driver.h"

module ogs_bfm_light
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

   use fabm_types
   use ogs_bfm_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ogs_bfm_light
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_EIR, id_parEIR, id_xEPS
      ! Identifiers for dependencies
      type (type_dependency_id)            :: id_dz, id_xEPSp, id_ESS
      type (type_horizontal_dependency_id) :: id_I_0
!      type (type_state_variable_id)        :: id_P1chl, id_P2chl, id_P3chl, id_P4chl
      type (type_dependency_id)            :: id_chla      
      type (type_state_variable_id)        :: id_X1c, id_X2c, id_X3c

      ! Parameters
      real(rk) :: EPSESSX,EPS0X,pEIR_eowX
      real(rk) :: pEPSCHL,pEPSCDOM
   contains
!     Model procedures
      procedure :: initialize
      procedure :: get_light
   end type type_ogs_bfm_light

contains

   subroutine initialize(self,configunit)


!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_light),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%EPS0X,    'EPS0r',   '1/m',   'background shortwave attenuation', default=4.E-2_rk)
      call self%get_parameter(self%EPSESSX,  'EPSESS',  'm^2/mg','specific shortwave attenuation of silt', default=4.E-5_rk)
      call self%get_parameter(self%pEIR_eowX,'pEIR_eow','-',     'photosynthetically active fraction of shortwave radiation', default=.4_rk)
      call self%get_parameter(self%pEPSCHL,  'pEPSCHL', 'm^2/mg Chl', 'attenuation by Chlorophyll', default=.0088_rk)
      call self%get_parameter(self%pEPSCDOM, 'pEPSCDOM','m^2/mgC', 'attenuation by CDOM', default=.0646_rk)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_EIR,'EIR','uE m-2 d-1','shortwave radiation', &
              standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_parEIR,'parEIR','uE m-2 d-1','photosynthetically active radiation', &
              standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_xEPS,'xEPS','1/m','attenuation coefficient of shortwave flux', &
              source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_xEPSp,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
!      call self%register_state_dependency(self%id_P1chl,'P1chl','mg chl/m^3', 'Diatoms chlorophyll')
!      call self%register_state_dependency(self%id_P2chl,'P2chl','mg chl/m^3', 'Flagellates chlorophyll')
!      call self%register_state_dependency(self%id_P3chl,'P3chl','mg chl/m^3', 'PicoPhytoplankton chlorophyll')
!      call self%register_state_dependency(self%id_P4chl,'P4chl','mg chl/m^3', 'DinoFlagellates chlorophyll')
!     call self%register_dependency(self%id_ESS, type_bulk_standard_variable(name='mass_concentration_of_silt'))
      call self%register_dependency(self%id_chla, total_chlorophyll)      
      call self%register_state_dependency(self%id_X1c,'X1c','mg c/m^3', 'labile CDOM')
      call self%register_state_dependency(self%id_X2c,'X2c','mg c/m^3', 'semi-labile CDOM')
      call self%register_state_dependency(self%id_X3c,'X3c','mg c/m^3', 'semi-refractory CDOM')      
   end subroutine

   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_ogs_bfm_light),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: buffer,dz,xEPS,xtnc,EIR,ESS
!      real(rk) :: P1chl, P2chl, P3chl, P4chl
      real(rk) :: Pchla      
      real(rk) :: X1c, X2c, X3c


      _GET_HORIZONTAL_(self%id_I_0,buffer)

      if (buffer.lt.0._rk) buffer=0._rk

      _VERTICAL_LOOP_BEGIN_
!         _GET_(self%id_P1chl,P1chl)
!         _GET_(self%id_P2chl,P2chl)
!         _GET_(self%id_P3chl,P3chl)
!         _GET_(self%id_P4chl,P4chl)
         _GET_(self%id_chla,Pchla)         
         _GET_(self%id_X1c,X1c)
         _GET_(self%id_X2c,X2c)
         _GET_(self%id_X3c,X3c)

         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_xEPSp,xEPS) ! Extinction coefficient of shortwave radiation, due to particulate organic material (m-1)
!        _GET_(self%id_ESS,ESS)   ! Suspended silt
         xEPS = self%EPS0X + self%pEPSCHL * (Pchla) + self%pEPSCDOM * (X1c + X2c + X3c)
!        xEPS = xEPS + self%EPS0X + self%EPSESSX*ESS
         xtnc = xEPS*dz
         EIR = buffer/xtnc*(1.0_rk-exp(-xtnc))*WtoQuanta*SEC_PER_DAY  ! [uE m-2 d-1]  Note: this computes the vertical average, not the value at the layer centre.
         buffer = buffer*exp(-xtnc)
         _SET_DIAGNOSTIC_(self%id_EIR,EIR)                     ! Local shortwave radiation
         _SET_DIAGNOSTIC_(self%id_parEIR,EIR*self%pEIR_eowX)   ! Local photosynthetically active radiation
         _SET_DIAGNOSTIC_(self%id_xEPS,xEPS)                   ! Vertical attenuation of shortwave radiation
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module
