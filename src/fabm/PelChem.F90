#include "fabm_driver.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelChem
!
! DESCRIPTION
!       This process describes the additional dynamics of dissolved
!       compounds in the watercolumn. Parameterized processes are:
!       - nitrification
!       - denitrification
!       - reoxidation of reduction equivalents
!       - dissolution of biogenic silica
!       This function also calls the carbonate system dynamics
!       (INCLUDE_PELCO2) and iron dynamics (INCLUDE_PELFE)
!       if activated
!
! !INTERFACE
 module bfm_PelChem

   use fabm_types
   use ogs_bfm_shared
   use ogs_bfm_pelagic_base
  !  
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  private

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_PelChem
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving
      ! from type_ogs_bfm_pelagic_base!
      
      ! Identifiers for state variables of other models
      type (type_state_variable_id) :: id_O3c,id_O2o,id_O3h,id_O4n          !  dissolved inorganic carbon, oxygen, total alkalinity, N2
      type (type_state_variable_id) :: id_N3n,id_N4n,id_N5s,id_N6r          !  nutrients: phosphate, nitrate, ammonium, silicate, iron
      type (type_state_variable_id) :: id_R1c,id_R1p,id_R1n,id_R2c,id_R3c   !  dissolved organic carbon (R1: labile, R2: semi-labile, R3: semi-refractory)
      type (type_state_variable_id) :: id_R6c,id_R6p,id_R6n,id_R6s          !  particulate organic carbon
      type (type_state_variable_id) :: id_O5c                               !  Free calcite (liths) - used by calcifiers only
      type (type_state_variable_id) :: id_X1c, id_X2c, id_X3c               !  CDOM
      ! Environmental dependencies
      type (type_dependency_id)    :: id_parEIR,id_ETW   ! PAR and temperature
      type (type_dependency_id)    :: id_PAR_tot

      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_eo      ! Oxygen limitation factor MM
      type (type_diagnostic_variable_id) :: id_er      ! Reduction equiv. limitation factor MM
      type (type_diagnostic_variable_id) :: id_flN4N3n ! flux of nitrification
      type (type_diagnostic_variable_id) :: id_flN4N3n_o2  ! consume of O2 for nitrification
      type (type_diagnostic_variable_id) :: id_flN3O4n ! denitrification
      type (type_diagnostic_variable_id) :: id_flN3O4n_N6r ! impact of denitrification on reduction equivalent
      type (type_diagnostic_variable_id) :: id_fN6O2r  ! reoxydation of reduction equivalents
      type (type_diagnostic_variable_id) :: id_fR6N5s  ! dissolution of biogenic silicate
      type (type_diagnostic_variable_id) :: id_degX1c  ! photodegradation cdom X1c
      type (type_diagnostic_variable_id) :: id_degX2c  ! photodegradation cdom X2c
      type (type_diagnostic_variable_id) :: id_degX3c  ! photodegradation cdom X3c
      type (type_diagnostic_variable_id) :: id_remX3c  ! remineralization of cdom X3c
      type (type_diagnostic_variable_id) :: id_remR3c  ! remineralization of dom R3c      
      type (type_diagnostic_variable_id) :: id_varO3h_for_nitr ! variation of O3h for nitrification (-1 mole of NH4 -> + 2 mole of alk)
      type (type_diagnostic_variable_id) :: id_varO3h_for_denitr ! variation of O3h for denitrification (-1 mole of NO3 (consumed) -> + 1 mole of alk)
      ! Parameters (described in subroutine initialize, below)
      real(rk) :: p_clO2o, p_clN6r, p_sN4N3, p_q10N4N3, p_qon_nitri, p_qro
      real(rk) :: p_sN3O4n, p_rPAo, p_qon_dentri, p_rOS, p_sR6N5, p_q10R6N5      
      real(rk) :: p_bX1c, p_bX2c, p_bX3c, p_IX1, p_IX2, p_IX3, p_rX3c, p_q10X    
      integer :: p_Esource
  
    contains     

   ! Model procedures
      procedure :: initialize
      procedure :: do 
     end type type_ogs_bfm_PelChem


contains

  subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_PelChem),intent(inout),target :: self
      integer,                              intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      real(rk) :: pippo1
      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are
      ! given in d-1.
!      self%dt = 86400._rk
!EOP
!-----------------------------------------------------------------------
!BOC
! Initialize pelagic base model (this also sets the time unit to per day,instead of the default per second)
      call self%initialize_bfm_base

     ! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used
      ! by FABM (or its host)
      ! to present parameters to the user for configuration (e.g., through a
      ! GUI)
      
      call self%get_parameter(self%p_clO2o,       'p_clO2o',     '[mmolO2/m3]' ,    'Half-saturation O2 concentration for nitrification and reoxidation')
      call self%get_parameter(self%p_clN6r,       'p_clN6r',     '[mmolHS/m3]' ,    'Half-saturation concentration of reduction equivalents for denitrification[mmolO2/m3]')
      call self%get_parameter(self%p_sN4N3,       'p_sN4N3',     '[1/d]'       ,    'Specific nitrification rate at 10 degC')
      call self%get_parameter(self%p_q10N4N3,     'p_q10N4N3',   '[-]'         ,    'Q10 factor for nitrification/denit')
      call self%get_parameter(self%p_qon_nitri,   'p_qon_nitri', '[mmolO2/mmolN]',  'Stoichiometric coefficient for nitrification')
      call self%get_parameter(self%p_qro,         'p_qro',       '[mmolHS-/mmolO2]','Stoichiometric coefficient for anaerobic reactions')
      call self%get_parameter(self%p_sN3O4n,      'p_sN3O4n',    '[1/d]'      ,     'Specific denitrification rate')
      call self%get_parameter(self%p_rPAo,        'p_rPAo',      '[mmolO2/m3/d]' ,  'Reference anoxic mineralization rate')
      call self%get_parameter(self%p_qon_dentri,  'p_qon_dentri','[mmolO2/mmolN]',  'Stoichiometric coefficient for denitrification')
      call self%get_parameter(self%p_rOS,         'p_rOS',       '[1/d]',           'Specific reoxidation rate of reduction equivalents')
      call self%get_parameter(self%p_sR6N5,       'p_sR6N5',     '[1/d]',           'Specific remineralization rate of biogenic silica')
      call self%get_parameter(self%p_q10R6N5,     'p_q10R6N5',   '[-]',             'Q10 factor for biogenic silica')
      call self%get_parameter(self%p_bX1c,        'p_bX1c',      '[1/d]',           'photodegradation (bleaching) rate X1c')
      call self%get_parameter(self%p_bX2c,        'p_bX2c',      '[1/d]',           'photodegradation (bleaching) rate X2c')
      call self%get_parameter(self%p_bX3c,        'p_bX3c',      '[1/d]',           'photodegradation (bleaching) rate X3c')
      call self%get_parameter(self%p_IX1,         'p_IX1',       '[uE m-2 s-1]',    'light threshold for X1c bleaching')
      call self%get_parameter(self%p_IX2,         'p_IX2',       '[uE m-2 s-1]',    'light threshold for X2c bleaching')
      call self%get_parameter(self%p_IX3,         'p_IX3',       '[uE m-2 s-1]',    'light threshold for X3c bleaching')
      call self%get_parameter(self%p_rX3c,        'p_rX3c',      '[1/d]',           'degradation rate X3c')
      call self%get_parameter(self%p_q10X,        'p_q10X',      '[-]',             'q10 for CDOM degradation')
      call self%get_parameter(self%p_Esource,     'p_Esource',   '5-6',             'source of light for CDOM bleaching')

      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_O3h,'O3h','mmol /m^3','alkalinity')
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O2/m^3','dissolved oxygen')
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3','nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
      call self%register_state_dependency(self%id_N5s,'N5s','mmol Si/m^3','silicate')
      call self%register_state_dependency(self%id_N6r,'N6r','mmol HS/m^3','reduction equivalent')
      call self%register_state_dependency(self%id_R3c,'R3c','mgC /m^3','semi-refractory DOC')
      call self%register_state_dependency(self%id_R6s,'R6s','mmol Si/m^3','biogenic silicate')
      call self%register_state_dependency(self%id_X1c,'X1c','mgC/m^3','labile cdom')
      call self%register_state_dependency(self%id_X2c,'X2c','mgC/m^3','semilabile cdom')
      call self%register_state_dependency(self%id_X3c,'X3c','mgC/m^3','semi refractory cdom')
      call self%register_state_dependency(self%id_O4n,'O4n','mmolN/m^3','N2')

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_parEIR,standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      ! Dependency from multispectral model
      call self%register_dependency(self%id_PAR_tot,type_bulk_standard_variable(name='PAR_tot'))      

      call self%register_diagnostic_variable(self%id_eo,        'eo'  ,       '-',           'oxygen regulating factor with MichelisMenten',output=output_none)
      call self%register_diagnostic_variable(self%id_er,        'er'  ,       '-',           'reduction equiv. regulating factor with MichelisMenten',output=output_none)
      call self%register_diagnostic_variable(self%id_flN4N3n,   'flN4N3n' ,   'mmolN/m3/d',  'nitrification',output=output_none)
      call self%register_diagnostic_variable(self%id_flN4N3n_o2,'flN4N3n_o2', 'mmolO/m3/d',  'consume of O2 for nitrification',output=output_none)
      call self%register_diagnostic_variable(self%id_flN3O4n,   'flN3O4n',    'mmolN/m3/d',  'denitrification',output=output_none)
      call self%register_diagnostic_variable(self%id_flN3O4n_N6r,'flN3O4n_N6r','mmolHS/m3/d',  'impact of denitrification on reduction equivalent',output=output_none)
      call self%register_diagnostic_variable(self%id_fN6O2r,    'fN6O2r',     'mmolHS/m3/d',  'reoxydation of reduction equivalents',output=output_none)
      call self%register_diagnostic_variable(self%id_fR6N5s,    'fR6N5s',     'mmolSi/m3/d',  'dissolution of biogenic silicate',output=output_none)
      call self%register_diagnostic_variable(self%id_degX1c,    'degX1c',     'mgC/m3/d',  'photodegradation of cdom X1c')
      call self%register_diagnostic_variable(self%id_degX2c,    'degX2c',     'mgC/m3/d',  'photodegradation of cdom X2c')
      call self%register_diagnostic_variable(self%id_degX3c,    'degX3c',     'mgC/m3/d',  'photodegradation of cdom X3c')
      call self%register_diagnostic_variable(self%id_remX3c,    'remX3c',     'mgC/m3/d',  'remineralization of cdom X3c')
      call self%register_diagnostic_variable(self%id_remR3c,    'remR3c',     'mgC/m3/d',  'remineralization of dom R3c')      
!     call self%register_diagnostic_dependency(self%id_flPTN6r,'flPTN6r','mmolHS/m3/d','total rate of formation of reduction equivalent') ! from PelBac
      call self%register_diagnostic_variable(self%id_varO3h_for_nitr,'varO3h_for_nitr','mmol/m3/d','O3h increase for nitrification',output=output_none)
      call self%register_diagnostic_variable(self%id_varO3h_for_denitr,'varO3h_for_denitr','mmol/m3/d','O3h increase for denitrification',output=output_none)

  end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ogs_bfm_PelChem),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW, parEIR
      real(rk) :: N5s,N3n,N4n,O2o,N6r, R3c, R6s, O4n, O3h
      real(rk) :: X1c, X2c, X3c
      real(rk) :: eo, er
      real(rk) :: flN4N3n, flN4N3n_o2, flN3O4n, flN3O4n_N6r, fN6O2r, rPAo,fR6N5s
      real(rk) :: degX1c, degX2c, degX3c
      real(rk) :: remX3c, remR3c


     ! Enter spatial loops (if any)
      _LOOP_BEGIN_

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         ! Retrieve ambient nutrient concentrations
         _GET_(self%id_N3n,N3n)
         _GET_(self%id_N4n,N4n)
         _GET_(self%id_N5s,N5s)
         _GET_(self%id_N6r,N6r)
         _GET_(self%id_R3c,R3c)
         _GET_(self%id_R6s,R6s)

         _GET_(self%id_O2o,O2o)
         _GET_(self%id_O3h,O3h)
         _GET_(self%id_O4n,O4n)
         _GET_(self%id_X1c,X1c)
         _GET_(self%id_X2c,X2c)
         _GET_(self%id_X3c,X3c)

         ! Retrieve environmental dependencies (water temperature,
         ! photosynthetically active radation)
         _GET_(self%id_ETW,ETW)
         ! From where to get the light
         ! Both parEIR and PAR_tot are in uE m-2 d-1, 6=parEIR from light, 5=PAR_tot from light_spectral
         select case (self%p_Esource)
         case (5)
            _GET_(self%id_PAR_tot,   parEIR)   ! uE m-2 d-1
         case (6)
            _GET_(self%id_parEIR,    parEIR)   ! uE m-2 d-1
         end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Regulating factors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  eo  =   MM(  max(p_small,O2o),  self%p_clO2o)
  er  =   MM(  N6r,  self%p_clN6r)

     _SET_DIAGNOSTIC_(self%id_eo,eo)
     _SET_DIAGNOSTIC_(self%id_er,er)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Nitrification in the water
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  flN4N3n =  max(ZERO,self%p_sN4N3* N4n* eTq(  ETW,  self%p_q10N4N3) * eo)

  _SET_DIAGNOSTIC_(self%id_flN4N3n,flN4N3n)

  flN4N3n_o2 =  -( flN4N3n* self%p_qon_nitri)

  _SET_DIAGNOSTIC_(self%id_flN4N3n_o2,flN4N3n_o2) ! consume of o2 for nitrafication

  
! call flux_vector( iiPel, ppN4n,ppN3n, flN4N3n(:) )
 _SET_ODE_(self%id_N3n,flN4N3n)
 _SET_ODE_(self%id_N4n,-flN4N3n)  
! call flux_vector( iiPel, ppO2o,ppO2o,-( flN4N3n(:)* p_qon_nitri) )
 _SET_ODE_(self%id_O2o,flN4N3n_o2)
 _SET_ODE_(self%id_O3h, -2._rk*flN4N3n)  ! Alkalinity contributions: +1 for NH4, -1 for nitrate

 _SET_DIAGNOSTIC_(self%id_varO3h_for_nitr, -2*flN4N3n)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Denitrification in the water
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  _GET_(self%id_flPTN6r,flPTN6r)
  rPAo  =   flPTN6r/ self%p_qro    !!!! ATTENTION flPTN6r will be defined in computed in PelBac
  flN3O4n = max(ZERO,self%p_sN3O4n* eTq( ETW, self%p_q10N4N3)* er* rPAo/ self%p_rPAo* &
               N3n)
  
   _SET_DIAGNOSTIC_(self%id_flN3O4n,flN3O4n) ! denitrification
 
   flN3O4n_N6r = -( self%p_qro* flN3O4n*self%p_qon_dentri*insw(-(O2o-N6r/self%p_qro)))

   _SET_DIAGNOSTIC_(self%id_flN3O4n_N6r,flN3O4n_N6r) ! impact of denitrification on reduction equivalent
               
!  call flux_vector( iiPel, ppN3n,ppO4n, flN3O4n(:) )
 _SET_ODE_(self%id_N3n,-flN3O4n)
 _SET_ODE_(self%id_O4n,flN3O4n)  
!  call flux_vector( iiPel, ppN6r,ppN6r,-( p_qro* flN3O4n(:)* p_qon_dentri* &
!      insw( -( O2o(:)- N6r(:)/ p_qro))) )
 _SET_ODE_(self%id_N6r,flN3O4n_N6r)

 _SET_ODE_(self%id_O3h,flN3O4n)   ! - 1 mole of N3n -> + 1 mole of alkalinity
 _SET_DIAGNOSTIC_(self%id_varO3h_for_denitr,flN3O4n)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Reoxidation of reduction equivalents
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  fN6O2r  =   self%p_rOS* N6r* eo
    _SET_DIAGNOSTIC_(self%id_fN6O2r,fN6O2r) ! reoxydation of reduction equivalents


! call flux_vector( iiPel, ppN6r,ppN6r,-( fN6O2r) )
 _SET_ODE_(self%id_N6r,-fN6O2r)       
! call flux_vector( iiPel, ppO2o,ppO2o,-( fN6O2r/ p_qro) )
 _SET_ODE_(self%id_O2o,-fN6O2r/ self%p_qro)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Dissolution of biogenic silicate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  fR6N5s  =   self%p_sR6N5* eTq(  ETW,  self%p_q10R6N5)* R6s
   _SET_DIAGNOSTIC_(self%id_fR6N5s,fR6N5s) ! dissolution of biogenic silicate

! call flux_vector( iiPel, ppR6s,ppN5s, fR6N5s )
 _SET_ODE_(self%id_N5s,fR6N5s)
 _SET_ODE_(self%id_R6s,-fR6N5s)

!GP #ifdef INCLUDE_PELFE
!GP   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!GP   !  Dissolved Iron Chemistry (dissolution and scavenging)
!GP   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!GP   ! Linear regeneration of bioavailable iron
!GP   fR1N7f(:)  =  p_sR1N7* eTq(  ETW(:),  p_q10R6N7)* R1f(:)
!GP   call flux_vector( iiPel, ppR1f, ppN7f, fR1N7f(:) )
!GP 
!GP  fR6N7f(:)  =  p_sR6N7* eTq(  ETW(:),  p_q10R6N7)* R6f(:)
!GP  call flux_vector( iiPel, ppR6f, ppN7f, fR6N7f(:) )
!GP
!GP  ! Scavenging of free dissolved Iron
!GP  ! Adsorption onto Particulate Organic Matter (Eq.11 in Parekh et al.,2015 - GBC )
!GP  fscavN7f(:) = max ( ZERO, p_scavOrg * N7f * (R6c(:)**0.58_RLEN) )
!GP  
!GP  ! Inorganic component ( Linear relaxation to the Iron Ligand concentration )
!GP  fscavN7f(:) = fscavN7f(:) + max(ZERO,p_scavIng*(N7f-p_N7fLigand))
!GP call flux_vector( iiPel, ppN7f, ppN7f, -fscavN7f(:) )
!GP
!GP#endif

!GP TO BE MOVED IN A NEW FILE 
!GP #ifdef INCLUDE_PELCO2
!GP  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!GP   ! Carbonate chemistry
!GP   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!GP   call PelagicCSYS()
!GP #endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Degradation of cDOM --> photodegradation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!GP  BIOPTIMOD T1
!GP  degCDOM = R1l * ( eTq( ETW(:), 2.95D0 ) * 0.003D0 + 0.167D0 * min(PAR(:)/60.0D0,1.0D0) ) ! Eq 13
!GP  BIOPTIMOD T2
!GP   PAR(:) =EIR(:)

! Check unit of measure of PAR here!   parEIR is in uE m-2 d-1, p_IXn in uE m-2 s-1                             
  degX1c = X1c * ( self%p_bX1c * min(parEIR/(self%p_IX1*SEC_PER_DAY),1.0_rk) ) ! Eq 13
  degX2c = X2c * ( self%p_bX2c * min(parEIR/(self%p_IX2*SEC_PER_DAY),1.0_rk) ) ! Eq 13
  degX3c = X3c * ( self%p_bX3c * min(parEIR/(self%p_IX3*SEC_PER_DAY),1.0_rk) ) ! Eq 13
  
!EA  degR1l = R1l * ( 0.167D0 * min(PAR(:)/60.0D0,1.0D0) )
!EA  degR2l = R2l * ( 0.167D0 * min(PAR(:)/60.0D0,1.0D0) )
!EA  degR3l = R3l * ( eTq( ETW(:), 2.95D0 ) * 0.00003D0 + 0.167D0 * min(PAR(:)/60.0D0,1.0D0) )
  
  _SET_DIAGNOSTIC_(self%id_degX1c,degX1c) ! photodegradation of labile CDOM
  _SET_DIAGNOSTIC_(self%id_degX2c,degX2c) ! photodegradation of semi-labile CDOM
  _SET_DIAGNOSTIC_(self%id_degX3c,degX3c) ! photodegradation of semi-refractory CDOM
  
! call flux_vector( iiPel, ppR1l, ppR3c, degR1l )
 _SET_ODE_(self%id_R3c,degX1c)
 _SET_ODE_(self%id_X1c,-degX1c)
! call flux_vector( iiPel, ppR2l, ppR3c, degR2l )
 _SET_ODE_(self%id_R3c,degX2c)
 _SET_ODE_(self%id_X2c,-degX2c)
! call flux_vector( iiPel, ppR3l, ppR3c, degR3l )
 _SET_ODE_(self%id_R3c,degX3c)
 _SET_ODE_(self%id_X3c,-degX3c)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Remineralization of semi-recalcitrant DOC
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  remX3c = X3c * ( eTq( ETW, self%p_q10X ) * self%p_rX3c )
  remR3c = R3c * ( eTq( ETW, self%p_q10X ) * self%p_rX3c )
    
  _SET_DIAGNOSTIC_(self%id_remX3c,remX3c) ! remineralization of semi-refractory CDOM
  _SET_DIAGNOSTIC_(self%id_remR3c,remR3c) ! remineralization of semi-refractory DOC
  
 _SET_ODE_(self%id_O3c, remX3c)
 _SET_ODE_(self%id_X3c,-remX3c)

 _SET_ODE_(self%id_O3c, remR3c)
 _SET_ODE_(self%id_R3c,-remR3c)

     ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do
end module

! GP !EOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
