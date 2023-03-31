#include "fabm_driver.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelBac
!
! DESCRIPTION
!   This process describes the dynamics of all phytoplankton
!    groups. The differences in behaviour
!    are expressed by differences in parameter-values only.
!    
! !INTERFACE
 module bfm_PelBac

   use fabm_types
   use ogs_bfm_shared
!  use fabm_particle

   use ogs_bfm_pelagic_base

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   implicit none

   private
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_pelagic_bacteria
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving
      ! from type_ogs_bfm_pelagic_base!

      ! Identifiers for state variables of other models
      type (type_state_variable_id) :: id_O3c,id_O2o,id_O3h                 !  dissolved inorganic carbon, oxygen, total alkalinity
      type (type_state_variable_id) :: id_N1p,id_N3n,id_N4n                 !  nutrients: phosphate, nitrate, ammonium, silicate, iron
      type (type_state_variable_id) :: id_R1c,id_R1p,id_R1n,id_R2c,id_R3c   !  dissolved organic carbon (R1: labile, R2: semi-labile)
      type (type_state_variable_id) :: id_R6c,id_R6p,id_R6n,id_R6s          !  particulate organic carbon
      type (type_state_variable_id) :: id_X1c,id_X2c,id_X3c                 !  CDOM
      type (type_state_variable_id) :: id_N6r              
      type (type_state_variable_id) :: id_O5c                               !  Free calcite (liths) - used by calcifiers only
      ! Environmental dependencies
      type (type_dependency_id)            :: id_ETW   ! temperature
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_eT    ! temperature q10 factor
      type (type_diagnostic_variable_id) :: id_ETWd  ! temperature Celsius
      type (type_diagnostic_variable_id) :: id_eO2   ! Oxygen limitation
      type (type_diagnostic_variable_id) :: id_eN4n  ! External nutrient limitation
      type (type_diagnostic_variable_id) :: id_eN1p  ! External nutrient limitation
      type (type_diagnostic_variable_id) :: id_rd    ! Mortality
      type (type_diagnostic_variable_id) :: id_rum   ! Potential uptake by bacteria
      type (type_diagnostic_variable_id) :: id_cuR1  ! correction of organic material quality
      type (type_diagnostic_variable_id) :: id_cuR6  ! correction of organic material quality
      type (type_diagnostic_variable_id) :: id_iNIn  ! internal quota nitrogen limitation 
      type (type_diagnostic_variable_id) :: id_iN1p  ! internal quota phosphorus limitation 
      type (type_diagnostic_variable_id) :: id_iN    ! P and N limitation
      type (type_diagnostic_variable_id) :: id_ruR1c ! Labile DOC realized uptake
      type (type_diagnostic_variable_id) :: id_ruR2c ! Semi-Labile DOC realized uptake
      type (type_diagnostic_variable_id) :: id_ruR3c ! Semi-Refractory DOC realized uptake
      type (type_diagnostic_variable_id) :: id_ruR6c ! POC realized uptake
      type (type_diagnostic_variable_id) :: id_ruX1c ! Labile CDOM realized uptake
      type (type_diagnostic_variable_id) :: id_ruX2c ! Semi-Labile CDOM realized uptake
      type (type_diagnostic_variable_id) :: id_ruX3c ! Semi-Refractory CDOM realized uptake
      type (type_diagnostic_variable_id) :: id_rut   ! realized substrate uptake
      type (type_diagnostic_variable_id) :: id_rug   ! 
      type (type_diagnostic_variable_id) :: id_ruR1n ! Organic Nitrogen uptake
      type (type_diagnostic_variable_id) :: id_ruR6n ! Organic Nitrogen uptake
      type (type_diagnostic_variable_id) :: id_ruR1p ! Organic Nitrogen uptake
      type (type_diagnostic_variable_id) :: id_ruR6p ! Organic Nitrogen uptake
      type (type_diagnostic_variable_id) :: id_rrc   ! total respiration
      type (type_diagnostic_variable_id) :: id_flN6rPBA ! electron acceptor in the respiration process
      type (type_diagnostic_variable_id) :: id_runr    !Net production
      type (type_diagnostic_variable_id) :: id_renh4   !Ammonium remineralization
      type (type_diagnostic_variable_id) :: id_rumn3   ! Nitrate uptake
      type (type_diagnostic_variable_id) :: id_rumn4   ! Ammonium uptake
      type (type_diagnostic_variable_id) :: id_rumn    ! Nitrogen uptake
      type (type_diagnostic_variable_id) :: id_ren     ! actual Inorganic Nitrogen uptake
      type (type_diagnostic_variable_id) :: id_rep     ! Direct uptake of phosphate
      type (type_diagnostic_variable_id) :: id_repo4   ! Phosphate remineralization
      type (type_diagnostic_variable_id) :: id_reR2c   ! Carbon excretion as Semi-Labile
      type (type_diagnostic_variable_id) :: id_reR3c   ! Excess carbon
      type (type_diagnostic_variable_id) :: id_run     ! Net production
      type (type_diagnostic_variable_id) :: id_B_N_O3h   ! variation of alk due to net N uptake/release
      type (type_diagnostic_variable_id) :: id_B_P_O3h   ! variation of alk due to net P uptake/release
      type (type_diagnostic_variable_id) :: id_varO3h_reminNH4 ! variation of alk due to NH4 remineralization

      ! Parameters (described in subroutine initialize, below)
      real(rk) :: p_q10, p_chdo, p_sd, p_sd2, p_suhR1, p_sulR1, p_suR2
      real(rk) :: p_suR3, p_suR6, p_sum, p_pu_ra, p_pu_ra_o, p_srs, p_qncPBA
      real(rk) :: p_qpcPBA, p_qlnc, p_qlpc, p_qun, p_qup, p_chn, p_chp
      real(rk) :: p_ruen, p_ruep, p_rec, p_pu_ea_R3,p_qro
      real(rk) :: p_pe_R1c, p_pe_R1n, p_pe_R1p
      real(rk) :: p_fX1b, p_fX2b, p_fX3b      
      integer :: p_version

   contains

      ! Model procedures
      procedure :: initialize
      procedure :: do

   end type type_ogs_bfm_pelagic_bacteria

   ! Constants

contains

  subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_pelagic_bacteria),intent(inout),target :: self
      integer,                              intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer, parameter ::  BACT1=1,BACT2=2,BACT3=3
!EOP
!-----------------------------------------------------------------------
!BOC
     ! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used
      ! by FABM (or its host)
      ! to present parameters to the user for configuration (e.g., through a
      ! GUI)
      call self%get_parameter(self%p_version, 'p_version', '-',       'Switch for bacteria parameterization')
      call self%get_parameter(self%p_q10,     'p_q10',     '-',       'Characteristic Q10 coefficient')
      call self%get_parameter(self%p_chdo,    'p_chdo',    'mmol/m3', 'Half-saturation constant for O2 limitation')
      call self%get_parameter(self%p_sd,      'p_sd',      '1/d',     'Specific mortality rate')
      call self%get_parameter(self%p_sd2,     'p_sd2',     '1/d',     'Density dependent specific mortality rate')
      call self%get_parameter(self%p_suhR1,   'p_suhR1',   '1/d',     'Specific potential uptake for nutrient-rich DOM')
      call self%get_parameter(self%p_sulR1,   'p_sulR1',   '1/d',     'Specific potential uptake for nutrient-poor DOM')
      call self%get_parameter(self%p_suR2,    'p_suR2',    '1/d',     'Specific potential uptake for semi-labile DOC')
      call self%get_parameter(self%p_suR3,    'p_suR3',    '1/d',     'Specific potential uptake for semi-refractory DOC')
      call self%get_parameter(self%p_suR6,    'p_suR6',    '1/d',     'Specific potential uptake for POM (1/d)')
      call self%get_parameter(self%p_sum,     'p_sum',     '1/d',     'Potential specific growth rate')
      call self%get_parameter(self%p_pu_ra,   'p_pu_ra',     '-',     'Activity respiration fraction')
      call self%get_parameter(self%p_pu_ra_o, 'p_pu_ra_o',   '-',     'Additional respiration fraction at low O2 conc')
      call self%get_parameter(self%p_srs,     'p_srs',     '1/d',     'Specific rest respiration')
      call self%get_parameter(self%p_qncPBA, 'p_qncPBA', 'mmolN/mgC', 'Optimal N/C ratio')
      call self%get_parameter(self%p_qpcPBA, 'p_qpcPBA', 'mmolP/mgC', 'Optimal P/C ratio')
      call self%get_parameter(self%p_qlnc,   'p_qlnc'  , 'mmolN/mgC', 'Minimal N/C ratio')
      call self%get_parameter(self%p_qlpc,   'p_qlpc'  , 'mmolP/mgC', 'Minimal P/C ratio')
      call self%get_parameter(self%p_qun,   'p_qun'  , 'mmolN/mgC/day', ' Membrane affinity for N')
      call self%get_parameter(self%p_qup,   'p_qup'  , 'mmolP/mgC/day', ' Membrane affinity for P')
      call self%get_parameter(self%p_chn,   'p_chn'  , 'mmolN/m3', 'Half saturation ammonium conc. for uptake ')
      call self%get_parameter(self%p_chp,   'p_chp'  , 'mmolP/m3', 'Half saturation phosphate conc. for uptake ')
      call self%get_parameter(self%p_ruen,   'p_ruen'  , '1/d', 'Relaxation timescale for N uptake/remin. ')
      call self%get_parameter(self%p_ruep,   'p_ruep'  , '1/d', 'Relaxation timescale for P uptake/remin. ')
      call self%get_parameter(self%p_rec,      'p_rec'  , '1/d', 'Relaxation timescale for semi-labile excretion')
      call self%get_parameter(self%p_pu_ea_R3, 'p_pu_ea_R3'  , '-', 'Excretion of semi-refractory DOC')
      call self%get_parameter(self%p_qro, 'p_qro'  , 'mmolHS-/mmolO2', 'Stoichiometric coefficient for anaerobic reactions')
      call self%get_parameter(self%p_pe_R1c, 'p_pe_R1c'  , '-', 'Fractional content of C in cytoplasm')
      call self%get_parameter(self%p_pe_R1n, 'p_pe_R1n'  , '-', 'Fractional content of N in cytoplasm')
      call self%get_parameter(self%p_pe_R1p, 'p_pe_R1p'  , '-', 'Fractional content of P in cytoplasm')
!              --------- Flux partition CDOM parameters ------------
      call self%get_parameter(self%p_fX1b,   'p_fX1b',  '-',  'colored fraction in labile dissolved organic carbon', default=0.02_rk)
      call self%get_parameter(self%p_fX2b,   'p_fX2b',  '-',  'colored fraction in semi-labile dissolved organic carbon', default=0.02_rk)
      call self%get_parameter(self%p_fX3b,   'p_fX3b',  '-',  'colored fraction in semi-refractory dissolved organic carbon', default=0.02_rk)
      
! Register state variables (handled by type_bfm_pelagic_base)
      call self%initialize_bfm_base()
      call self%add_constituent('c',1.e-4_rk)
      call self%add_constituent('n',1.26e-6_rk)
      call self%add_constituent('p',4.288e-8_rk)
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O2/m^3','dissolved oxygen')
      call self%register_state_dependency(self%id_O3c,'O3c','mg C/m^3','dissolved organic carbon')
      call self%register_state_dependency(self%id_O3h,'O3h','mmol /m^3','alkalinity')
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3','nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
      call self%register_state_dependency(self%id_N6r,'N6r','mmol Eq/m^3','Reduction Equivalent')
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3','labile DOC')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','labile DOP')
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','labile DON')
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','semi-labile DOC')
      call self%register_state_dependency(self%id_R3c,'R3c','mg C/m^3','semi-refractory DOC')
      call self%register_state_dependency(self%id_R6c,'R6c','mg C/m^3','POC')
      call self%register_state_dependency(self%id_R6p,'R6p','mmol P/m^3','POP')
      call self%register_state_dependency(self%id_R6n,'R6n','mmol N/m^3','PON')
      call self%register_state_dependency(self%id_X1c,'X1c','mg C/m^3','labile CDOM')
      call self%register_state_dependency(self%id_X2c,'X2c','mg C/m^3','semi-labile CDOM')
      call self%register_state_dependency(self%id_X3c,'X3c','mg C/m^3','semi-refractory CDOM')
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      ! Register diagnostic variables (i.e., model outputs)
!      call self%register_diagnostic_variable(self%id_ETWd, 'ETW',  'C','temperature Celsius',output=output_none)
      call self%register_diagnostic_variable(self%id_et,   'et',   '-','temperature factor',output=output_none)
      call self%register_diagnostic_variable(self%id_eO2, 'eO2',   '-','Oxygen limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_eN4n, 'eN4n', '-','External nutrient limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_eN1p, 'eN1p', '-','External nutrient limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_rd,  'rd',  'mgC/m3/d', 'Mortality')
      call self%register_diagnostic_variable(self%id_rum, 'rum', 'mgC/m3/d', 'Potential uptake by bacteria')
      call self%register_diagnostic_variable(self%id_cuR1, 'cuR1', '-','correction of organic material quality',output=output_none)
      call self%register_diagnostic_variable(self%id_cuR6, 'cuR6', '-','correction of organic material quality',output=output_none)
      call self%register_diagnostic_variable(self%id_iNIn, 'iNIn', '-','Nutrient limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_iN1p, 'iN1p', '-','Nutrient limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_iN,   'iN',   '-','Total Nutrient limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_ruR1c,'ruR1c', 'mgC/m3/d','Labile DOC realized uptake')
      call self%register_diagnostic_variable(self%id_ruR2c,'ruR2c', 'mgC/m3/d','Semi-Labile DOC realized uptake')
      call self%register_diagnostic_variable(self%id_ruR3c,'ruR3c', 'mgC/m3/d','Semi-Refractory DOC realized uptake')
      call self%register_diagnostic_variable(self%id_ruR6c,'ruR6c', 'mgC/m3/d','POC realized uptake')
      call self%register_diagnostic_variable(self%id_ruX1c,'ruX1c', 'mgC/m3/d','labile-CDOM realized uptake')
      call self%register_diagnostic_variable(self%id_ruX2c,'ruX2c', 'mgC/m3/d','semilabile-CDOM realized uptake')
      call self%register_diagnostic_variable(self%id_ruX3c,'ruX3c', 'mgC/m3/d','semirefractory-CDOM realized uptake')
      call self%register_diagnostic_variable(self%id_rut,  'rut', 'mgC/m3/d','realized substrate uptake')
      call self%register_diagnostic_variable(self%id_rug,  'rug', 'mgC/m3/d','Actual uptake by bacteria')
      call self%register_diagnostic_variable(self%id_ruR1n,'ruR1n', 'mmolN/m3/d',' Organic Nitrogen uptake')
      call self%register_diagnostic_variable(self%id_ruR6n,'ruR6n', 'mmolN/m3/d',' Organic Nitrogen uptake')
      call self%register_diagnostic_variable(self%id_ruR1p,'ruR1p', 'mmolP/m3/d',' Organic Phosphorus uptake')
      call self%register_diagnostic_variable(self%id_ruR6p,'ruR6p', 'mmolP/m3/d',' Organic Phosphorus uptake',output=output_none)
      call self%register_diagnostic_variable(self%id_rrc,'rrc', 'mgC/m3/d','Aerobic and anaerobic respiration',output=output_none)
      call self%register_diagnostic_variable(self%id_flN6rPBA,'flN6rPBA', '-','electron acceptor in the respiration process',output=output_none)
      select case (self%p_version)
         case ( BACT1 )
          call self%register_diagnostic_variable(self%id_renh4,'ren', 'mmolN/m3/d','Direct uptake of ammonium',output=output_none)
          call self%register_diagnostic_variable(self%id_rep,'rep', 'mmolP/m3/d','Direct uptake of phosphate',output=output_none)
         case ( BACT2 )
          call self%register_diagnostic_variable(self%id_run,'run', 'mgC/m3/d','Net production',output=output_none)
          call self%register_diagnostic_variable(self%id_renh4,'renh4', 'mmolN/m3/d','Ammonium remineralization',output=output_none)
          call self%register_diagnostic_variable(self%id_rumn3,'rumn3', 'mmolN/m3/d','Nitrate uptake',output=output_none)
          call self%register_diagnostic_variable(self%id_rumn4,'rumn4', 'mmolN/m3/d','Ammonium uptake',output=output_none)
          call self%register_diagnostic_variable(self%id_rumn,'rumn', 'mmolN/m3/d','Nitrogen uptake',output=output_none)
          call self%register_diagnostic_variable(self%id_ren,'ren', 'mmolN/m3/d','actual Inorganic Nitrogen uptake',output=output_none)
          call self%register_diagnostic_variable(self%id_rep,'rep', 'mmolP/m3/d','Direct uptake of phosphate',output=output_none)
          call self%register_diagnostic_variable(self%id_repo4,'repo4', 'mmolP/m3/d','Phosphate remineralization',output=output_none)
          call self%register_diagnostic_variable(self%id_reR3c, 'reR3c', 'mmolC/m3/d','Excess carbon')
          call self%register_diagnostic_variable(self%id_B_P_O3h,'varO3h_for_Putil', 'mmol/m3/d','variation of alk due to bacteria P utilization',output=output_none)
          call self%register_diagnostic_variable(self%id_B_N_O3h,'varO3h_for_Nutil', 'mmol/m3/d','variation of alk due to bacteria N utilization',output=output_none)
          call self%register_diagnostic_variable(self%id_varO3h_reminNH4,'varO3h_reminNH4', 'mmol/m3/d','variation of alk due to NH4 remineralization',output=output_none)
         case ( BACT3 ) 
          call self%register_diagnostic_variable(self%id_reR2c,'reR2c', 'mgC/m3/d','Carbon excretion as Semi-Labile')
          call self%register_diagnostic_variable(self%id_reR3c,'reR3c', 'mgC/m3/d','Carbon excretion as Semi-Refractory')
          call self%register_diagnostic_variable(self%id_ren,'ren', 'mmolN/m3/d','Dissolved Nitrogen dynamics',output=output_none)
          call self%register_diagnostic_variable(self%id_rep,'rep', 'mmolP/m3/d','Dissolved Phosphorus dynamics',output=output_none)
  end select
   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ogs_bfm_pelagic_bacteria),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! LOCAL VARIABLES:
      integer, parameter ::  BACT1=1,BACT2=2,BACT3=3
      real(rk) :: bacc, bacp, bacn
      real(rk) :: O2o,N1p,N3n,N4n
      real(rk) :: R1c,R1p,R1n
      real(rk) :: R2c
      real(rk) :: R3c
      real(rk) :: X1c,X2c,X3c
      real(rk) :: R6c,R6p,R6n
      real(rk) :: N6r
      real(rk) :: ETW,et,eO2
      real(rk) :: eN4n,eN1p
      real(rk) :: qpcPBA,qncPBA
      real(rk) :: qpcR1,qncR1,qpcR6,qncR6
      real(rk) :: rd
      real(rk) :: rum
      real(rk) :: cuR1,cuR6
      real(rk) :: iNIn,iN1p, iN
      real(rk) :: ruR1c,ruR2c,ruR3c,ruR6c
      real(rk) :: ruX1c,ruX2c,ruX3c
      real(rk) :: rut,rug
      real(rk) :: rrc
      real(rk) :: flN6rPBA
      real(rk) :: ruR1n,ruR6n
      real(rk) :: ruR1p,ruR6p
      real(rk) :: huln,ren, rep, rep1
      real(rk) :: run
      real(rk) :: rump,hulp
      real(rk) :: rumn,rumn3,rumn4
      real(rk) :: r
      real(rk) :: reR2c, reR3c

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         ! Retrieve local biomass (carbon, phosphorus, nitrogen
         ! concentrations).

         ! Concentrations excluding background (used in sink terms)
         _GET_(self%id_c,bacc)
         _GET_(self%id_p,bacp)
         _GET_(self%id_n,bacn)

         ! Retrieve ambient nutrient concentrations
         _GET_(self%id_O2o,O2o)
         _GET_(self%id_N1p,N1p)
         _GET_(self%id_N3n,N3n)
         _GET_(self%id_N4n,N4n)

         ! Reduction equivalent
         _GET_(self%id_N6r,N6r)

         ! Retrieve ambient organic matter concentrations
         _GET_(self%id_R1c,R1c)
         _GET_(self%id_R1p,R1p)
         _GET_(self%id_R1n,R1n)

         _GET_(self%id_R2c,R2c)

         _GET_(self%id_R3c,R3c)

         _GET_(self%id_R6c,R6c)
         _GET_(self%id_R6p,R6p)
         _GET_(self%id_R6n,R6n)

         _GET_(self%id_X1c,X1c)
         _GET_(self%id_X2c,X2c)
         _GET_(self%id_X3c,X3c)

         ! Retrieve environmental dependencies (water temperature,
         ! photosynthetically active radation)
         _GET_(self%id_ETW,ETW)
  ! Quota collectors
         qpcPBA = bacp/(bacc+p_small) ! add some epsilon (add in shared) to avoid divide by 0
         qncPBA = bacn/(bacc+p_small)

         qpcR1  = R1p/(R1c+p_small)
         qncR1  = R1n/(R1c+p_small)

         qpcR6  = R6p/(R6c+p_small)
         qncR6  = R6n/(R6c+p_small)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et = eTq(ETW, self%p_q10)

 _SET_DIAGNOSTIC_(self%id_et,et)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Oxygen environment: bacteria are both aerobic and anaerobic
  ! To provide a faster switching between the two metabolic pathways the
  ! oxygen regulating factor eO2 is cubic
  ! (eq. 19 in Vichi et al., 2004)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2 = MM_power(max(p_small,O2o),  self%p_chdo,3)

 _SET_DIAGNOSTIC_(self%id_eO2,eO2)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! External nutrient limitation (used by some parametrizations)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eN4n = MM(N4n, self%p_chn)
  eN1p = MM(N1p, self%p_chp)

 _SET_DIAGNOSTIC_(self%id_eN4n,eN4n)
 _SET_DIAGNOSTIC_(self%id_eN1p,eN1p)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Mortality:
  !   1. first order mortality: p_sd 
  !   2. density dependent mortality due to virus infection: p_sd2
  !
  !   It is assumed that mortality is distributed in the same way over
  !   DOC (R1) and detritus (R6) as for phytoplankton and microzooplankton
  !   using the p_pe_R1x parameters defined in Param
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd  =  ( self%p_sd*et + self%p_sd2*bacc ) * bacc

 _SET_DIAGNOSTIC_(self%id_rd,rd)
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppbacc, ppR6c, rd*(ONE-p_pe_R1c)              , tfluxC)
  _SET_ODE_(self%id_R6c,rd*(ONE-self%p_pe_R1c))
  _SET_ODE_(self%id_c,-(rd*(ONE-self%p_pe_R1c)))
!SEAMLESS call quota_flux(iiPel, ppbacn, ppbacn, ppR6n, rd*qncPBA(bac,:)*(ONE-p_pe_R1n), tfluxN)
  _SET_ODE_(self%id_R6n,rd*qncPBA*(ONE-self%p_pe_R1n))
  _SET_ODE_(self%id_n,-(rd*qncPBA*(ONE-self%p_pe_R1n)))
!SEAMLESS call quota_flux(iiPel, ppbacp, ppbacp, ppR6p, rd*qpcPBA(bac,:)*(ONE-p_pe_R1p), tfluxP)
  _SET_ODE_(self%id_R6p,rd*qpcPBA*(ONE-self%p_pe_R1p))
  _SET_ODE_(self%id_p,-(rd*qpcPBA*(ONE-self%p_pe_R1p)))
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppbacc, ppR1c, 0.98D0*rd*p_pe_R1c , tfluxC) ! flux to non CDOM
  _SET_ODE_(self%id_R1c,(1.00D0-self%p_fX1b)*rd*self%p_pe_R1c)
  _SET_ODE_(self%id_c,-(1.00D0-self%p_fX1b)*rd*self%p_pe_R1c)
!SEAMLESS call quota_flux(iiPel, ppbacc, ppbacc, ppR1l, 0.02D0*rd*p_pe_R1c , tfluxC)  ! flux to CDOM
  _SET_ODE_(self%id_X1c,self%p_fX1b*rd*self%p_pe_R1c)
  _SET_ODE_(self%id_c,-self%p_fX1b*rd*self%p_pe_R1c)
!SEAMLESS  call quota_flux(iiPel, ppbacn, ppbacn, ppR1n, rd*qncPBA(bac,:)*p_pe_R1n, tfluxN) 
  _SET_ODE_(self%id_R1n,rd*qncPBA*self%p_pe_R1n)
  _SET_ODE_(self%id_n,-rd*qncPBA*self%p_pe_R1n)
!SEAMLESS  call quota_flux(iiPel, ppbacp, ppbacp, ppR1p, rd*qpcPBA(bac,:)*p_pe_R1p, tfluxP)
  _SET_ODE_(self%id_R1p,rd*qpcPBA*self%p_pe_R1p)
  _SET_ODE_(self%id_p,-rd*qpcPBA*self%p_pe_R1p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Substrate availability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( self%p_version )

    case ( BACT3 )  ! Polimene et al. (2006)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      ! Note: oxygen control in eq 5 (Polimene et al. 2006) is not included
      !        as bacteria are both aerobic and anaerobic
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  =  self%p_sum*et*bacc

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! No correction of organic material quality
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      cuR1 = ONE
      cuR6 = ONE

    case ( BACT1,BACT2 )  ! Vichi et al. (2004,2007), Lazzari et al. (2012) 

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
      ! Nutrient limitation (intracellular, eq. 51 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      iNIn = min(ONE, max(ZERO, qncPBA/self%p_qncPBA))  !Nitrogen
      iN1p = min(ONE, max(ZERO, qpcPBA/self%p_qpcPBA))  !Phosphorus
      iN   = min(iN1p, iNIn)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria (eq. 50 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  = iN*et*self%p_sum*bacc

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! correction of substrate quality depending on nutrient content
      ! (eq. 52 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      cuR1 = min(ONE, qpcR1/self%p_qpcPBA, qncR1/self%p_qncPBA)
      cuR6 = min(ONE, qpcR6/self%p_qpcPBA, qncR6/self%p_qncPBA)

  end select

 _SET_DIAGNOSTIC_(self%id_iN1p,iN1p)
 _SET_DIAGNOSTIC_(self%id_iNIn,iNIn)
 _SET_DIAGNOSTIC_(self%id_iN,iN)
 _SET_DIAGNOSTIC_(self%id_rum,rum)
 _SET_DIAGNOSTIC_(self%id_cuR1,cuR1)
 _SET_DIAGNOSTIC_(self%id_cuR6,cuR6)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate the realized substrate uptake rate depending on the
  ! type of detritus and quality (cuRx)
  ! See eq 27 in Vichi et al., 2004 for R2
  ! and eq 6 in Polimene et al., 2006 for R3 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ruR1c = (self%p_suhR1*cuR1 + self%p_sulR1*(ONE-cuR1))*R1c
  ruR2c = self%p_suR2*R2c
  ruR3c = self%p_suR3*R3c
  ruR6c = self%p_suR6*cuR6*R6c

  ruX1c = (self%p_suhR1*cuR1 + self%p_sulR1*(ONE-cuR1))*X1c
  ruX2c = self%p_suR2*X2c
  ruX3c = self%p_suR3*X3c

  rut   = p_small + ruR1c + ruX1c + ruR2c + ruX2c + ruR6c + ruR3c + ruX3c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria (eq. 50 Vichi et al. 2007)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rug = min( rum, rut )

 _SET_DIAGNOSTIC_(self%id_rug,rug)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR1c = rug*ruR1c/rut
  ruR2c = rug*ruR2c/rut
  ruR3c = rug*ruR3c/rut
  ruR6c = rug*ruR6c/rut

  ruX1c = rug*ruX1c/rut 
  ruX2c = rug*ruX2c/rut
  ruX3c = rug*ruX3c/rut

!SEAMLESS  call quota_flux(iiPel, ppbacc, ppR1c, ppbacc, ruR1c, tfluxC)
  _SET_ODE_(self%id_c,ruR1c)
  _SET_ODE_(self%id_R1c,-ruR1c)
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppR2c, ppbacc, ruR2c, tfluxC)
  _SET_ODE_(self%id_c,ruR2c)
  _SET_ODE_(self%id_R2c,-ruR2c)
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppR3c, ppbacc, ruR3c, tfluxC)
  _SET_ODE_(self%id_c,ruR3c)
  _SET_ODE_(self%id_R3c,-ruR3c)
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppR6c, ppbacc, ruR6c, tfluxC)
  _SET_ODE_(self%id_c,ruR6c)
  _SET_ODE_(self%id_R6c,-ruR6c)
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppR1l, ppbacc, ruR1l, tfluxC)
  _SET_ODE_(self%id_c,ruX1c)
  _SET_ODE_(self%id_X1c,-ruX1c)
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppR2l, ppbacc, ruR2l, tfluxC)
  _SET_ODE_(self%id_c,ruX2c)
  _SET_ODE_(self%id_X2c,-ruX2c)
!SEAMLESS  call quota_flux(iiPel, ppbacc, ppR3l, ppbacc, ruR3l, tfluxC)
  _SET_ODE_(self%id_c,ruX3c)
  _SET_ODE_(self%id_X3c,-ruX3c)

 _SET_DIAGNOSTIC_(self%id_ruR1c,ruR1c)
 _SET_DIAGNOSTIC_(self%id_ruR2c,ruR2c)
 _SET_DIAGNOSTIC_(self%id_ruR3c,ruR3c)
 _SET_DIAGNOSTIC_(self%id_ruR6c,ruR6c)

 _SET_DIAGNOSTIC_(self%id_ruX1c,ruX1c)
 _SET_DIAGNOSTIC_(self%id_ruX2c,ruX2c)
 _SET_DIAGNOSTIC_(self%id_ruX3c,ruX3c)

 _SET_DIAGNOSTIC_(self%id_rut,rut)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake CDOM uptake included
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR1n = qncR1*ruR1c
  ruR6n = qncR6*ruR6c

!SEAMLESS  call quota_flux(iiPel, ppbacn, ppR1n, ppbacn, ruR1n, tfluxN)
  _SET_ODE_(self%id_n,ruR1n)
  _SET_ODE_(self%id_R1n,-ruR1n)
!SEAMLESS  call quota_flux(iiPel, ppbacn, ppR6n, ppbacn, ruR6n, tfluxN)
  _SET_ODE_(self%id_n,ruR6n)
  _SET_ODE_(self%id_R6n,-ruR6n)

  ruR1p = qpcR1*ruR1c
  ruR6p = qpcR6*ruR6c

!SEAMLESS  call quota_flux(iiPel, ppbacp, ppR1p, ppbacp, ruR1p, tfluxP)
  _SET_ODE_(self%id_p,ruR1p)
  _SET_ODE_(self%id_R1p,-ruR1p)
!SEAMLESS  call quota_flux(iiPel, ppbacp, ppR6p, ppbacp, ruR6p, tfluxP)
  _SET_ODE_(self%id_p,ruR6p)
  _SET_ODE_(self%id_R6p,-ruR6p)

 _SET_DIAGNOSTIC_(self%id_ruR1n,ruR1n)
 _SET_DIAGNOSTIC_(self%id_ruR6n,ruR6n)

 _SET_DIAGNOSTIC_(self%id_ruR1p,ruR1p)
 _SET_DIAGNOSTIC_(self%id_ruR6p,ruR6p)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Aerobic and anaerobic respiration 
  ! Pelagic bacteria are a wide functional group comprising both aerobic and
  ! anaerobic bacteria. At (very) low Oxygen concentrations bacteria use
  ! N6r as electron acceptor in the respiration process. 
  ! However, the carbon cost is higher and an additional term is used.
  ! If nitrate is present, the rate of consumption of N6r is converted to N3n
  ! consumption (eq 19 Vichi et al., 2004 and PelChem.F90)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrc = (self%p_pu_ra+ self%p_pu_ra_o*(ONE-eO2) )*rug + self%p_srs* bacc* et
!SEAMLESS  call quota_flux( iiPel, ppbacc, ppbacc, ppO3c, rrc, tfluxC) 
  _SET_ODE_(self%id_O3c,rrc)
  _SET_ODE_(self%id_c,-rrc)
!SEAMLESS  call flux_vector( iiPel, ppO2o, ppO2o, -eO2*rrc/MW_C )
  _SET_ODE_(self%id_O2o,-eO2*rrc/MW_C )

  flN6rPBA = (ONE- eO2)*rrc/ MW_C* self%p_qro
!SEAMLESS  call flux_vector( iiPel, ppN6r, ppN6r, flN6rPBA )
  _SET_ODE_(self%id_N6r, flN6rPBA )

  ! Update the total rate of formation of reduction equivalent
  flPTN6r = flN6rPBA
! flPTN6r = flPTN6r + flN6rPBA

 _SET_DIAGNOSTIC_(self%id_rrc,rrc)
 _SET_DIAGNOSTIC_(self%id_flN6rPBA, flN6rPBA)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !             Fluxes from bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case (self%p_version)

    case ( BACT1 ) ! Vichi et al. 2007
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! There is no Carbon excretion in Vichi et al. 2007
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
      ! This rate is assumed to occur with a timescale p_ruen=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ren  =  (qncPBA - self%p_qncPBA)*bacc*self%p_ruen
!SEAMLESS      call quota_flux(iiPel, ppbacn, ppbacn, ppN4n,       ren*insw( ren), tfluxN)
  _SET_ODE_(self%id_N4n,ren*insw( ren))
  _SET_ODE_(self%id_n, -ren*insw( ren))
!SEAMLESS      call quota_flux(iiPel, ppbacn, ppN4n, ppbacn, -eN4n*ren*insw(-ren), tfluxN)
  _SET_ODE_(self%id_n,   -eN4n*ren*insw(-ren))
  _SET_ODE_(self%id_N4n, -(-eN4n*ren*insw(-ren)))

 _SET_DIAGNOSTIC_(self%id_ren,ren)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
      ! This rate is assumed to occur with a timescale of p_ruep=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rep  =  (qpcPBA - self%p_qpcPBA)*bacc*self%p_ruep
!SEAMLESS      call quota_flux(iiPel, ppbacp, ppbacp, ppN1p,       rep*insw( rep), tfluxP)
  _SET_ODE_(self%id_N1p,rep*insw( rep))
  _SET_ODE_(self%id_p, -rep*insw( rep))
!SEAMLESS      call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -eN1p*rep*insw(-rep), tfluxP)
  _SET_ODE_(self%id_p, -eN1p*rep*insw(-rep))
  _SET_ODE_(self%id_N1p, -(eN1p*rep*insw(-rep)))

 _SET_DIAGNOSTIC_(self%id_rep,rep)

    case ( BACT2 ) ! Vichi et al. 2004
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Net Production
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      run = rug - rrc
  
 _SET_DIAGNOSTIC_(self%id_run,run)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium remineralization  (Eq. 28 Vichi et al. 2004, note that there 
      ! is a bug in the paper as there should be no division by B1c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      huln = (ruR6n + ruR1n) - self%p_qncPBA*run
      ren  = huln*insw(huln)
!SEAMLESS      call quota_flux(iiPel, ppbacn, ppbacn, ppN4n, ren, tfluxN)
  _SET_ODE_(self%id_N4n, ren)
  _SET_ODE_(self%id_n,  -ren)

 _SET_DIAGNOSTIC_(self%id_renh4,ren)
 _SET_ODE_(self%id_O3h,ren)
 _SET_DIAGNOSTIC_(self%id_varO3h_reminNH4,ren)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Nitrogen uptake  (Eq. 29 Vichi et al. 2004, there is a bug
      ! here as well, the min should be a max as the numbers are negative!)
      ! from Ammonium and Nitrate when N is not balanced (huln<0)
      ! (nitrate uptake with ammonium inhibition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rumn3 = self%p_qun*N3n*bacc*(ONE-eN4n)
      rumn4 = self%p_qun*N4n*bacc
      rumn  = rumn3 + rumn4
      ren   = max(-rumn,huln)*insw(-huln)
!SEAMLESS      call quota_flux(iiPel, ppbacn, ppN4n, ppbacn, -ren*rumn4/rumn, tfluxN)
  _SET_ODE_(self%id_n, -ren*rumn4/rumn)
  _SET_ODE_(self%id_N4n, -(-ren*rumn4/rumn))
!SEAMLESS      call quota_flux(iiPel, ppbacn, ppN3n, ppbacn, -ren*rumn3/rumn, tfluxN)
  _SET_ODE_(self%id_n, -ren*rumn3/rumn)
  _SET_ODE_(self%id_N3n, -(-ren*rumn3/rumn))

 _SET_DIAGNOSTIC_(self%id_rumn3,rumn3)
 _SET_DIAGNOSTIC_(self%id_rumn4,rumn4)
 _SET_DIAGNOSTIC_(self%id_rumn,rumn)
 _SET_DIAGNOSTIC_(self%id_ren,ren)

! uptake of no3 -> increase alk and viceversa
! uptake of nh4 -> decrease alk and viceversa
 _SET_ODE_(self%id_O3h, (-ren*rumn3/rumn) + (ren*rumn4/rumn))
 _SET_DIAGNOSTIC_(self%id_B_N_O3h,(-ren*rumn3/rumn) + ren*rumn4/rumn)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Phosphate remineralization  (Eq. 28 Vichi et al. 2004, note that there 
      ! is an error in the paper as there should be no division by B1c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      hulp = (ruR6p + ruR1p) - self%p_qpcPBA*run
      rep  = hulp*insw(hulp)
!SEAMLESS      call quota_flux(iiPel, ppbacp, ppbacp, ppN1p, rep, tfluxP)
  _SET_ODE_(self%id_N1p, rep)
  _SET_ODE_(self%id_p, -rep)
  rep1=rep
 _SET_DIAGNOSTIC_(self%id_repo4,rep)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Phosphorus uptake  (Eq. 29 Vichi et al. 2004, there is a bug
      ! here as well, the min should be a max as the numbers are negative!)
      ! from dissolved phosphate whene P is not balanced (hulp<0)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rump = self%p_qup*N1p*bacc
      rep  = max(-rump,hulp)*insw(-hulp)
!SEAMLESS      call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -rep, tfluxP)
  _SET_ODE_(self%id_p, -rep)
  _SET_ODE_(self%id_N1p, -(-rep))

 _SET_DIAGNOSTIC_(self%id_rep,rep)  ! rep contiene valori negativi (uptake) 

! rep1 (remineraliz contiene valori positivi) => -ALK ; rep (uptake, ma la variabile contiene valori negativi) => +ALK
! _SET_ODE_(self%id_O3h,-(rep1) + -rep )
 _SET_DIAGNOSTIC_(self%id_B_P_O3h,-(rep1) + -rep)
  
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Excess carbon (also considering dissolved nutrient uptake ren and rep) 
      ! is released as R3c, no other excretion (reR2c=0)
      ! (eq. 30 Vichi et al. 2004, unfortunately there is another error in 
      ! the paper, the flux of dissolved nutrient is not written)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      r     = min(run, (ruR6n+ruR1n-ren)/self%p_qlnc)
      reR3c = run - min(r, (ruR6p+ruR1p-rep)/self%p_qlpc)
      reR3c = max(ZERO, reR3c)
!SEAMLESS      call quota_flux( iiPel, ppbacc, ppbacc,ppR3c, 0.98D0*reR3c ,tfluxC) ! flux to non-CDOM
  _SET_ODE_(self%id_R3c, (1.00D0-self%p_fX3b)*reR3c)
  _SET_ODE_(self%id_c,  -(1.00D0-self%p_fX3b)*reR3c)
!SEAMLESS      call quota_flux( iiPel, ppbacc, ppbacc,ppR3l, 0.02D0*reR3c ,tfluxC) ! flux to CDOM
  _SET_ODE_(self%id_X3c, self%p_fX3b*reR3c)
  _SET_ODE_(self%id_c,  -self%p_fX3b*reR3c)
 _SET_DIAGNOSTIC_(self%id_reR3c,reR3c)

    case ( BACT3 ) ! Polimene et al. (2006)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Carbon excretion as Semi-Labile (R2) and Semi-Refractory (R3) DOC
      ! The R2 rate is assumed to occur with a timescale of 1 day
      ! (eq 8 Polimene et al., 2006)
      ! The renewal of capsular material is a constant rate, equivalent
      ! to about 1/4 of the respiration rate, ~5% of uptake
      ! (Stoderegger and Herndl, 1998)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      reR2c = max((ONE-(qpcPBA/self%p_qpcPBA)), &
              (ONE-(qncPBA/self%p_qncPBA)))*self%p_rec
      reR2c = max(ZERO,reR2c)*bacc
      reR3c = rug*(ONE-self%p_pu_ra)*(self%p_pu_ra*self%p_pu_ea_R3)
 _SET_DIAGNOSTIC_(self%id_reR2c,reR2c)
 _SET_DIAGNOSTIC_(self%id_reR3c,reR3c)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
      ! This rate is assumed to occur with a timescale p_ruen=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ren = (qncPBA - self%p_qncPBA)*bacc*self%p_ruen
!SEAMLESS      call quota_flux(iiPel, ppbacn, ppbacn, ppN4n,       ren*insw( ren), tfluxN)
  _SET_ODE_(self%id_N4n, ren*insw( ren))
  _SET_ODE_(self%id_n,  -(ren*insw( ren)))
!SEAMLESS      call quota_flux(iiPel, ppbacn, ppN4n, ppbacn, -eN4n*ren*insw(-ren), tfluxN)
  _SET_ODE_(self%id_n, -eN4n*ren*insw(-ren))
  _SET_ODE_(self%id_N4n, -(-eN4n*ren*insw(-ren)))

 _SET_DIAGNOSTIC_(self%id_ren,ren)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
      ! This rate is assumed to occur with a timescale of p_ruep=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rep  =  (qpcPBA- self%p_qpcPBA)*bacc*self%p_ruep
!SEAMLESS      call quota_flux(iiPel, ppbacp, ppbacp, ppN1p,       rep*insw( rep), tfluxP)
  _SET_ODE_(self%id_N1p, rep*insw( rep))
  _SET_ODE_(self%id_p,  -rep*insw( rep))
!SEAMLESS      call quota_flux(iiPel, ppbacp, ppN1p, ppbacp, -eN1p*rep*insw(-rep), tfluxP)
  _SET_ODE_(self%id_p,  -eN1p*rep*insw(-rep))
  _SET_ODE_(self%id_N1p, -(-eN1p*rep*insw(-rep)))

  _SET_DIAGNOSTIC_(self%id_rep,rep)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Excretion fluxes (only losses to R2 and R3)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS      call quota_flux( iiPel, ppbacc, ppbacc, ppR2c, 0.98D0*reR2c, tfluxC) ! flux to non CDOM
  _SET_ODE_(self%id_R2c, (1.00D0-self%p_fX2b)*reR2c)
  _SET_ODE_(self%id_c,  -(1.00D0-self%p_fX2b)*reR2c)
!SEAMLESS      call quota_flux( iiPel, ppbacc, ppbacc, ppR2l, 0.02D0*reR2c, tfluxC) ! flux to CDOM
  _SET_ODE_(self%id_X2c, self%p_fX2b*reR2c)
  _SET_ODE_(self%id_c,  -self%p_fX2b*reR2c)

!SEAMLESS      call quota_flux( iiPel, ppbacc, ppbacc, ppR3c, 0.98D0*reR3c, tfluxC) ! flux to non CDOM
  _SET_ODE_(self%id_R3c, (1.00D0-self%p_fX3b)*reR3c)
  _SET_ODE_(self%id_c,  -(1.00D0-self%p_fX3b)*reR3c)
!SEAMLESS      call quota_flux( iiPel, ppbacc, ppbacc, ppR3l, 0.02D0*reR3c, tfluxC) ! flux to CDOM
  _SET_ODE_(self%id_X3c, self%p_fX3b*reR3c)
  _SET_ODE_(self%id_c,  -self%p_fX3b*reR3c)
  
  end select

      _LOOP_END_

   end subroutine do
end module
!SEAMLESS!EOC
!SEAMLESS!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS! MODEL  BFM - Biogeochemical Flux Model 
!SEAMLESS!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
