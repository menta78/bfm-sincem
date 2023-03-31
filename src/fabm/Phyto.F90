#include "fabm_driver.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   This process describes the dynamics of all phytoplankton
!    groups. The differences in behaviour
!    are expressed by differences in parameter-values only.
!    
! !INTERFACE
 module bfm_Phyto

   use fabm_types
   use ogs_bfm_shared
!  use fabm_particle

   use ogs_bfm_pelagic_base

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   implicit none

   private
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_primary_producer
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving
      ! from type_ogs_bfm_pelagic_base!

      ! Identifiers for state variables of other models
      type (type_state_variable_id) :: id_O3c,id_O2o,id_O3h                 !  dissolved inorganic carbon, oxygen, total alkalinity
      type (type_state_variable_id) :: id_N1p,id_N3n,id_N4n,id_N5s          !  nutrients: phosphate, nitrate, ammonium, silicate, iron
      type (type_state_variable_id) :: id_R1c,id_R1p,id_R1n,id_R2c          !  dissolved organic carbon (R1: labile, R2: semi-labile)
      type (type_state_variable_id) :: id_R6c,id_R6p,id_R6n,id_R6s          !  particulate organic carbon
      type (type_state_variable_id) :: id_X1c,id_X2c                        !  coloured dissolved organic carbon
      type (type_state_variable_id) :: id_O5c                               !  Free calcite (liths) - used by calcifiers only
      ! Environmental dependencies
      type (type_dependency_id)            :: id_ETW   ! PAR and temperature
!     type (type_dependency_id)            :: id_parEIR,id_ETW   ! PAR and temperature
!     type (type_dependency_id)            :: id_PAR_tot         ! scalar spectral PAR
!     type (type_dependency_id)            :: id_par_dia, id_par_flag, id_par_pico, id_par_dino
      type (type_dependency_id)            :: id_PAR
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_iN1p  ! internal quota phosphorus limitation 
      type (type_diagnostic_variable_id) :: id_iNIn  ! internal quota nitrogen limitation 
      type (type_diagnostic_variable_id) :: id_iN5s  ! internal quota silicon limitation 
      type (type_diagnostic_variable_id) :: id_eN5s  ! external silicate limitation 
      type (type_diagnostic_variable_id) :: id_iN    ! P and N limitation
      type (type_diagnostic_variable_id) :: id_tN    ! total nutrient limitation
      type (type_diagnostic_variable_id) :: id_ETWd  ! temperature Celsius
      type (type_diagnostic_variable_id) :: id_EIRd  ! light diagnostic
      type (type_diagnostic_variable_id) :: id_eT    ! temperature q10 factor
      type (type_diagnostic_variable_id) :: id_rr    ! light limitation exponent
      type (type_diagnostic_variable_id) :: id_eiPPY ! light limitation
      type (type_diagnostic_variable_id) :: id_sum   ! growth time scale
      type (type_diagnostic_variable_id) :: id_sdo   ! nutrient stress lysis
      type (type_diagnostic_variable_id) :: id_sea   ! activity excretion
      type (type_diagnostic_variable_id) :: id_seo   ! nutrient stress excrection
      type (type_diagnostic_variable_id) :: id_rr1c  ! lysis fraction to labile DOC
      type (type_diagnostic_variable_id) :: id_rrc   ! total respiration
      type (type_diagnostic_variable_id) :: id_rugc  ! gross primary production
      type (type_diagnostic_variable_id) :: id_flPIR2c  ! release to semi-labile transparent DOC
      type (type_diagnostic_variable_id) :: id_flPIR2c_act  ! activity release to semi-labile DOC
      type (type_diagnostic_variable_id) :: id_flPIR2c_tot  ! total release to semi-labile DOC      
      type (type_diagnostic_variable_id) :: id_f2cdom  ! fraction to CDOM    
      type (type_diagnostic_variable_id) :: id_run   ! net primary production
      type (type_diagnostic_variable_id) :: id_sadap  ! adaptation
      type (type_diagnostic_variable_id) :: id_cqun3  ! preference for ammonia 
      type (type_diagnostic_variable_id) :: id_rumn3  ! max pot. uptake of N3
      type (type_diagnostic_variable_id) :: id_rumn4  ! max pot. uptake of N4
      type (type_diagnostic_variable_id) :: id_rumn   ! max pot. uptake of DIN
      type (type_diagnostic_variable_id) :: id_rump   ! max pot. uptake of PO4
      type (type_diagnostic_variable_id) :: id_netgrowth ! netgrowth
      type (type_diagnostic_variable_id) :: id_sunPPY    ! Specific net growth rate 
      type (type_diagnostic_variable_id) :: id_misn      ! Intracellular missing amount of N
      type (type_diagnostic_variable_id) :: id_rupn      !  N uptake based on net assimilat. C
      type (type_diagnostic_variable_id) :: id_runn      !  actual uptake of NI
      type (type_diagnostic_variable_id) :: id_runn3     !  actual uptake of N3 
      type (type_diagnostic_variable_id) :: id_runn4     !  actual uptake of N4
      type (type_diagnostic_variable_id) :: id_fR1n      !  flux to R1n
      type (type_diagnostic_variable_id) :: id_misp      !  Intracellular missing amount of P
      type (type_diagnostic_variable_id) :: id_rupp      !  P uptake based on C uptake
      type (type_diagnostic_variable_id) :: id_runp      !  actual uptake
      type (type_diagnostic_variable_id) :: id_fR1p      !  flux to R1p
      type (type_diagnostic_variable_id) :: id_rr6n      !  Excretion to PON
      type (type_diagnostic_variable_id) :: id_rr1n      !  Excretion to DON
      type (type_diagnostic_variable_id) :: id_rr6p      !  Excretion to POP
      type (type_diagnostic_variable_id) :: id_rr1p      !  Excretion to DOP
      type (type_diagnostic_variable_id) :: id_rums      !  max pot. uptake of  S
      type (type_diagnostic_variable_id) :: id_miss      !  Intracellular missing amount of S
      type (type_diagnostic_variable_id) :: id_rups      !  S uptake based on C uptake
      type (type_diagnostic_variable_id) :: id_runs      !  actual uptake
      type (type_diagnostic_variable_id) :: id_rho_Chl   !  Chlorophyll to Carbon ration
      type (type_diagnostic_variable_id) :: id_rate_Chl  !  Chlorophyll production per unit of carbon
      type (type_diagnostic_variable_id) :: id_O3hconume_for_CaCO3prec !consume of alk for caco3 precipitation
      type (type_diagnostic_variable_id) :: id_Putil_O3h !  variation of O3h due to net utilization of P (uptake-release)
      type (type_diagnostic_variable_id) :: id_Nutil_O3h !  variation of O3h due to net utilization of N (uptake-release)
 
      ! Parameters (described in subroutine initialize, below)
      real(rk) :: p_q10,p_temp,p_sum,p_srs,p_sdmo,p_thdo,p_seo,p_sheo,p_pu_ea,p_pu_ra
      real(rk) :: p_qun,p_lN4, p_qnlc, p_qncPPY, p_xqn, p_qup, p_qplc,p_qpcPPY, p_xqp
      real(rk) :: p_chPs, p_Contois, p_qus, p_qslc ,p_qscPPY
      real(rk) :: p_esNI,p_res
      real(rk) :: p_caco3r
      real(rk) :: p_sdchl, p_alpha_chl, p_quantum_yield, p_qlcPPY, p_epsChla, p_tochl_relt,p_EpEk_or
      real(rk) :: p_iswLtyp, p_chELiPPY, p_clELiPPY, p_ruELiPPY,p_addepth
      real(rk) :: p_rPIm
      real(rk) :: p_fX1p, p_fX2p
      integer :: p_switchDOC, p_switchSi,p_limnut,p_switchChl,p_Esource
      logical :: use_Si,p_netgrowth
      logical :: use_CaCO3
   contains

      ! Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement
      procedure :: get_sinking_rate

   end type type_ogs_bfm_primary_producer

   ! Constants
   real(rk),parameter :: pippo = 0.002_rk

contains

  subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_primary_producer),intent(inout),target :: self
      integer,                              intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      real(rk) :: pippo1
!EOP
!-----------------------------------------------------------------------
!BOC
     ! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used
      ! by FABM (or its host)
      ! to present parameters to the user for configuration (e.g., through a
      ! GUI)
      call self%get_parameter(self%p_q10,   'p_q10',     '-', 'Characteristic Q10 coefficient')
      call self%get_parameter(self%p_temp,  'p_temp',    '-', ' Cut-off threshold for temperature factor')
      call self%get_parameter(self%p_sum,   'p_sum',     '1/d',        'Maximal productivity at 10 degrees C')
      call self%get_parameter(self%p_srs,   'p_srs',     '1/d', 'Respiration rate at 10 degrees C')
      call self%get_parameter(self%p_sdmo,  'p_sdmo',    '1/d',        ' Max.specific nutrient-stress lysis rate')
      call self%get_parameter(self%p_thdo,  'p_thdo',    '-',        'Half saturation constant for nutrient stress lysis')
      call self%get_parameter(self%p_seo,   'p_seo',     '1/d',      'Extra lysis rate (biomass density-dependent)')
      call self%get_parameter(self%p_sheo,  'p_sheo',    'mgC/3',    'Half saturation constant for extra lysis')
      call self%get_parameter(self%p_pu_ea, 'p_pu_ea',    '-',        'Excreted fraction of primary production')
      call self%get_parameter(self%p_pu_ra, 'p_pu_ra',    '-',        'Activity respiration fraction')
      call self%get_parameter(self%p_switchDOC,  'p_switchDOC',    '[1-3]', 'Switch for the type of DOC excretion')
!                             This choice must be consistent with bacteria
!                             1. All DOC is released as R1c (Vichi et al., 2007)
!                             2. Activity DOC is released as R2c (Vichi et al.,
!                             2004)
!                                (there is no nutrient-stress excretion)
!                             3. All DOC is released as R2c (Polimene et al.,
!                             2006)

!              --------- Nutrient parameters in phytoplankton -----------------
      call self%get_parameter(self%p_netgrowth,  'p_netgrowth',    '[T or F]','Logical switch for nutrient-limited growth')
!                             .T. nutrient-balanced growth (Vichi et al.2004)
!                             .F. nutrient-stress carbon excretion
!                               (Baretta-Bekker et al.1995 and Vichi et al.2007)

      call self%get_parameter(self%p_limnut,  'p_limnut',    '[0-2]', 'Switch for N-P co-limitation')
!                             0. Geometric mean
!                             1. Threshold (Liebig-like)
!                             2. Combined

!                   ---- N limitation control ----
      call self%get_parameter(self%p_qun,    'p_qun'   ,'m3/mgC/d', 'Membrane affinity for N')
      call self%get_parameter(self%p_lN4,    'p_lN4'   ,'mmolN/m3', 'Half saturation constant for NH4 uptake preference over NO3')
      call self%get_parameter(self%p_qnlc,   'p_qnlc'  ,'mmolN/mgC','Minimum quotum N:C ')
      call self%get_parameter(self%p_qncPPY, 'p_qncPPY','mmolN/mgC','reference quotum N:C')
      call self%get_parameter(self%p_xqn,    'p_xqn'   ,'-', 'Multiplication factor for luxury storage')
!                   ---- P limitation control ----
      call self%get_parameter(self%p_qup,    'p_qup'   ,'m3/mgC/d', 'Membrane affinity for P')
      call self%get_parameter(self%p_qplc,   'p_qplc'  ,   'mmolP/mgC', ' Minimum quotum P:C')
      call self%get_parameter(self%p_qpcPPY, 'p_qpcPPY', 'mmolP/mgC', 'reference quotum P:C')
      call self%get_parameter(self%p_xqp,    'p_xqp'   ,     '-',   'Multiplication factor for luxury storage')
!                   ---- Si limitation control ----
      call self%get_parameter(self%use_Si,   'use_Si','',          'use silicate',default=.false.)
      if (self%use_Si) then 
          call self%get_parameter(self%p_switchSi, 'p_switchSi',   '[1-2]',    'Switch for Silica limitation')
!                             1. Si limitation is controlled by external Si
!                             conc.
!                             2. Si limitation is controlled by internal quota

          call self%get_parameter(self%p_chPs,'p_chPs', 'mmolSi/m3', 'Half saturation conc. for dissolved Si limitation')
          call self%get_parameter(self%p_Contois,  'p_Contois', '>=0', ' If >0, use Contois formulation')
          call self%get_parameter(self%p_qus,      'p_qus',  'm3/mgC/d', 'membrane affinity for Si')
          call self%get_parameter(self%p_qslc,     'p_qslc', 'mmolSi/mgC','minimum quotum for Si:C')
          call self%get_parameter(self%p_qscPPY,   'p_qscPPY','mmolSi/mgC',  'reference quotum Si:C')
      endif
!                   ---- nutrient stressed sinking ----
      call self%get_parameter(self%p_esNI,  'p_esNI',       '-', 'Nutrient stress threshold for sinking')
      call self%get_parameter(self%p_res,   'p_res',       'm/d', 'Maximum Sinking velocity ')
!              --------- Chlorophyll parameters -----------
      call self%get_parameter(self%p_switchChl,     'p_switchChl',    '1-4',    'Switch for Chla-a synthesis')
      call self%get_parameter(self%p_sdchl,         'p_sdchl',        '1/d',    'Specific turnover rate for Chla')
      call self%get_parameter(self%p_alpha_chl,     'p_alpha_chl',    'mgC s m2/mgChl/uE','Initial slope of the P-E curve')
      call self%get_parameter(self%p_quantum_yield, 'p_quantum_yield','mgC/uE', 'Photochemical efficiency')
      call self%get_parameter(self%p_Esource,       'p_Esource',      '1-6',    'source of light for PP')
      
      call self%get_parameter(self%p_qlcPPY,   'p_qlcPPY',  'mgChla/mgC','reference quotum Chla:C')
      call self%get_parameter(self%p_epsChla,   'p_epsChla',  'm2/mgChla', 'Chla-specific extinction coefficient')
      call self%get_parameter(self%p_tochl_relt,   'p_tochl_relt',  '1/d', 'Relaxation rate towards maximum Chla:C')
      call self%get_parameter(self%p_EpEk_or,   'p_EpEk_or',  '-',    'Optimal value of E_PAR/E_K')
!              --------- Light parameters ERSEM-II -----------
      call self%get_parameter(self%p_iswLtyp,   'p_iswLtyp',  '0-6',    'Shape of the productivity function')
      call self%get_parameter(self%p_chELiPPY,   'p_chELiPPY',  'W/m2', 'Maximum Iopt')
      call self%get_parameter(self%p_clELiPPY,   'p_clELiPPY',  'W/m2', 'Minimum Iopt')
      call self%get_parameter(self%p_ruELiPPY,   'p_ruELiPPY',  '1/d',  'Maximum daily shift in Iopt (1/d)')
      call self%get_parameter(self%p_addepth,   'p_addepth',  'm', 'Adaptation depth. Meaningless with high-res models')
!              --------- Sinking parameters -----------
      call self%get_parameter(self%p_rPIm,   'p_rPIm',  'm/d', 'Phytoplankton background sinking rate')
!              --------- Calcite parameters only for P2 ------------
      call self%get_parameter(self%use_CaCO3,   'use_CaCO3','',          'use calcite',default=.false.)
      call self%get_parameter(self%p_caco3r,   'p_caco3r',  '-',  'reference PIC:POC rain ratio', default=0._rk)
!              --------- Flux partition CDOM parameters ------------
      call self%get_parameter(self%p_fX1p,   'p_fX1p',  '-',  'colored fraction in labile dissolved organic carbon', default=0.02_rk)
      call self%get_parameter(self%p_fX2p,   'p_fX2p',  '-',  'colored fraction in semi-labile dissolved organic carbon', default=0.02_rk)
      
! Register state variables (handled by type_bfm_pelagic_base)
      call self%initialize_bfm_base()
      call self%add_constituent('c',1.e-4_rk)
      call self%add_constituent('n',1.26e-6_rk)
      call self%add_constituent('p',4.288e-8_rk)
      call self%add_constituent('f',5.e-6_rk)  ! NB this does nothing if iron support is disabled.
      call self%add_constituent('chl',3.e-6_rk)
      if (self%use_Si) call self%add_constituent('s',1.e-6_rk)
!     call self%add_constituent('c',1.e-4_rk,   c0)
!     call self%add_constituent('n',1.26e-6_rk, c0*qnrpicX)
!     call self%add_constituent('p',4.288e-8_rk,c0*qprpicX)
!     call self%add_constituent('f',5.e-6_rk,   0.0_rk)  ! NB this does nothing if iron support is disabled.
!     call self%add_constituent('chl',3.e-6_rk, c0*self%phim)
!     if (self%use_Si) call self%add_constituent('s',1.e-6_rk,c0*self%qsc)
      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_O3c,'O3c','mg C/m^3','dissolved organic carbon')
      call self%register_state_dependency(self%id_O3h,'O3h','mmol/m^3','alkalinity')
      call self%register_state_dependency(self%id_O5c,'O5c','mg C/m^3','particulate inorganic carbon')
      call self%register_state_dependency(self%id_O2o,'O2o','mmolO2/m^3','dissolved oxygen')
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      call self%register_state_dependency(self%id_N3n,'N3n','mmol N/m^3','nitrate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')
      if (self%use_Si) call self%register_state_dependency(self%id_N5s,'N5s','mmol Si/m^3','silicate')
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3','labile DOC')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','labile DOP')
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','labile DON')
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','semi labile DOC')
      call self%register_state_dependency(self%id_R6c,'R6c','mg C/m^3','POC')
      call self%register_state_dependency(self%id_R6p,'R6p','mmol P/m^3','POP')
      call self%register_state_dependency(self%id_R6n,'R6n','mmol N/m^3','PON')
      if (self%use_Si) call self%register_state_dependency(self%id_R6s,'R6s','mmol Si/m^3','POS')
      call self%register_state_dependency(self%id_X1c,'X1c','mg C/m^3','labile CDOM')
      call self%register_state_dependency(self%id_X2c,'X2c','mg C/m^3','semilabile CDOM')
      ! Register environmental dependencies (temperature, shortwave radiation)
!     call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      ! Dependency from multispectral model
!        select case (self%p_Esource)
!        case (1)
!           call self%register_dependency(self%id_par_dia,type_bulk_standard_variable(name='PAR_dia'))
!           call self%register_dependency(self%id_par,type_bulk_standard_variable(name='PAR_dia'))
!        case (2)
!           call self%register_dependency(self%id_par_flag,type_bulk_standard_variable(name='PAR_flag'))
!        case (3)
!           call self%register_dependency(self%id_par_pico,type_bulk_standard_variable(name='PAR_pico'))
!        case (4)
!            call self%register_dependency(self%id_par_dino,type_bulk_standard_variable(name='PAR_dino'))
!        case (5)
!            call self%register_dependency(self%id_PAR_tot,type_bulk_standard_variable(name='PAR_tot'))
!        case (6)
!            call self%register_dependency(self%id_PAR_tot,type_bulk_standard_variable(name='PAR'))
!           write(*,*) 'Use monochromatic light'
!        end select
         call self%register_dependency(self%id_PAR, 'PAR', '?????', 'effective PAR')

      
      ! Register diagnostic variables (i.e., model outputs)
      call self%register_diagnostic_variable(self%id_iN1p, 'iN1p', '-','internal quota phosphorus limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_iNIn, 'iNIn', '-','internal quota nitrogen limitation',output=output_none)
      if (self%use_Si) then 
         select case (self%p_switchSi) 
             case (1)  ! external control
                 call self%register_diagnostic_variable(self%id_eN5s, 'eN5s', '-','external silicate limitation',output=output_none)
             case(2)
                 call self%register_diagnostic_variable(self%id_iN5s, 'iN5s', '-','internal quota silicon limitation',output=output_none)
         end select
      end if
      call self%register_diagnostic_variable(self%id_iN,   'iN'  , '-','N and P  nutrient limitation')
      call self%register_diagnostic_variable(self%id_tN,   'tN',   '-','total nutrient limitation')
      call self%register_diagnostic_variable(self%id_ETWd, 'ETW',  'C','temperature Celsius',output=output_none)
      call self%register_diagnostic_variable(self%id_et,   'et',   '-','temperature factor',output=output_none)
      call self%register_diagnostic_variable(self%id_EIRd, 'EIR',  'uE/m2/s','PAR')
      call self%register_diagnostic_variable(self%id_rr,   'r',    '-',  'light limitation exponent')
      call self%register_diagnostic_variable(self%id_eiPPY,'eiPPY','-',  'light limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_sum,  'sum',  '1/d','growth time scale')
      call self%register_diagnostic_variable(self%id_sdo,  'sdo', 'mgC/m3/d','nutrient stress lysis')
      call self%register_diagnostic_variable(self%id_sea,  'sea', 'mgC/m3/d','activity excretion')
      call self%register_diagnostic_variable(self%id_seo,  'seo', 'mgC/m3/d','nutrient stress excretion')
      call self%register_diagnostic_variable(self%id_rr1c, 'rr1c','mgC/m3/d','lysis fraction to labile DOC')
      call self%register_diagnostic_variable(self%id_rrc,  'rrc', 'mgC/m3/d','total respiration')
      call self%register_diagnostic_variable(self%id_rugc, 'rugc','mgC/m3/d','Gross primary production')
      call self%register_diagnostic_variable(self%id_flPIR2c_tot,'flPIR2c_tot', 'mgC/m3/d', 'total flux to semilabile DOC')      
      call self%register_diagnostic_variable(self%id_flPIR2c_act,'flPIR2c_act', 'mgC/m3/d', 'activity flux to semilabile DOC')
      call self%register_diagnostic_variable(self%id_flPIR2c,    'flPIR2c',     'mgC/m3/d', 'flux to transparent semilabile DOC')
      call self%register_diagnostic_variable(self%id_f2cdom,  'f2cdom', '-', 'fraction to semilabile CDOM')      
      call self%register_diagnostic_variable(self%id_run,   'run',   'mgC/m3/d','net primary production')
      call self%register_diagnostic_variable(self%id_sadap, 'sadap', 'mgC/m3/d',' adaptation',output=output_none)
      call self%register_diagnostic_variable(self%id_cqun3, 'cqun3', '-',' preference for ammonia',output=output_none)
      call self%register_diagnostic_variable(self%id_rumn3, 'rumn3', '-',' max pot. uptake of N3',output=output_none)
      call self%register_diagnostic_variable(self%id_rumn4, 'rumn4', '-',' max pot. uptake of N4',output=output_none)
      call self%register_diagnostic_variable(self%id_rumn,  'rumn', '?',' max pot. uptake of DIN',output=output_none)
      call self%register_diagnostic_variable(self%id_rump,  'rump', '?',' max pot. uptake of DIN',output=output_none)
      call self%register_diagnostic_variable(self%id_netgrowth,  'netgrowth', 'mgC/m3/d',' netgrowth')
      call self%register_diagnostic_variable(self%id_sunPPY,  'sunPPY', '?',' Specific net growth rate',output=output_none)
      call self%register_diagnostic_variable(self%id_misn,  'misn', '?',' Intracellular missing amount of N',output=output_none)
      call self%register_diagnostic_variable(self%id_rupn,  'rupn', '?',' N uptake based on net assimilat. C',output=output_none)
      call self%register_diagnostic_variable(self%id_runn,  'runn', '?',' actual uptake of NI',output=output_none)
      call self%register_diagnostic_variable(self%id_runn3, 'runn3', '?',' actual uptake of N3',output=output_none)
      call self%register_diagnostic_variable(self%id_runn4, 'runn4', '?',' actual uptake of N4',output=output_none)
      call self%register_diagnostic_variable(self%id_fR1n,  'fR1n', '?',' flux to R1n',output=output_none)

      call self%register_diagnostic_variable(self%id_misp,  'misp', '?',' Intracellular missing amount of P',output=output_none)
      call self%register_diagnostic_variable(self%id_rupp,  'rupp', '?',' P uptake based on C uptake',output=output_none)
      call self%register_diagnostic_variable(self%id_runp,  'runp', '?',' actual P uptake ',output=output_none)
      call self%register_diagnostic_variable(self%id_fR1p,  'fR1p', '?',' flux to R1p  ',output=output_none)
      call self%register_diagnostic_variable(self%id_rr6n,  'rr6n', '?',' Excretion to PON  ',output=output_none)
      call self%register_diagnostic_variable(self%id_rr1n,  'rr1n', '?',' Excretion to DON  ',output=output_none)
      call self%register_diagnostic_variable(self%id_rr6p,  'rr6p', '?',' Excretion to POP  ',output=output_none)
      call self%register_diagnostic_variable(self%id_rr1p,  'rr1p', '?',' Excretion to DOP  ',output=output_none)
      if (self%use_Si) then 
         call self%register_diagnostic_variable(self%id_rums,  'rums', '?',' max pot. uptake of S',output=output_none)
         call self%register_diagnostic_variable(self%id_miss,  'miss', '?',' Intracellular missing amount of S',output=output_none)
         call self%register_diagnostic_variable(self%id_rups,  'rups', '?',' S uptake based on C uptake',output=output_none)
         call self%register_diagnostic_variable(self%id_runs,  'runs', '?',' actual S uptake ',output=output_none)
      endif
      call self%register_diagnostic_variable(self%id_rho_Chl,  'rho_Chl', 'mgChl/mgC','Chlorophyll production per unit of carbon ')
      call self%register_diagnostic_variable(self%id_rate_Chl,  'rate_Chl', 'mgChl/m3/d',' Chlorophyll production ')
      if (self%use_CaCO3) then
         call self%register_diagnostic_variable(self%id_O3hconume_for_CaCO3prec,'consO3h_caco3','mmol/m3/d','consume of O3h for CaCO3 precipitation',output=output_none)
      endif
         call self%register_diagnostic_variable(self%id_Nutil_O3h,'varO3h_for_Nutil','mmol/m3/d','variation of O3h due to N uptake/release',output=output_none)
         call self%register_diagnostic_variable(self%id_Putil_O3h,'varO3h_for_Putil','mmol/m3/d','variation of O3h due to P uptake/release',output=output_none)
   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ogs_bfm_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

   ! !LOCAL VARIABLES:
      real(rk) :: ETW,et, parEIR
      real(rk) :: phytoc, phytop, phyton, phytol, phytos
      real(rk) :: N5s,N1p,N3n,N4n,O2o
      real(rk) :: R1c,R1n,R1p
      real(rk) :: R2c
      real(rk) :: R6c,R6p,R6n,R6s
      real(rk) :: X1c,X2c
      real(rk) :: iNIn,iN1p,eN5s,iN5s,iNf,iNI
      real(rk) :: iN,tN
      real(rk) :: qpcPPY,qncPPY,qlcPPY,qscPPY
      real(rk) :: fpplim
      real(rk) :: r
      real(rk) :: eiPPY,photochem
      real(rk) :: sum
      real(rk) :: sdo, sea, seo
      real(rk) :: pe_R6, rr1c, rr6c
      real(rk) :: sra, srs, srt, rrc
      real(rk) :: rugc, slc, flPIR2c, flPIR2c_tot, f2cdom
      real(rk) :: run, sadap 
      real(rk) :: cqun3, rumn3, rumn4, rumn, rump
      real(rk) :: netgrowth, sunPPY  
      real(rk) :: tmp
      real(rk) :: misn, rupn, runn,  runn3, runn4, fR1n
      real(rk) :: misp, rupp, runp,  fR1p  
      real(rk) :: rr6n, rr1n, rr6p, rr1p
      real(rk) :: rums, miss, rups, runs
      real(rk) :: rho_Chl, rate_Chl, chl_opt

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         ! Retrieve local biomass (carbon, phosphorus, nitrogen, chlorophyll
         ! concentrations).

         ! Concentrations excluding background (used in sink terms)
         _GET_(self%id_c,phytoc)
         _GET_(self%id_p,phytop)
         _GET_(self%id_n,phyton)
         _GET_(self%id_chl,phytol)
         if (self%use_Si) then
            _GET_(self%id_s,phytos)
         endif

         ! Retrieve ambient nutrient concentrations
         _GET_(self%id_N1p,N1p)
         _GET_(self%id_N3n,N3n)
         _GET_(self%id_N4n,N4n)

         ! Retrieve ambient oxygen concentrations
         _GET_(self%id_O2o,O2o)

         ! Retrieve environmental dependencies (water temperature,
         ! photosynthetically active radation)
         _GET_(self%id_ETW,ETW)
         ! From where to get the light
!        select case (self%p_Esource)
!        case (1)
!           _GET_(self%id_par_dia,  parEIR)   ! uE mgChl-1 d-1
!        case (2)
!           _GET_(self%id_par_flag, parEIR)   ! uE mgChl-1 d-1
!        case (3)
!           _GET_(self%id_par_pico, parEIR)   ! uE mgChl-1 d-1       
!        case (4)
!           _GET_(self%id_par_dino, parEIR)   ! uE mgChl-1 d-1
!        case (5)
!           _GET_(self%id_PAR_tot,  parEIR)   ! uE m-2 d-1
!        case (6)
!           _GET_(self%id_parEIR,   parEIR)   ! uE m-2 d-1
!        end select

            _GET_(self%id_PAR,   parEIR)   ! uE m-2 d-1
         
  ! Quota collectors
         qpcPPY = phytop/(phytoc+p_small) ! add some epsilon (add in shared) to avoid divide by 0
         qncPPY = phyton/(phytoc+p_small)
         qlcPPY = phytol/(phytoc+p_small)
         if (self%use_Si) then
            qscPPY=phytos/(phytoc+p_small)
         endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitations (intracellular and extracellular)
  ! fpplim is the combined non-dimensional factor limiting photosynthesis
  ! Note for silicate limitation:
  !  p_switchSi =1 : external regulation of silica limitation 
  !  p_switchSi =2 : internal regulation of silica limitation 
  ! The standard Michaelis-Menten formulation contains the Contois parameter
  ! p_Contois=0: standard Michaelis Menten Formulation
  ! 0<p_Contois<=1: The Contois formulation is active. 
  !                 The limiting role of the population size (intraspecific 
  !                 competition) can be tuned by increasing p_Contois 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!        write(*,*) 'ONE', ONE
!        write(*,*) 'p_small', p_small
!        write(*,*) 'qpcPPY', qpcPPY
!        write(*,*) 'self%p_qplc', self%p_qplc
!        write(*,*) 'self%p_qpcPPY', self%p_qpcPPY
  iN1p = min( ONE, max( p_small, ( qpcPPY &
         - self%p_qplc)/( self%p_qpcPPY- self%p_qplc)))
!        write(*,*) 'iN1p', iN1p
  iNIn = min( ONE, max( p_small, ( qncPPY &
         - self%p_qnlc)/( self%p_qncPPY- self%p_qnlc)))
  if (self%use_Si) then
!    _GET_WITH_BACKGROUND_(self%id_N5s,N5s)
     _GET_(self%id_N5s,N5s)
     select case (self%p_switchSi) 
       case (1)  ! external control
         eN5s = min( ONE, N5s/(N5s + self%p_chPs+(self%p_Contois*phytos)))
         fpplim = eN5s
         iN5s   =  ONE
         _SET_DIAGNOSTIC_(self%id_eN5s,eN5s)

       case (2) ! internal control
         iN5s = min(ONE, max( p_small, ( qscPPY &
                - self%p_qslc)/( self%p_qscPPY- self%p_qslc)))
         fpplim = iN5s
         eN5s   = ONE
         _SET_DIAGNOSTIC_(self%id_iN5s,iN5s)

     end select
  else 
     iN5s   = ONE
     eN5s   = ONE
     fpplim = ONE
  end if
      _SET_DIAGNOSTIC_(self%id_iN1p,iN1p)
      _SET_DIAGNOSTIC_(self%id_iNIn,iNIn)
!SEAMLESS#ifdef INCLUDE_PELFE
!SEAMLESS  if (ppphytof > 0) then
!SEAMLESS     iN7f = min( ONE, max( p_small, ( qfcPPY(phyto,:) &
!SEAMLESS            - p_qflc(phyto))/( p_qfcPPY(phyto)- p_qflc(phyto))))
!SEAMLESS     fpplim = fpplim*iN7f
!SEAMLESS  else 
!SEAMLESS     iN7f   = ONE  
!SEAMLESS  end if
!SEAMLESS#endif
!SEAMLESS
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Multiple nutrient limitation
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
select case ( self%p_limnut)
    case ( 0 )
      iN  =   (iN1p* iNIn)**(0.5_rk)  ! geometric mean

    case ( 1 )
      iN  =   min(  iN1p,  iNIn)  ! Liebig rule

    case ( 2 )
      iN  =   2.0_rk/( ONE/ iN1p+ ONE/ iNIn)  ! combined

 end select

! tN only controls sedimentation of phytoplankton (Liebig)
  tN= min(iN,fpplim)

    _SET_DIAGNOSTIC_(self%id_iN,iN)
    _SET_DIAGNOSTIC_(self%id_tN,tN)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Phytoplankton
  ! Include cut-off at low temperature if p_temp>0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Retrieve environmental dependencies (water temperature, photosynthetically active radation)

  et  =   eTq(  ETW, self%p_q10)
  et  =   max(ZERO,et - self%p_temp)

    _SET_DIAGNOSTIC_(self%id_ETWd,ETW)
    _SET_DIAGNOSTIC_(self%id_et,et)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Light limitation with Chl dynamics
  ! If Chl is a diagnostic variable the limiting factor has been 
  ! computed in Light/PhotoAvailableRadiation.F90
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  _SET_DIAGNOSTIC_(self%id_EIRd,parEIR)

! Irradiance EIR is in uE m-2 d-1, 
! Irr is  average irradiance in uE m-2 day-1
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Compute exponent E_PAR/E_K = alpha0/PBmax
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     if (self%p_Esource<5) then
       photochem = self%p_quantum_yield
     else
       photochem = self%p_alpha_chl
     end if
  r = qlcPPY*photochem/self%p_sum * parEIR
  
select case ( LightPeriodFlag)
   case ( 1 ) ! instantaneous light
! no other factors needed
   case ( 2 ) ! daylight average is used
! recompute r and photsynthesis limitation using daylight scaling
      fpplim  =   fpplim*SUNQ/HOURS_PER_DAY
      r = r*HOURS_PER_DAY/SUNQ
   case ( 3 ) ! on-off
!     fpplim  =   fpplim*ThereIsLight
end select
  _SET_DIAGNOSTIC_(self%id_rr,r)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Light limitation factor according to Platt
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   eiPPY = ( ONE- exp( - r))

  _SET_DIAGNOSTIC_(self%id_eiPPY,eiPPY)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Total photosynthesis
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 sum  =   self%p_sum*et*eiPPY*fpplim

  _SET_DIAGNOSTIC_(self%id_sum,sum)

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Lysis and excretion
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 sdo  =  ( self%p_thdo/( iN + self%p_thdo))* self%p_sdmo  ! nutr. -stress lysis
! extra lysis for high-density
 sdo  =   sdo+ self%p_seo* MM(phytoc, self%p_sheo)

 sea  =   sum* self%p_pu_ea  ! activity excretion

 if (self%p_netgrowth) then
     seo = ZERO
 else 
! nutrient stress excretion
     seo = sum*(ONE-self%p_pu_ea)*(ONE- iN) 
 end if

 _SET_DIAGNOSTIC_(self%id_sdo,sdo)
 _SET_DIAGNOSTIC_(self%id_sea,sea)
 _SET_DIAGNOSTIC_(self%id_seo,seo)

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Apportioning over R1 and R6:
 ! Cell lysis generates both DOM and POM.
 ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
 ! Assuming that this structural part is not easily degradable,
 ! at least a fraction equal to the minimum quota is released as POM.
 ! Therefore, nutrients (and C) in the structural part go to R6.
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pe_R6 = min( self%p_qplc/( qpcPPY + p_small), self%p_qnlc/ &
          ( qncPPY+ p_small))
  pe_R6 = min(  ONE,  pe_R6)
  rr6c  =     pe_R6     * sdo * phytoc
  rr1c  = (ONE - pe_R6) * sdo * phytoc

 _SET_DIAGNOSTIC_(self%id_rr1c,rr1c)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Respiration rate
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 sra  =   self%p_pu_ra*( sum - sea - seo)  ! activity
 srs  =   et* self%p_srs                   ! basal
 srt  =   sra+ srs                         ! total
 rrc  =   srt* phytoc                      ! total actual respiration

 _SET_DIAGNOSTIC_(self%id_rrc,rrc)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Production, productivity and C flows
 ! The release of DOC is controlled by a specific switch.
 ! Beware that this switch must be consistent with the utilization of DOC 
 ! by Bacteria. If DOC is released in a form that is not used by 
 ! Bacteria, it will accumulate infinitely removing carbon from the system
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 rugc  =   sum* phytoc  ! gross production
 slc   =   sea + seo + srt+ sdo  ! specific loss terms
select case (self%p_switchDOC)
    case (1)
! All activity excretions are assigned to R1
        rr1c = rr1c + sea*phytoc + seo*phytoc
        flPIR2c = ZERO
    case (2)
! Activity excretion is only assigned to R2
        flPIR2c = sea* phytoc
!       BFM1D_exR2ac(phyto,:) = sea* phytoc
    case (3)
! Activity and Nutrient-stress excretions are assigned to R2
        flPIR2c  =  seo*phytoc + sea*phytoc
end select
 _SET_DIAGNOSTIC_(self%id_rugc,rugc)

!SEAMLESS  call quota_flux( iiPel, ppphytoc ,ppO3c,ppphytoc, rugc, tfluxC )  
  _SET_ODE_(self%id_c,rugc)
  _SET_ODE_(self%id_O3c,-rugc)
!SEAMLESS  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR1c, 0.98D0 * rr1c, tfluxC ) !  flux is partitioned to non CDOM
  _SET_ODE_(self%id_c,-(1.00D0-self%p_fX1p) * rr1c)
  _SET_ODE_(self%id_R1c,(1.00D0-self%p_fX1p) * rr1c)
!SEAMLESS  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR1l, 0.02D0 * rr1c, tfluxC ) !  flux is partitioned to CDOM
  _SET_ODE_(self%id_c,-self%p_fX1p * rr1c)
  _SET_ODE_(self%id_X1c, self%p_fX1p * rr1c)

!CEA Activity excretion produces CDOM, nutrient-stress excretion dont  
!SEAMLESS  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR2l, 0.02D0 * flPIR2c, tfluxC ) ! flux to CDOM
  f2cdom = self%p_fX2p * ( (phytol/phytoc)/self%p_qlcPPY ) 
  _SET_ODE_(self%id_c, -f2cdom * flPIR2c)
  _SET_ODE_(self%id_X2c,f2cdom * flPIR2c)

  _SET_DIAGNOSTIC_(self%id_flPIR2c_act, flPIR2c)
  _SET_DIAGNOSTIC_(self%id_f2cdom, f2cdom)
  
!SEAMLESS  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR6c, rr6c, tfluxC )
  _SET_ODE_(self%id_c,-rr6c)
  _SET_ODE_(self%id_R6c,rr6c)
!SEAMLESS
!SEAMLESS  call quota_flux( iiPel, ppphytoc, ppphytoc,ppO3c, rrc, tfluxC )
  _SET_ODE_(self%id_c,-rrc)
  _SET_ODE_(self%id_O3c,rrc)
!SEAMLESS  call flux_vector( iiPel, ppO2o,ppO2o,-( rrc/ MW_C) )
  _SET_ODE_(self%id_O2o,-(rrc/MW_C))
!SEAMLESS  call flux_vector( iiPel, ppO2o,ppO2o, rugc/ MW_C ) 
  _SET_ODE_(self%id_O2o,rugc/MW_C)
!SEAMLESS
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Potential-Net prim prod. (mgC /m3/d)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 if (self%p_netgrowth) then
     sadap  =   max(  0.05_rk,  sum- slc)
 else
     sadap  =   et*self%p_sum
 end if

run  =   max(  ZERO, ( sum- slc)* phytoc)  ! net production

 _SET_DIAGNOSTIC_(self%id_sadap, sadap)
 _SET_DIAGNOSTIC_(self%id_run,run)

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Nutrient Uptake: calculate maximum uptake of N, P
 ! based on affinity
 ! Ammonium preference is considered if p_lN4 /= 0
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cqun3  =   self%p_lN4/( p_small + self%p_lN4 + N4n)
  rumn3  =   self%p_qun * N3n * phytoc * cqun3  ! max pot. uptake of N3
  rumn4  =   self%p_qun * N4n * phytoc          ! max pot. uptake of N4
  rumn  =   rumn3 + rumn4                       ! max pot. uptake of DIN
  rump  =   self%p_qup * N1p * phytoc           ! max pot. uptake of PO4

 _SET_DIAGNOSTIC_(self%id_cqun3, cqun3)
 _SET_DIAGNOSTIC_(self%id_rumn3, rumn3)
 _SET_DIAGNOSTIC_(self%id_rumn4, rumn4)
 _SET_DIAGNOSTIC_(self%id_rumn, rumn)
 _SET_DIAGNOSTIC_(self%id_rump, rump)

  if (self%p_netgrowth) then
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Check which fraction of fixed C can be used for new biomass
   ! given the internal storage.
   ! N and P uptake are compared sequentially
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth = min( run, ( rumn+ max( ZERO, 0.05_rk * &
      rugc * ( qncPPY - self%p_qnlc)))/ self%p_qnlc)
      netgrowth = min( netgrowth, ( rump+ max( ZERO, &
       0.05_rk * rugc*( qpcPPY - self%p_qplc)))/ self%p_qplc)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Excrete C that cannot be used for growth as carbo-hydrates:
   ! Correct the net C uptake
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth  =   max(  netgrowth,  ZERO)
      flPIR2c_tot  =   flPIR2c+ run- netgrowth
      flPIR2c  =  ((1.00D0-f2cdom)*flPIR2c) + run - netgrowth
      run  =   netgrowth
  end if

 _SET_DIAGNOSTIC_(self%id_netgrowth, netgrowth)

!!SEAMLESS  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR2c, 0.98D0 * flPIR2c, tfluxC ) ! flux to non CDOM
!  _SET_ODE_(self%id_c, -(1.00D0-f2cdom) * flPIR2c_tot)
!  _SET_ODE_(self%id_R2c,(1.00D0-f2cdom) * flPIR2c_tot)
!!SEAMLESS  call quota_flux( iiPel, ppphytoc, ppphytoc,ppR2l, 0.02D0 * flPIR2c, tfluxC ) ! flux to CDOM
!  _SET_ODE_(self%id_c, -f2cdom * flPIR2c_tot)
!  _SET_ODE_(self%id_X2c,f2cdom * flPIR2c_tot)

! CEA 98% of activity excretion + nutrient estress excretion produce only R2c
  _SET_ODE_(self%id_c,-flPIR2c)
  _SET_ODE_(self%id_R2c,flPIR2c)

 _SET_DIAGNOSTIC_(self%id_flPIR2c, flPIR2c)
 _SET_DIAGNOSTIC_(self%id_flPIR2c_tot, flPIR2c_tot)
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Specific net growth rate (d-1)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sunPPY  =   run/( p_small+ phytoc)

 _SET_DIAGNOSTIC_(self%id_sunPPY, sunPPY)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misn  =   sadap*( self%p_xqn * self%p_qncPPY * phytoc - phyton)  ! Intracellular missing amount of N
  rupn  =   self%p_xqn * self%p_qncPPY* run  ! N uptake based on net assimilat. C
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of NI

  r      =   insw(runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of Nn
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of Nn
!SEAMLESS  call quota_flux( iiPel, ppphyton, ppN3n,ppphyton, runn3, tfluxN )  ! source/sink.n
  _SET_ODE_(self%id_n,runn3)
  _SET_ODE_(self%id_N3n,-runn3)
!SEAMLESS  call quota_flux( iiPel, ppphyton, ppN4n,ppphyton, runn4, tfluxN )  ! source/sink.n
  _SET_ODE_(self%id_n,runn4)
  _SET_ODE_(self%id_N4n,-runn4)
!SEAMLESS GC: Alkalinity contributions: +1 for NO3 and  -1 for NH4: uptake of 1
!mole of NO3 increases alkalinity; uptake of 1 mole of NH4 decreases alkalinity (Wolf-Gladrow etal., 2007)
  _SET_ODE_(self%id_O3h,runn3-runn4)
  tmp = - runn*( ONE- r)
!SEAMLESS  call quota_flux( iiPel, ppphyton, ppphyton,ppR1n,tmp, tfluxN)  ! source/sink.n
  _SET_ODE_(self%id_n,-tmp)
  _SET_ODE_(self%id_R1n,tmp)

 _SET_DIAGNOSTIC_(self%id_misn, misn)
 _SET_DIAGNOSTIC_(self%id_rupn, rupn)
 _SET_DIAGNOSTIC_(self%id_runn, runn)
 _SET_DIAGNOSTIC_(self%id_runn3, runn3)
 _SET_DIAGNOSTIC_(self%id_runn4, runn4)
 _SET_DIAGNOSTIC_(self%id_fR1n, 0.0_rk)
 _SET_DIAGNOSTIC_(self%id_Nutil_O3h, runn3-runn4)

  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  misp  =   sadap*( self%p_xqp * self%p_qpcPPY * phytoc- phytop)  ! intracellular missing amount of P
  rupp  =   self%p_xqp * run* self%p_qpcPPY  ! P uptake based on C uptake
#ifdef EXTRACOST
  rupp  =   self%p_xqp * run* self%p_qpcPPY - ( sdo+ srs)* phytop  ! P uptake based on C uptake
#endif
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

 _SET_DIAGNOSTIC_(self%id_misp, misp)
 _SET_DIAGNOSTIC_(self%id_rupp, rupp)

  r  =   insw(runp)
  tmp = runp*r
 _SET_DIAGNOSTIC_(self%id_runp, tmp)
!SEAMLESS  call quota_flux(iiPel, ppphytop, ppN1p,ppphytop, tmp, tfluxP)  ! source/sink.p
  _SET_ODE_(self%id_p,tmp)
  _SET_ODE_(self%id_N1p,-tmp)
!SEAMLESS GC: Alkalinity contributions: +1 for PO4 (i.e., uptake of 1 mole of PO4 
! increases alkalinity by 1 mole  (Wolf-Gladrow etal., 2007)
!  _SET_ODE_(self%id_O3h,tmp)
  _SET_DIAGNOSTIC_(self%id_Putil_O3h,tmp)

  tmp = - runp*( ONE- r)
!SEAMLESS  call quota_flux(iiPel, ppphytop, ppphytop,ppR1p, tmp, tfluxP)  ! source/sink.p
  _SET_ODE_(self%id_p,-tmp)
  _SET_ODE_(self%id_R1p,tmp)

 _SET_DIAGNOSTIC_(self%id_fR1p, tmp)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr6n =     pe_R6     * sdo * phyton
  rr1n = (ONE - pe_R6) * sdo * phyton

  rr6p =     pe_R6     * sdo * phytop
  rr1p = (ONE - pe_R6) * sdo * phytop

 _SET_DIAGNOSTIC_(self%id_rr6n, rr6n)
 _SET_DIAGNOSTIC_(self%id_rr1n, rr1n)

 _SET_DIAGNOSTIC_(self%id_rr6p, rr6p)
 _SET_DIAGNOSTIC_(self%id_rr1p, rr1p)
!SEAMLESS  call quota_flux( iiPel, ppphyton, ppphyton,ppR6n, rr6n, tfluxN )  ! source/sink.n
  _SET_ODE_(self%id_n,-rr6n)
  _SET_ODE_(self%id_R6n,rr6n)
!SEAMLESS  call quota_flux( iiPel, ppphyton, ppphyton,ppR1n, rr1n, tfluxN )  ! source/sink.n
  _SET_ODE_(self%id_n,-rr1n)
  _SET_ODE_(self%id_R1n,rr1n)
!SEAMLESS  call quota_flux( iiPel, ppphytop, ppphytop,ppR6p, rr6p, tfluxP )  ! source/sink.p
  _SET_ODE_(self%id_p,-rr6p)
  _SET_ODE_(self%id_R6p,rr6p)
!SEAMLESS  call quota_flux( iiPel, ppphytop, ppphytop,ppR1p, rr1p, tfluxP )  ! source/sink.p
  _SET_ODE_(self%id_p,-rr1p)
  _SET_ODE_(self%id_R6p,rr1p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( self%use_Si )  then
    select case (self%p_switchSi)
    case (1)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Gross uptake of silicate excluding respiratory costs
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      runs = max(ZERO, self%p_qscPPY * (sum-srs) * phytoc)
    case (2)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Silicate uptake based on intracellular needs (note, no luxury)
      !  There can be efflux of dissolved silicate (M-J et al., 2000)
      !  however this generates fake remineralization and it is not implemented
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rums  =   self%p_qus * N5s * phytoc  ! max pot uptake based on affinity
      miss  =   sadap*max(ZERO, self%p_qscPPY * phytoc - phytos) ! intracellular missing Si
      rups  =   run* self%p_qscPPY  ! Si uptake based on net C uptake
      runs  =   min(  rums,  rups+ miss)  ! actual uptake
    end select
              
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Uptake and Losses of Si (only lysis)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS    call flux_vector( iiPel, ppN5s,ppphytos, runs)
  _SET_ODE_(self%id_N5s,-runs)
  _SET_ODE_(self%id_s,runs)
!SEAMLESS    call flux_vector( iiPel, ppphytos, ppR6s, sdo*phytos )
  _SET_ODE_(self%id_s,-sdo*phytos)
  _SET_ODE_(self%id_R6s,sdo*phytos)

 _SET_DIAGNOSTIC_(self%id_rums, rums)
 _SET_DIAGNOSTIC_(self%id_miss, miss)
 _SET_DIAGNOSTIC_(self%id_rups, rups)
 _SET_DIAGNOSTIC_(self%id_runs, runs)
  endif
!SEAMLESS
!SEAMLESS#ifdef INCLUDE_PELFE
!SEAMLESS  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS  ! Nutrient dynamics: IRON
!SEAMLESS  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS
!SEAMLESS  if (ppphytof > 0) then
!SEAMLESS     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS     ! Net uptake
!SEAMLESS     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS     rumf  =   p_quf(phyto)* N7f(:)* phytoc  ! max potential uptake
!SEAMLESS     ! intracellular missing amount of Fe
!SEAMLESS     misf  =   sadap*max(ZERO,p_xqf(phyto)*p_qfcPPY(phyto)*phytoc - phytof)  
!SEAMLESS     rupf  =   p_xqf(phyto)* run* p_qfcPPY(phyto)  ! Fe uptake based on C uptake
!SEAMLESS     runf  =   min(  rumf,  rupf+ misf)  ! actual uptake
!SEAMLESS     r  =   insw(runf)
!SEAMLESS     ! uptake from inorganic if shortage
!SEAMLESS     call flux_vector( iiPel, ppN7f,ppphytof, runf* r )
!SEAMLESS     ! release to dissolved organic to keep the balance if excess
!SEAMLESS     call flux_vector(iiPel, ppphytof,ppR1f,- runf*( ONE- r))
!SEAMLESS   
!SEAMLESS     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS     ! Losses of Fe
!SEAMLESS     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS     rr6f  =   rr6c* p_qflc(phyto)
!SEAMLESS     rr1f  =   sdo* phytof- rr6f
!SEAMLESS     call flux_vector( iiPel, ppphytof,ppR1f, rr1f )
!SEAMLESS     call flux_vector( iiPel, ppphytof,ppR6f, rr6f )
!SEAMLESS  end if
!SEAMLESS#endif
!SEAMLESS
!SEAMLESS if ( self%ChlDynamicsFlag== 2) then
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Chl-a synthesis and photoacclimation
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   select case (self%p_switchChl)
     case (1) ! PELAGOS
          rho_Chl = self%p_qlcPPY* min(ONE, self%p_sum * eiPPY* phytoc/( &
                    photochem*( phytol+ p_small)* parEIR))
          rate_Chl = rho_Chl*(sum - seo - sea - sra) * phytoc - sdo*phytol
     case (2) ! OPATM-BFM
          rho_Chl  = self%p_qlcPPY * sum/( photochem* qlcPPY * parEIR * et + p_small)
          rate_Chl = iN* rho_Chl* run- max( self%p_sdchl * ( ONE - iN), sdo)* &
              phytol+ min( ZERO, sum- slc+ sdo)* max( ZERO, phytol- self%p_qlcPPY * phytoc)
     case (3) ! UNIBO
          rho_Chl = self%p_qlcPPY*min(ONE,          &
                    (sum-seo-sea-sra) *phytoc /          &
                    (photochem*(phytol+p_small) *parEIR))
          ! The "optimal" chl concentration corresponds to the chl that
          ! (given the actual C biomass) would give (Epar/Ek)=p_EpEk
          chl_opt = self%p_EpEk_or * self%p_sum*phytoc/  &
                    (photochem*parEIR+p_small)
          !  Actual chlorophyll concentration exceeding the "optimal" value is 
          !  discarded with a p_tochl_relt relaxation.
          rate_Chl = rho_Chl*(sum-seo-sea-sra)*phytoc-(sdo+srs)*phytol - &
                     max(ZERO,(phytol-chl_opt))*self%p_tochl_relt
     case (4) ! NIOZ
         ! total synthesis, only when there is net production (run > 0)
         ! The fixed loss rate due to basal respiration is introduced to have 
         ! chl loss in the absence of light (< 1 uE/m2/s)
          rho_Chl = self%p_qlcPPY * min(ONE, self%p_sum * eiPPY* phytoc/( &
                    photochem * ( phytol+ p_small)* parEIR))
             rate_Chl = rho_Chl*run - self%p_sdchl * phytol*max( ZERO, ( self%p_thdo-tN)) &
                    -srs * phytol * ONE/(parEIR+ONE)
   end select
!SEAMLESS    call flux_vector( iiPel, ppphytol,ppphytol, rate_Chl )
  _SET_ODE_(self%id_chl,rate_Chl)
!SEAMLESS end if

 _SET_DIAGNOSTIC_(self%id_rho_Chl, rho_Chl)
 _SET_DIAGNOSTIC_(self%id_rate_Chl, rate_Chl)
!SEAMLESS#if defined INCLUDE_PELCO2
!SEAMLESS  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS  ! PIC (calcite/aragonite) production
!SEAMLESS  ! The idea in PISCES is that the calcite flux exists only when associated
!SEAMLESS  ! to a carbon release from phytoplankton (there is no calcite storage in phyto)
!SEAMLESS  ! First compute the realized rain ratio for each phytoplankton species
!SEAMLESS  ! The presence of PIC in phytoplankton group is controlled by p_caco3r
!SEAMLESS  ! with the following regulating factors:
!SEAMLESS  !  - nutrient limitation
!SEAMLESS  !  - temperature enhancement
!SEAMLESS  !  - density enhancement
!SEAMLESS  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( self%use_CaCO3 )  then   ! used only in P2 -> it is very crapy solution
   if ( self%p_caco3r > ZERO ) then
     qccPPY = min(0.8_rk,self%p_caco3r*tN*et*MM(phytoc, self%p_sheo))
     qccPPY = max(0.02_rk,qccPPY)
! Calcite production represented as a flux between DIC and PIC, impacting ALK
!SEAMLESS     call flux_vector( iiPel, ppO3c,ppO5c, qccPPY(phyto, :)*rr6c )
     _SET_ODE_(self%id_O3c,-qccPPY*rr6c) 
     _SET_ODE_(self%id_O5c,+qccPPY*rr6c) 
!SEAMLESS     call flux_vector( iiPel, ppO3h,ppO3h, -C2ALK*qccPPY(phyto, :)*rr6c )
     _SET_ODE_(self%id_O3h, -2_rk*qccPPY*rr6c)
     _SET_DIAGNOSTIC_(self%id_O3hconume_for_CaCO3prec, -2_rk*qccPPY*rr6c)
  endif
  endif

!SEAMLESS#endif

 ! End of computation section for process PhytoDynamics

!SEAMLESS  end subroutine PhytoDynamics

      ! Leave spatial loops (if any)
      _SET_DIAGNOSTIC_(self%id_flPIR2c,flPIR2c)

      _LOOP_END_

   end subroutine do
   function get_sinking_rate(self,_ARGUMENTS_LOCAL_) result(sediPPY)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! Sedimentation
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS  if ( p_res(phyto)> ZERO) then
!SEAMLESS    sediPPY(phyto,:) = sediPPY(phyto,:) &
!SEAMLESS                   + p_res(phyto)* max( ZERO, ( p_esNI(phyto)-tN))
!SEAMLESS  end if
!SEAMLESS
   class (type_ogs_bfm_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_LOCAL_

      real(rk) :: phytoc,phytop,phyton,phytos
      real(rk) :: qpcPPY,qncPPY,qscPPY
      real(rk) :: N5s
      real(rk) :: sediPPY
      real(rk) :: iN1p, iNIn, eN5s, fpplim, iN5s, iN, tN  

       _GET_(self%id_c,phytoc)
       _GET_(self%id_p,phytop)
       _GET_(self%id_n,phyton)

       if (self%use_Si) then
           _GET_(self%id_s,phytos)
       endif

  ! Quota collectors
       qpcPPY = phytop/(phytoc+p_small) ! add some epsilon (add in shared) to avoid divide by 0
       qncPPY = phyton/(phytoc+p_small)
       if (self%use_Si) then
           qscPPY=phytos/(phytoc+p_small)
       endif

       iN1p = min( ONE, max( p_small, ( qpcPPY - self%p_qplc)/( self%p_qpcPPY- self%p_qplc)))
       iNIn = min( ONE, max( p_small, ( qncPPY - self%p_qnlc)/( self%p_qncPPY- self%p_qnlc)))

       if (self%use_Si) then
           _GET_(self%id_N5s,N5s)
           select case (self%p_switchSi) 
               case (1)  ! external control
                   eN5s = min( ONE, N5s/(N5s + self%p_chPs+(self%p_Contois*phytos)))
                   fpplim = eN5s
                   iN5s   =  ONE
               case (2) ! internal control
                   iN5s = min(ONE, max( p_small, ( qscPPY &
                        - self%p_qslc)/( self%p_qscPPY- self%p_qslc)))
                   fpplim = iN5s
                   eN5s   = ONE
           end select
       else 
           iN5s   = ONE
           eN5s   = ONE
           fpplim = ONE
       end if


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Multiple nutrient limitation
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       select case ( self%p_limnut)
           case ( 0 )
               iN  =   (iN1p* iNIn)**(0.5_rk)  ! geometric mean

           case ( 1 )
               iN  =   min(  iN1p,  iNIn)  ! Liebig rule

           case ( 2 )
               iN  =   2.0_rk/( ONE/ iN1p+ ONE/ iNIn)  ! combined

       end select

! tN only controls sedimentation of phytoplankton (Liebig)
       tN= min(iN,fpplim)

      sediPPY = self%rm + self%p_res * MAX(ZERO,( self%p_esNI-tN))
!SEAMLESS  if ( p_res(phyto)> ZERO) then
!SEAMLESS    sediPPY(phyto,:) = sediPPY(phyto,:) &
!SEAMLESS                   + p_res(phyto)* max( ZERO, ( p_esNI(phyto)-tN))
!SEAMLESS  end if
   end function

   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_ogs_bfm_primary_producer),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: sediPPY

      _LOOP_BEGIN_

         ! Retrieve local state
         sediPPY = -get_sinking_rate(self,_ARGUMENTS_LOCAL_)
         _SET_VERTICAL_MOVEMENT_(self%id_c,sediPPY)
         _SET_VERTICAL_MOVEMENT_(self%id_p,sediPPY)
         _SET_VERTICAL_MOVEMENT_(self%id_n,sediPPY)
         if (self%use_Si) _SET_VERTICAL_MOVEMENT_(self%id_s,sediPPY)
         _SET_VERTICAL_MOVEMENT_(self%id_chl,sediPPY)
         if (use_iron) _SET_VERTICAL_MOVEMENT_(self%id_f,sediPPY)

      _LOOP_END_

   end subroutine get_vertical_movement
end module
!SEAMLESS!EOC
!SEAMLESS!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!SEAMLESS! MODEL  BFM - Biogeochemical Flux Model 
!SEAMLESS!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
