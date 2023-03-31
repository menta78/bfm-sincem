#include "fabm_driver.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MesoZoo
!
! DESCRIPTION
!   This submodel describes the carbon dynamics and associated
!   nutrient dynamics in mesozooplankton 
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
! !INTERFACE
 module bfm_Mesozoo

   use fabm_types
   use fabm_particle

   use ogs_bfm_shared
   use ogs_bfm_pelagic_base

!
! !USES:

!  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  ! Modules (use of ONLY is strongly encouraged!)
!  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!  use global_mem, ONLY:RLEN, ONE, ZERO
!#ifdef NOPOINTERS
!  use mem
!#else
!  use mem, ONLY: D3STATE, O2o, N1p, N4n, R6c, R6p, R2c, &
!    R6n, PhytoPlankton, MicroZooPlankton, MesoZooPlankton
!  use mem, ONLY: Depth, ppO2o, ppMMHg, ppN1p, ppN4n, ppR6c, ppR6n, ppR6p, ppR6s, &
!    ppPhytoPlankton, ppMicroZooPlankton, ppMesoZooPlankton, ETW, &
!    qncPPY, qpcPPY, qlcPPY, qscPPY, qncMIZ, qpcMIZ, qncMEZ, qpcMEZ, iiPhytoPlankton, &
!    iiMicroZooPlankton, iiMesoZooPlankton, iiC, iiN, iiP, iiL, iiS, NO_BOXES, &
!    iiBen, iiPel, flux_vector, quota_flux
!#ifdef INCLUDE_PELCO2
!  use mem, ONLY: ppO3c, ppO5c, ppO3h, qccPPY
!#endif
!#ifdef INCLUDE_PELFE
!  use mem, ONLY: iiF, qfcPPY, ppR6f
!#endif
!use mem, ONLY: iiH, qhcPPY, qhcMIZ,qhcMEZ
!#endif
!#ifdef BFM_GOTM
!  use mem, ONLY: jnetMeZc
!#endif
!  use mem_Param,  ONLY: p_small
!  use bfm_error_msg, ONLY: bfm_error
!  use constants,ONLY: MIN_VAL_EXPFUN, MW_C, C2ALK
!  use mem_MesoZoo

!  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  ! The following vector functions are used: eTq, MM, MM_power, nutlim
!  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  use mem_globalfun,   ONLY: eTq, MM, MM_power, nutlim

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  implicit none

  private

! !INPUT:
!  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  integer,intent(IN)  :: zoo

!  
!
! !AUTHORS
!   First ERSEM version by N. Broekhuizen and A.D. Bryant
!   Additional parametrizations by P. Ruardij and M. Vichi 
!   Dynamical allocation by G. Mattia 
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi
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
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!

   type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_mesozoo
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving
      ! from type_ogs_bfm_pelagic_base!

      !! Identifiers for state variables of other models

      !! In ERSEM
      !type (type_dependency_id),    allocatable,dimension(:) :: id_preyc,id_preyn,id_preyp,id_preys,id_preyf,id_preyl

      type (type_state_variable_id), allocatable,dimension(:) :: id_preyc     !  carbon prey
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyn
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyp
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyl
      type (type_state_variable_id), allocatable,dimension(:) :: id_preys
      ! type (type_state_variable_id),    allocatable,dimension(:) :: id_preyf
      type (type_model_id),      allocatable,dimension(:) :: id_prey

      type (type_state_variable_id)      :: id_O2o, id_O3c                  !  disolved oxygen, dissolved inorganic carbon
      type (type_state_variable_id)      :: id_N1p,id_N4n                   !  nutrients: phosphate, ammonium
      type (type_state_variable_id)      :: id_R6c,id_R6p,id_R6n            !  particulate organic carbon
      type (type_state_variable_id)      :: id_R6s                          !  biogenic silica
      type (type_state_variable_id)      :: id_O3h                          !  alkalinity
      type (type_state_variable_id)      :: id_O5c                          !  calcite
  
      !! Environmental dependencies
      type (type_dependency_id)          :: id_ETW                          ! temperature

      !! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_qncMEZ     ! N:C quontum
      type (type_diagnostic_variable_id) :: id_qpcMEZ     ! P:C quontum      
      type (type_diagnostic_variable_id) :: id_prey2l     ! test      
      type (type_diagnostic_variable_id) :: id_prey1x     ! test      

      type (type_diagnostic_variable_id) :: id_ETWd   ! temperature Celsius
      type (type_diagnostic_variable_id) :: id_et     ! physiological temperature response
      type (type_diagnostic_variable_id) :: id_eo     ! oxygen temperature response

      type (type_diagnostic_variable_id) :: id_rumc   ! total potential food
      type (type_diagnostic_variable_id) :: id_rugc   ! total food uptake rate (eq 38 Vichi et al. 2007)
      type (type_diagnostic_variable_id) :: id_sut    ! specific uptake rate considering potentially available food

      type (type_diagnostic_variable_id) :: id_rut_c  ! carbon ingestion rate
      type (type_diagnostic_variable_id) :: id_rut_n  ! nitrogen ingestion rate
      type (type_diagnostic_variable_id) :: id_rut_p  ! phosphorus ingestion rate
      type (type_diagnostic_variable_id) :: id_rrc    ! respiration rate

      type (type_diagnostic_variable_id) :: id_rdo_c  ! low oxygen mortality rate
      type (type_diagnostic_variable_id) :: id_rd_c   ! density dependent mortality rate
      type (type_diagnostic_variable_id) :: id_rq6c   ! carbon egestion rate
      type (type_diagnostic_variable_id) :: id_rq6n   ! nitrogen egestion rate
      type (type_diagnostic_variable_id) :: id_rq6p   ! phosphorus egestion rate

      type (type_diagnostic_variable_id) :: id_ren    ! ammonium remineralization rate
      type (type_diagnostic_variable_id) :: id_rep    ! phosphate remineralization rate

      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_ruPPYc    ! prey-specific grazing
      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preyld    ! prey chl for diagnostic
      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_CaCO3precip ! precipitation of PIC
      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_CaCO3_to_O3h ! consume of alk due to precipitation of PIC

      type (type_diagnostic_variable_id) :: id_temp_p    ! 
      type (type_diagnostic_variable_id) :: id_temp_n    ! 
      type (type_diagnostic_variable_id) :: id_limit     ! which constituent is limiting
      type (type_diagnostic_variable_id) :: id_pe_R6c    ! rate removal C
      type (type_diagnostic_variable_id) :: id_pe_N1p    ! rate removal P
      type (type_diagnostic_variable_id) :: id_pe_N4n    ! rate removal N
      type (type_diagnostic_variable_id) :: id_varO3h_Nutil ! variation of O3h due to NH44 utilization by Zoo
      type (type_diagnostic_variable_id) :: id_varO3h_Putil ! variation of O3h due to PO4 utilization by Zoo
 
      !! Parameters (described in subroutine initialize, below)
      integer  :: nprey
      real(rk), allocatable :: p_pa(:)
      integer, allocatable :: p_isP2(:)
!      logical, allocatable :: p_pl(:)
!      logical, allocatable :: p_ps(:)
      real(rk) :: p_q10, p_srs, p_sum, p_sd
      real(rk) :: p_vum, p_puI, p_peI, p_sdo, p_sds
      real(rk) :: p_pecaco3, p_qpcMEZ, p_qncMEZ, p_clO2o
!      Examples
!      real(rk) :: p_paPPY, p_paMIZ, p_paMEZ   ! diet matrix
!      integer :: p_switchDOC, p_switchSi,p_limnut,p_switchChl
!      logical :: use_Si,p_netgrowth

   contains

      ! Model procedures
      procedure :: initialize
      procedure :: do

   end type type_ogs_bfm_mesozoo

!   ! Constants
!   real(rk),parameter :: pippo = 0.002_rk

contains

  subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_mesozoo),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer           :: iprey
      character(len=16) :: index
      real(rk) :: pippo1
!      logical  :: preyisphyto, preyisdiat

!EOP
!-------------------------------------------------------------------------!
!BOC
!
! Obtain the values of all model parameters from FABM.
      ! Specify the long name and units of the parameters, which could be used
      ! by FABM (or its host) to present parameters
      ! to the user for configuration (e.g., through a GUI)
      call self%get_parameter(self%p_q10,    'p_q10',     '-',          'Q10 value for physiological rates')
      call self%get_parameter(self%p_srs,    'p_srs',     '1/d',        'respiration rate at reference temperature')
      call self%get_parameter(self%p_sum,    'p_sum',     '1/d',        'maximal productivity at reference temperature')
      call self%get_parameter(self%p_sd,     'p_sd',      '1/d',        'backgroung natural mortality')
      call self%get_parameter(self%p_puI,    'p_puI',     '-',          'assimilation efficiency')
      call self%get_parameter(self%p_peI,    'p_peI',     '-',          'fraction of phaeces production')
      call self%get_parameter(self%p_vum,    'p_vum',     'm3/mgC/d',   'specific search volume')
      call self%get_parameter(self%p_sdo,    'p_sdo',     'm3/mgC/d',   'specific density-dependent mortality')
      call self%get_parameter(self%p_sds,    'p_sds',     '-',          'exponent of density-dependent mortality')
      call self%get_parameter(self%p_pecaco3,'p_pecaco3', '-',          'portion of egested calcified shells during grazing')
      call self%get_parameter(self%p_qpcMEZ, 'p_qpcMEZ',  'mmolP/mgC',  'maximum quotum P:C')
      call self%get_parameter(self%p_qncMEZ, 'p_qncMEZ',  'mmolN/mgC',  'maximum quotum N:C')
      call self%get_parameter(self%p_clO2o,  'p_clO2o',   'mmolO2/m3',  'half-saturation oxygen concentration')

! Register state variables (handled by type_bfm_pelagic_base)
!     call self%initialize_ogs_bfm_base(sedimentation=.true.)
      call self%initialize_bfm_base()
      call self%add_constituent('c',1.e-4_rk)
      call self%add_constituent('n',1.26e-6_rk)
      call self%add_constituent('p',4.288e-8_rk)
!      call self%add_constituent('f',5.e-6_rk)  ! NB this does nothing if iron support is disabled.

! Register links to external preys
      call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
      ! Get prey-specific parameters.
      allocate(self%p_pa(self%nprey))     !Availability of nprey for predator
      allocate(self%p_isP2(self%nprey))   !is P2? [=1 for P2 and 0 otherwise]
!      allocate(self%p_pl(self%nprey))     !Does the prey have Chl?
!      allocate(self%p_ps(self%nprey))     !Does the prey have Silica?
      allocate(self%id_prey(self%nprey))
      allocate(self%id_preyc(self%nprey))
      allocate(self%id_preyn(self%nprey))
      allocate(self%id_preyp(self%nprey))
      allocate(self%id_preyl(self%nprey))
      allocate(self%id_preys(self%nprey))

      allocate(self%id_preyld(self%nprey))
      allocate(self%id_ruPPYc(self%nprey))

      allocate(self%id_CaCO3precip(self%nprey))
      allocate(self%id_CaCO3_to_O3h(self%nprey))

      do iprey=1,self%nprey
        write (index,'(i0)') iprey
        call self%get_parameter(self%p_pa(iprey),'suprey'//trim(index),'-','Availability for prey type '//trim(index))
        call self%get_parameter(self%p_isP2(iprey),'isP2'//trim(index),'-','identify P2 among the preys '//trim(index))
        
        call self%register_state_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mg C/m^3',   'prey '//trim(index)//' carbon')
        call self%register_state_dependency(self%id_preyn(iprey),'prey'//trim(index)//'n','mmol N/m^3', 'prey '//trim(index)//' nitrogen')
        call self%register_state_dependency(self%id_preyp(iprey),'prey'//trim(index)//'p','mmol P/m^3', 'prey '//trim(index)//' phosphorous')

        call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
        call self%request_coupling_to_model(self%id_preyc(iprey),self%id_prey(iprey),'c')
        call self%request_coupling_to_model(self%id_preyn(iprey),self%id_prey(iprey),'n')
        call self%request_coupling_to_model(self%id_preyp(iprey),self%id_prey(iprey),'p')


        call self%register_state_dependency(self%id_preyl(iprey),'prey'//trim(index)//'Chl','mg Chl/m^3', 'prey '//trim(index)//' chlorophyll')
        call self%request_coupling_to_model(self%id_preyl(iprey),self%id_prey(iprey),total_chlorophyll)


        call self%register_diagnostic_variable(self%id_preyld(iprey),'prey'//trim(index)//'Chld','mg Chl/m^3',  'prey '//trim(index)//' chlorophyll',output=output_none)
        call self%register_diagnostic_variable(self%id_ruPPYc(iprey),'prey'//trim(index)//'rate','mg C/mg C/d', 'prey '//trim(index)//' grazing rate',output=output_none)


        call self%register_state_dependency(self%id_preys(iprey),'prey'//trim(index)//'s','mmol Si/m^3', 'prey '//trim(index)//' silica')
        call self%request_coupling_to_model(self%id_preys(iprey), self%id_prey(iprey),standard_variables%total_silicate)

!#ifdef INCLUDE_PELFE
!        call self%register_state_dependency(self%id_preyc(iprey),'prey'//trim(index)//'f','umol Fe/m^3',   'prey '//trim(index)//' iron')
!        call self%request_coupling_to_model(self%id_preyf(iprey),self%id_prey(iprey),'f')
!#endif
     enddo



! Register environmental dependencies (temperature)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
 
! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O2/m^3','dissolved oxygen')
      call self%register_state_dependency(self%id_O3c,'O3c','mg C/m^3'   ,'dissolved inorganic carbon')
      call self%register_state_dependency(self%id_O3h,'O3h','mmol /m^3'  ,'alkalinity')
      call self%register_state_dependency(self%id_O5c,'O5c' ,'mgC/m^3'    ,'calcite'    ,required=.false.)
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3' ,'phosphate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3' ,'ammonium')
      call self%register_state_dependency(self%id_R6c,'R6c','mg C/m^3'   ,'POC')
      call self%register_state_dependency(self%id_R6p,'R6p','mmol P/m^3' ,'POP')
      call self%register_state_dependency(self%id_R6n,'R6n','mmol N/m^3' ,'PON')
      call self%register_state_dependency(self%id_R6s,'R6s','mmol Si/m^3','biogenic silica')

! Register diagnostic variables (i.e., model outputs)
      call self%register_diagnostic_variable(self%id_qncMEZ,'qncMEZ', 'mmolN/mgC', 'N:C quontum',output=output_none)
      call self%register_diagnostic_variable(self%id_qpcMEZ,'qpcMEZ', 'mmolP/mgC', 'P:C quontum',output=output_none)
      call self%register_diagnostic_variable(self%id_prey2l,'prey2l', 'Chl', 'test',output=output_none)
      call self%register_diagnostic_variable(self%id_prey1x,'prey1x', 's', 'test',output=output_none)
      call self%register_diagnostic_variable(self%id_ETWd,  'ETW',   'C',        'temperature Celsius',output=output_none)
      call self%register_diagnostic_variable(self%id_et,    'et',    '-',        'temperature factor',output=output_none)
      call self%register_diagnostic_variable(self%id_eo,    'eo',    '-',        'oxygen limitation',output=output_none)
      call self%register_diagnostic_variable(self%id_rumc,  'rumc',  'mgC/m3',   'total potential food',output=output_none)
      call self%register_diagnostic_variable(self%id_rugc,  'rugc',  'mgC/m3/d', 'total food uptake rate',output=output_none)
      call self%register_diagnostic_variable(self%id_sut,   'sut',   '1/d',      'specific uptake rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rut_c, 'rut_c', 'mgC/m3/d', 'carbon ingestion rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rut_n, 'rut_n', 'mgC/m3/d', 'nitrogen ingestion rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rut_p, 'rut_p', 'mgC/m3/d', 'phosphorus ingestion rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rrc,   'rrc',   'mgC/m3/d', 'respiration rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rdo_c, 'rdo_c', 'mgC/m3/d', 'low oxygen mortality rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rd_c,  'rd_c',  'mgC/m3/d', 'density dependent mortality rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rq6c,  'rq6c',  'mgC/m3/d', 'carbon egestion rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rq6n,  'rq6n',  'mgC/m3/d', 'nitrogen egestion rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rq6p,  'rq6p',  'mgC/m3/d', 'phosphorus egestion rate',output=output_none)
      call self%register_diagnostic_variable(self%id_ren,   'ren',   'mmolN/m3/d',  'ammonium remineralization rate',output=output_none)
      call self%register_diagnostic_variable(self%id_rep,   'rep',   'mmolP/m3/d',  'phosphate remineralization rate',output=output_none)
      call self%register_diagnostic_variable(self%id_varO3h_Nutil,'varO3h_Nutil','mmol/m3/d','variaz O3h due to N utiliz',output=output_none)
      call self%register_diagnostic_variable(self%id_varO3h_Putil,'varO3h_Putil','mmol/m3/d','variaz O3h due to P utiliz',output=output_none)

      call self%register_diagnostic_variable(self%id_temp_p,  'temp_p',   '-',  '-',output=output_none)
      call self%register_diagnostic_variable(self%id_temp_n,  'temp_n',   '-',  '-',output=output_none)
      call self%register_diagnostic_variable(self%id_limit,   'limit',    '-',  'limiting constituent',output=output_none)
      call self%register_diagnostic_variable(self%id_pe_R6c,  'pe_R6c',   'mgC/m3/d',    'removal of C',output=output_none)
      call self%register_diagnostic_variable(self%id_pe_N1p,  'pe_N1p',   'mmolP/m3/d',  'removal of P',output=output_none)
      call self%register_diagnostic_variable(self%id_pe_N4n,  'pe_N4n',   'mmolN/m3/d',  'removal of N',output=output_none)


     do iprey=1,self%nprey
       write (index,'(i0)') iprey
       if (self%p_isP2(iprey).eq.1) then
       call self%register_diagnostic_variable(self%id_CaCO3precip(iprey),'_'//trim(index)//'_CaCO3precip','mg C/m^3/d','prey '//trim(index)//' CaCO3precip',output=output_none)
       call self%register_diagnostic_variable(self%id_CaCO3_to_O3h(iprey),'_'//trim(index)//'_consumeO3h_for_CaCO3precip','mmol/m^3/d','prey '//trim(index)//' consumeO3h_for_CaCO3precip',output=output_none)
       endif
      end do

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_ogs_bfm_mesozoo),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_


!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  integer  :: i
!  integer  :: ppzooc, ppzoon, ppzoop,ppzooh
!  integer, save :: first =0
!  integer,dimension(NO_BOXES)  :: limit
!  real(RLEN),allocatable,save,dimension(:) :: sut,temp_p,temp_n,rumc,rugc,eo,  &
!                                       et,rrs_c,rrs_n,rrs_p,rut_c, &
!                                       rut_n,rut_p,rd_c,rd_n,rd_p,sdo,rdo_c,  &
!                                       rdo_n,rdo_p,ret_c,ret_n,ret_p,ru_c, &
!                                       ru_n,ru_p,pu_e_n,pu_e_p,prI,pe_R6c
!
!  real(RLEN),allocatable,save,dimension(:) :: pe_N1p,pe_N4n,ruPPYc,ruMIZc,ruMEZc,rq6c, &
!                                       rq6n,rq6p,rrc,ren,rep,tfluxC, tfluxN, tfluxP,tfluxH,     &
!                                       zooc,zoop,zoon
!  real(RLEN),allocatable,save,dimension(:,:) :: PPYc,MIZc,MEZc
!  real(RLEN),allocatable,save,dimension(:) :: net,r
!  integer :: AllocStatus
!#ifndef INCLUDE_PELCO2
!  integer,parameter :: ppO3c = 0
!#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    integer  :: iprey
    real(rk), dimension(self%nprey) :: preycP,preypP,preynP,preylP,preysP   !, preyfP
    real(rk), dimension(self%nprey) :: rupreyc, PPYc
    real(rk) :: preyP
    real(rk) :: zooc, zoop, zoon
    real(rk) :: ETW,et,eo
    real(rk) :: O2o
!    real(rk) :: O3c
!    real(rk) :: N1p,N4n
!    real(rk) :: R6c,R6p,R6n
    real(rk) :: rumc,rugc,sut
    real(rk) :: rut_c,rut_n,rut_p
    real(rk) :: qpcMEZ,qncMEZ 
    real(rk) :: ruPPYc
    real(rk) :: prI, rrc, rdo_c, rd_c
    real(rk) :: rq6c, rq6n, rq6p
    real(rk) :: rep, ren
    real(rk) :: temp_n, temp_p
    real(rk) :: ru_c, ru_p, ru_n
    real(rk) :: pu_e_n, pu_e_p
    integer       :: limit
    real(rk) :: pe_R6c, pe_N1p, pe_N4n
    

    ! Enter spatial loops (if any)
    _LOOP_BEGIN_



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  if (first==0) then
!     ALLOCATE ( PPYc(NO_BOXES,iiPhytoPlankton), MIZc(NO_BOXES,iiMicroZooPlankton), &
!        &       MEZc(NO_BOXES,iiMesoZooPlankton),  &
!        &       zooc(NO_BOXES), zoop(NO_BOXES), zoon(NO_BOXES),        &
!        &       rep(NO_BOXES), ren(NO_BOXES), rrc(NO_BOXES),           &
!        &       rq6p(NO_BOXES), rq6n(NO_BOXES), rq6c(NO_BOXES),        &
!        &       ruPPYc(NO_BOXES), ruMIZc(NO_BOXES), ruMEZc(NO_BOXES),  &
!        &       pe_N4n(NO_BOXES), pe_N1p(NO_BOXES), pe_R6c(NO_BOXES),  &
!        &       prI(NO_BOXES), pu_e_p(NO_BOXES), pu_e_n(NO_BOXES),     &
!        &       ru_p(NO_BOXES), ru_n(NO_BOXES), ru_c(NO_BOXES),        &
!        &       ret_p(NO_BOXES), ret_n(NO_BOXES), ret_c(NO_BOXES),     &
!        &       rdo_p(NO_BOXES), rdo_n(NO_BOXES), rdo_c(NO_BOXES),     &
!        &       rd_p(NO_BOXES) , rd_n(NO_BOXES) , rd_c(NO_BOXES) ,     &
!        &       rut_p(NO_BOXES), rut_n(NO_BOXES), rut_c(NO_BOXES),     &
!        &       eo(NO_BOXES), sdo(NO_BOXES), et(NO_BOXES), sut(NO_BOXES),   &
!        &       rumc(NO_BOXES), rugc(NO_BOXES), net(NO_BOXES), r(NO_BOXES), &
!        &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES),  &
!        &       temp_p(NO_BOXES), temp_n(NO_BOXES),tfluxH(NO_BOXES),   &
!        &      STAT = AllocStatus )

!     IF( AllocStatus /= 0 ) call bfm_error('MesoZooDynamics','Error allocating arrays')
!     first=1
!  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Copy state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

         ! Concentrations excluding background (used in sink terms)
!  ppzooc = ppMesoZooPlankton(zoo,iiC)
!  ppzoon = ppMesoZooPlankton(zoo,iiN)
!  ppzoop = ppMesoZooPlankton(zoo,iiP)

         _GET_(self%id_c,zooc)
         _GET_(self%id_p,zoop)
         _GET_(self%id_n,zoon)

!  zooc = D3STATE(ppzooc,:)
!  zoon = zooc * qncMEZ(zoo,:)
!  zoop = zooc * qpcMEZ(zoo,:)

         ! Retrieve environmental dependencies (water temperature,oxygen)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_O2o,O2o)
!         _GET_(self%id_O3c,O3c)

    ! Get prey concentrations
    do iprey = 1, self%nprey
      _GET_(self%id_preyc(iprey), preycP(iprey))
      _GET_(self%id_preyn(iprey), preynP(iprey))
      _GET_(self%id_preyp(iprey), preypP(iprey))
      _GET_(self%id_preyl(iprey), preylP(iprey))
      _GET_(self%id_preys(iprey), preysP(iprey))
!#ifdef INCLUDE_PELFE
      ! _GET_(self%id_preyf(iprey), preyfP(iprey))
!#endif
    enddo


!! Quota collectors
!  tfluxC = ZERO
!  tfluxN = ZERO
!  tfluxP = ZERO
   qpcMEZ = zoop/(zooc+p_small)      ! add some epsilon (add in shared) to avoid divide by 0
   qncMEZ = zoon/(zooc+p_small)

      _SET_DIAGNOSTIC_(self%id_qncMEZ,qncMEZ)
      _SET_DIAGNOSTIC_(self%id_qpcMEZ,qpcMEZ)

      _SET_DIAGNOSTIC_(self%id_prey2l,preylP(2))
      _SET_DIAGNOSTIC_(self%id_prey1x,preysP(1))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature and oxygen response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   eo = MM_power(max(p_small,O2o), self%p_clO2o, 3)
   et = eTq(ETW, self%p_q10)
      _SET_DIAGNOSTIC_(self%id_ETWd,ETW)
      _SET_DIAGNOSTIC_(self%id_et,et)
      _SET_DIAGNOSTIC_(self%id_eo,eo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total potential food given the non-dim prey availability
  ! with loops over all LFGs.
  ! There is no parameter for capture efficiency in mesozooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  rumc = ZERO
!  do i = 1, iiPhytoPlankton
!    PPYc(:,i) = p_paPPY(zoo,i)*PhytoPlankton(i,iiC)
!    rumc = rumc + PPYc(:,i)
!  end do
!  do i = 1, iiMicroZooPlankton
!    MIZc(:,i) = p_paMIZ(zoo,i)*MicroZooPlankton(i,iiC)
!    rumc = rumc + MIZc(:,i)
!  end do
!  do i = 1, iiMesoZooPlankton
!    MEZc(:,i) = p_paMEZ(zoo,i)*MesoZooPlankton(i,iiC)
!    rumc = rumc + MEZc(:,i)
!  end do

    rumc   = ZERO
    do iprey = 1, self%nprey
      PPYc(iprey) = self%p_pa(iprey)*preycP(iprey)
      rumc = rumc + PPYc(iprey)
    end do
    
      _SET_DIAGNOSTIC_(self%id_rumc,rumc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake rate (eq 38 Vichi et al. 2007) and 
  ! specific uptake rate considering potentially available food (sut)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   rugc  = et*self%p_sum*MM(self%p_vum*rumc, self%p_sum)*zooc
   sut = rugc/(p_small + rumc)

      _SET_DIAGNOSTIC_(self%id_rugc,rugc)
      _SET_DIAGNOSTIC_(self%id_sut,sut)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptakes from every LFG
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   rut_c = ZERO
   rut_n = ZERO
   rut_p = ZERO

  ! Phytoplankton
!  do i = 1, iiPhytoPlankton
    do iprey = 1, self%nprey
     ruPPYc = sut*PPYc(iprey)
!    call flux_vector(iiPel, ppPhytoPlankton(i,iiC), ppzooc, ruPPYc)
    _SET_ODE_(self%id_c,             ruPPYc)
    _SET_ODE_(self%id_preyc(iprey), -ruPPYc)
!    call flux_vector(iiPel, ppPhytoPlankton(i,iiN), ppzoon, ruPPYc*qncPPY(i,:))
    _SET_ODE_(self%id_n,             ruPPYc*(preynP(iprey)/(preycP(iprey)+p_small)))
    _SET_ODE_(self%id_preyn(iprey), -ruPPYc*(preynP(iprey)/(preycP(iprey)+p_small)))
!    call flux_vector(iiPel, ppPhytoPlankton(i,iiP), ppzoop, ruPPYc*qpcPPY(i,:))
    _SET_ODE_(self%id_p,             ruPPYc*(preypP(iprey)/(preycP(iprey)+p_small)))
    _SET_ODE_(self%id_preyp(iprey), -ruPPYc*(preypP(iprey)/(preycP(iprey)+p_small)))

     rut_c = rut_c + ruPPYc
     rut_n = rut_n + ruPPYc*(preynP(iprey)/(preycP(iprey)+p_small))
     rut_p = rut_p + ruPPYc*(preypP(iprey)/(preycP(iprey)+p_small))

!    if (self%p_pl(iprey)) then
    ! Chl is transferred to the infinite sink
!    call flux_vector(iiPel, ppPhytoPlankton(i,iiL), ppPhytoPlankton(i,iiL), -ruPPYc*qlcPPY(i,:))
    _SET_ODE_(self%id_preyl(iprey), -ruPPYc*(preylP(iprey)/(preycP(iprey)+p_small)))
!    end if

      _SET_DIAGNOSTIC_(self%id_ruPPYc(iprey), ruPPYc)
      _SET_DIAGNOSTIC_(self%id_preyld(iprey), preylP(iprey))


!    if (self%p_ps(iprey)) then
    ! silicon constituent is transferred to biogenic silicate
!    if ( ppPhytoPlankton(i,iiS) .gt. 0 ) & 
!       call flux_vector(iiPel, ppPhytoPlankton(i,iiS), ppR6s, ruPPYc*qscPPY(i,:))
    _SET_ODE_(self%id_R6s,           ruPPYc*(preysP(iprey)/(preycP(iprey)+p_small)))
    _SET_ODE_(self%id_preys(iprey), -ruPPYc*(preysP(iprey)/(preycP(iprey)+p_small)))
!    end if

!#ifdef INCLUDE_PELFE
!    ! Fe constituent is transferred to particulate iron
!    if ( ppPhytoPlankton(i,iiF) .gt. 0 ) & 
!       call flux_vector(iiPel, ppPhytoPlankton(i,iiF), ppR6f, ruPPYc*qfcPPY(i,:))
!#endif

!#if defined INCLUDE_PELCO2
! PIC (calcite/aragonite) production associated to the grazed biomass
! The idea in PISCES is that the calcite flux exists only when associated
! to a carbon release from phytoplankton (there is no calcite storage in phyto)
! Use the realized rain ratio for each phytoplankton species and assume
! that only a portion is egested
! Calcite production is parameterized as a flux between DIC and PIC
! that affects alkalinity
   if (self%p_isP2(iprey).eq.1) then
!    call flux_vector( iiPel, ppO3c,ppO5c, p_pecaco3(zoo)*ruPPYc*qccPPY(i,:))
    _SET_ODE_(self%id_O3c,-self%p_pecaco3*ruPPYc*qccPPY)        ! precipitation of CaCO3 consumes DIC
    _SET_ODE_(self%id_O5c,self%p_pecaco3*ruPPYc*qccPPY)         ! precipitation of CaCO3 produces PIC
!    call flux_vector( iiPel, ppO3h,ppO3h, -C2ALK*p_pecaco3(zoo)*ruPPYc*qccPPY(i,:))
    _SET_ODE_(self%id_O3h,-C2ALK*self%p_pecaco3*ruPPYc*qccPPY)  ! precipitation of CaCO3 consumes 2 alkalinity
   
    _SET_DIAGNOSTIC_(self%id_CaCO3precip(iprey), self%p_pecaco3*ruPPYc*qccPPY)
    _SET_DIAGNOSTIC_(self%id_CaCO3_to_O3h(iprey),-C2ALK*self%p_pecaco3*ruPPYc*qccPPY)
   endif
!#endif

     end do

      _SET_DIAGNOSTIC_(self%id_rut_c,rut_c)
      _SET_DIAGNOSTIC_(self%id_rut_n,rut_n)
      _SET_DIAGNOSTIC_(self%id_rut_p,rut_p)


  ! Microzooplankton
!  do i = 1, iiMicroZooPlankton
!    ruMIZc = sut*MIZc(:,i)
!    call flux_vector(iiPel, ppMicroZooPlankton(i,iiC), ppzooc, ruMIZc)
!    call flux_vector(iiPel, ppMicroZooPlankton(i,iiN), ppzoon, ruMIZc*qncMIZ(i,:))
!    call flux_vector(iiPel, ppMicroZooPlankton(i,iiP), ppzoop, ruMIZc*qpcMIZ(i,:))
!    call flux_vector(iiPel, ppMicroZooPlankton(i,iiH), ppzooh, ruMIZc*qhcMIZ(i,:))
!    rut_c = rut_c + ruMIZc
!    rut_n = rut_n + ruMIZc*qncMIZ(i,:)
!    rut_p = rut_p + ruMIZc*qpcMIZ(i,:)
!  end do

  ! Mesozooplankton
!  do i = 1, iiMesoZooPlankton
!    ruMEZc = sut*MEZc(:, i)
!    ! Note that intra-group predation (cannibalism) is not added as a flux
!    if ( i/= zoo ) then
!      call flux_vector(iiPel, ppMesoZooPlankton(i,iiC), ppzooc, ruMEZc)
!      call flux_vector(iiPel, ppMesoZooPlankton(i,iiN), ppzoon, ruMEZc*qncMEZ(i,:))
!      call flux_vector(iiPel, ppMesoZooPlankton(i,iiP), ppzoop, ruMEZc*qpcMEZ(i,:))
!      call flux_vector(iiPel, ppMesoZooPlankton(i,iiH), ppzooh, ruMEZc*qhcMEZ(i,:))
!    end if
!    rut_c = rut_c + ruMEZc
!    rut_n = rut_n + ruMEZc*qncMEZ(i,:)
!    rut_p = rut_p + ruMEZc*qpcMEZ(i,:)
!  end do


  ! Note that tfluxC include also intra-group predation
!  tfluxC = tfluxC + rut_c 
!  tfluxN = tfluxN + rut_n
!  tfluxP = tfluxP + rut_p


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity respiration and basal metabolism
  ! First compute the the energy cost of ingestion
  ! 1 - assimilation - egestion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   prI = ONE - self%p_puI - self%p_peI
   rrc = prI * rut_c + self%p_srs*et*zooc
!  call flux_vector(iiPel, ppO2o, ppO2o, -rrc/MW_C)
   _SET_ODE_(self%id_O2o,-(rrc/MW_C))
!  call quota_flux(iiPel, ppzooc, ppzooc, ppO3c, rrc, tfluxC)
   _SET_ODE_(self%id_c,-rrc)
   _SET_ODE_(self%id_O3c,rrc)

   _SET_DIAGNOSTIC_(self%id_rrc,rrc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Specific rates of low oxygen mortality
  ! and Density dependent mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   rdo_c = self%p_sdo*(ONE-eo)*et*zooc
   rd_c  = self%p_sd*zooc**self%p_sds

      _SET_DIAGNOSTIC_(self%id_rdo_c,rdo_c)
      _SET_DIAGNOSTIC_(self%id_rd_c, rd_c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total egestion including pellet production 
  ! Eq. 40 and 44 Vichi et al. 2007
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   rq6c = self%p_peI*rut_c + rdo_c + rd_c
   rq6n = self%p_peI*rut_n + qncMEZ*(rdo_c + rd_c)
   rq6p = self%p_peI*rut_p + qpcMEZ*(rdo_c + rd_c)

      _SET_DIAGNOSTIC_(self%id_rq6c,rq6c)
      _SET_DIAGNOSTIC_(self%id_rq6n,rq6n)
      _SET_DIAGNOSTIC_(self%id_rq6p,rq6p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient remineralization 
  ! basal metabolism + excess of non-limiting nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ren = self%p_srs*et*eo*zoon 
   rep = self%p_srs*et*eo*zoop 
      _SET_DIAGNOSTIC_(self%id_ren,ren)
      _SET_DIAGNOSTIC_(self%id_rep,rep)

!  call quota_flux(iiPel, ppzoop, ppzoop, ppN1p, rep, tfluxP)
   _SET_ODE_(self%id_p,  -rep)
   _SET_ODE_(self%id_N1p, rep)
!   _SET_ODE_(self%id_O3h, -rep)      ! release of 1 PO4 decreases 1 alkalinity

!  call quota_flux(iiPel, ppzoon, ppzoon, ppN4n, ren, tfluxN)
   _SET_ODE_(self%id_n,  -ren)
   _SET_ODE_(self%id_N4n, ren)
   _SET_ODE_(self%id_O3h, ren)      ! release of 1 NH4 increases 1 alkalinity

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes to particulate organic matter
  ! Add the correction term for organic carbon release in case of
  ! nutrient limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  call quota_flux(iiPel, ppzooc, ppzooc,ppR6c, rq6c, tfluxC)
   _SET_ODE_(self%id_c,  -rq6c)
   _SET_ODE_(self%id_R6c, rq6c)
!  call quota_flux(iiPel, ppzoon, ppzoon,ppR6n, rq6n, tfluxN)
   _SET_ODE_(self%id_n,  -rq6n)
   _SET_ODE_(self%id_R6n, rq6n)
!  call quota_flux(iiPel, ppzoop, ppzoop,ppR6p, rq6p, tfluxP)
   _SET_ODE_(self%id_p,  -rq6p)
   _SET_ODE_(self%id_R6p, rq6p)


!  if ( zoon > 0 .and. zoop > 0 ) then
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Check the assimilation rate for Carbon, Nitrogen and Phosphorus
     ! Note that activity respiration does not involve nutrient utilization
     ! so more nutrients than carbon are taken up.
     ! Then compute P:C and N:C ratios in the assimilation rate
     ! Eq 41 in Vichi et al. 2007 (there is an error in the denominator,
     ! the \Iota_c should be \Iota_i, with i=n,p)
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ru_c = self%p_puI*rut_c
      ru_n = (self%p_puI + prI)* rut_n
      ru_p = (self%p_puI + prI)* rut_p
      pu_e_n  =   ru_n/( p_small+ ru_c)
      pu_e_p  =   ru_p/( p_small+ ru_c)

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Eliminate the excess of the non-limiting constituent
     ! Determine whether C, P or N is the limiting element and assign the
     ! value to variable limit
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      limit = 1.0_rk
      !limit = c
      temp_p  = pu_e_p/qpcMEZ
      temp_n  = pu_e_n/qncMEZ

      _SET_DIAGNOSTIC_(self%id_temp_p, temp_p)
      _SET_DIAGNOSTIC_(self%id_temp_n, temp_n)

      if ( temp_p<temp_n .OR. abs(temp_p-temp_n)<p_small ) then 
          if ( pu_e_p< qpcMEZ ) then
            limit = 3.0_rk
            !limit = p
         end if
      else
          if ( pu_e_n<qncMEZ ) then
            limit = 2.0_rk
            !limit = n
          end if
      end if

      _SET_DIAGNOSTIC_(self%id_limit,  limit)


     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Compute the correction terms depending on the limiting constituent
     ! Eq. 42 Vichi et al 2007 for a combination of N and P limitation
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if ( limit == 1.0_rk ) then
!      if ( limit == iiC ) 
         pe_R6c = ZERO
         pe_N1p = max(ZERO, (ONE - self%p_peI)*rut_p - self%p_qpcMEZ*ru_c)
         pe_N4n = max(ZERO, (ONE - self%p_peI)*rut_n - self%p_qncMEZ*ru_c)
      else if ( limit == 3.0_rk ) then
!      else if ( limit == iiP )
         pe_R6c = max(ZERO, ru_c - (ONE - self%p_peI)*rut_p/self%p_qpcMEZ)
         pe_N1p = ZERO
         pe_N4n = max( ZERO, (ONE - self%p_peI)*rut_n - self%p_qncMEZ*(ru_c - pe_R6c))
      else if ( limit == 2.0_rk ) then
!      else if ( limit == iiN )
         pe_R6c = max(ZERO, ru_c - (ONE - self%p_peI)*rut_n/self%p_qncMEZ)
         pe_N1p = max(ZERO, (ONE - self%p_peI)*rut_p - self%p_qpcMEZ*(ru_c - pe_R6c))
         pe_N4n = ZERO
      end if   

      _SET_DIAGNOSTIC_(self%id_pe_R6c, pe_R6c)
      _SET_DIAGNOSTIC_(self%id_pe_N1p, pe_N1p)
      _SET_DIAGNOSTIC_(self%id_pe_N4n, pe_N4n)


!  else
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Eliminate the excess of the non-limiting constituent under fixed quota
     ! Determine whether C, P or N is limiting (Total Fluxes Formulation)
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!     limit = nutlim(tfluxc,tfluxn,tfluxp,qncMEZ(zoo,:),qpcMEZ(zoo,:),iiC,iiN,iiP)
   
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Compute the correction terms depending on the limiting constituent
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!     WHERE     ( limit == iiC )
!         pe_R6c = ZERO
!         pe_N1p = max(ZERO,tfluxp  - p_qpcMEZ(zoo)* tfluxc)
!         pe_N4n = max(ZERO,tfluxn  - p_qncMEZ(zoo)* tfluxc)
!     ELSEWHERE ( limit == iiP )
!         pe_N1p = ZERO
!         pe_R6c = max(ZERO, tfluxc  - tfluxp/p_qpcMEZ(zoo))
!         pe_N4n = max(ZERO, tfluxn  - tfluxp/p_qpcMEZ(zoo)*p_qncMEZ(zoo) )
!     ELSEWHERE ( limit == iiN )
!         pe_N4n = ZERO
!         pe_R6c = max(ZERO, tfluxc  - tfluxn/p_qncMEZ(zoo))
!         pe_N1p = max(ZERO, tfluxp  - tfluxn/p_qncMEZ(zoo)*p_qpcMEZ(zoo))
!     END WHERE

!#ifdef DEBUG
!     write(*,*) '+++++++++++++++'
!     if ( limit(1)==iiC ) then
!     write(*,*) 'tfluxp', tfluxp,'pe_N1p', tfluxp  - p_qpcMEZ(zoo)* tfluxc 
!     write(*,*) 'tfluxn', tfluxn,'pe_N4n', tfluxn  - p_qncMEZ(zoo)* tfluxc
!     write(*,*) 'tfluxc', tfluxc,'pe_R6c', ZERO 
!     write(*,*) 'ooooooooooooooo'
!     endif
   
!     if ( limit(1)==iiP ) then
!     write(*,*) 'tfluxp', tfluxp,'pe_N1p', ZERO 
!     write(*,*) 'tfluxn', tfluxn,'pe_N4n', tfluxn - tfluxp/p_qpcMEZ(zoo)*p_qncMEZ(zoo)
!     write(*,*) 'tfluxc', tfluxc,'pe_R6c', tfluxc - tfluxp/p_qpcMEZ(zoo)
!     write(*,*) 'ooooooooooooooo'
!     endif
   
!     if ( limit(1)==iiN ) then
!     write(*,*) 'tfluxp', tfluxp,'pe_N1p', tfluxp  - tfluxn/p_qncMEZ(zoo)*p_qpcMEZ(zoo)
!     write(*,*) 'tfluxn', tfluxn,'pe_N4n', ZERO
!     write(*,*) 'tfluxc', tfluxc,'pe_R6c', tfluxc  - tfluxn/p_qncMEZ(zoo) 
!     endif
!     write(*,*) '+++++++++++++++'
!#endif
!   
!  endif



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correction term for excess of non-limiting nutrients as organic carbon 
  ! release (POC) and nutrient remineralization (PO4 and NH)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  call flux_vector(iiPel, ppzooc, ppR6c, pe_R6c)
   _SET_ODE_(self%id_c,  -pe_R6c)
   _SET_ODE_(self%id_R6c, pe_R6c)
!  call flux_vector(iiPel, ppzoop, ppN1p, pe_N1p)
   _SET_ODE_(self%id_p,  -pe_N1p)
   _SET_ODE_(self%id_N1p, pe_N1p)
!   _SET_ODE_(self%id_O3h, -pe_N1p)      ! release of 1 PO4 decrease alkalinity
!  call flux_vector(iiPel, ppzoon, ppN4n, pe_N4n)
   _SET_ODE_(self%id_n,  -pe_N4n)
   _SET_ODE_(self%id_N4n, pe_N4n)
   _SET_ODE_(self%id_O3h, pe_N4n)      ! release of 1 NH4 increase alkalinity


   _SET_DIAGNOSTIC_(self%id_varO3h_Nutil, pe_N4n + ren )
   _SET_DIAGNOSTIC_(self%id_varO3h_Putil, -rep )
    _LOOP_END_
  end subroutine do
end module


!  end subroutine MesoZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
