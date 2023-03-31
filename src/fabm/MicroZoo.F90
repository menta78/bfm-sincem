#include "fabm_driver.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MicroZoo
!
! DESCRIPTION
!  Microzooplankton dynamics
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

  module bfm_MicroZoo

    use fabm_types
    use fabm_particle

    use ogs_bfm_shared
    use ogs_bfm_pelagic_base
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   use global_mem, ONLY:RLEN,ZERO,ONE,rnd_SEED
!   use constants,  ONLY:MW_C,C2ALK
! #ifdef NOPOINTERS
!   use mem
! #else
!   use mem, ONLY: D3STATE, O2o, R1c, R6c, R1n, R6n, &
!     R2c,R1p, R6p, N4n, N1p, PhytoPlankton, MicroZooPlankton, PelBacteria
!   use mem, ONLY: ppPelBacteria, ppO2o, ppR1c, ppR1l, ppR6c, ppR6s, Depth,&
!     ppR1n, ppR6n, ppR1p, ppR1l, ppR6p, ppN4n, ppN1p, ppPhytoPlankton, ppMicroZooPlankton, &
!     ETW, eO2mO2, qncPBA, qpcPBA, qncPPY, qpcPPY, qncMIZ, qpcMIZ, iiPelBacteria, &
!     qlcPPY, qscPPY, iiPhytoPlankton, iiMicroZooPlankton, iiC, iiN, iiP, iiL, iiS, &
!     NO_BOXES, iiBen, iiPel, flux_vector, quota_flux
! #ifdef INCLUDE_PELCO2
!   use mem, ONLY: ppO3c, ppO5c, ppO3h, qccPPY
! #endif
! #ifdef INCLUDE_PELFE
!   use mem, ONLY: iiF, qfcPPY, ppR6f
! #endif
! #endif
!   use mem_Param,     ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_small
!   use bfm_error_msg, ONLY: bfm_error
!   use mem_MicroZoo
!   use standalone, only: delt,maxdelt
!   use random_generator
!   use gaussian_generator


!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   ! The following vector functions are used: eTq, MM, nutlim
!   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   use mem_globalfun,   ONLY: eTq, MM, nutlim

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  implicit none

  private

!
! !AUTHORS
!   First ERSEM version by H. Baretta-Bekker and J.W. Baretta
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  type,extends(type_ogs_bfm_pelagic_base),public :: type_ogs_bfm_microzoo
      ! NB: own state variables (c,n,p,s,f,chl) are added implicitly by deriving
      ! from type_ogs_bfm_pelagic_base!
      ! Identifiers for preys
  
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyc
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyn
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyp
      type (type_state_variable_id), allocatable,dimension(:) :: id_preyl
      type (type_state_variable_id), allocatable,dimension(:) :: id_preys
      ! type (type_state_variable_id),    allocatable,dimension(:) :: id_preyf
      type (type_model_id),      allocatable,dimension(:) :: id_prey



      ! Identifiers for state variables of other models
      ! type (type_state_variable_id) :: id_O3c,id_O2o,id_O3h                  !  dissolved inorganic carbon, oxygen, total alkalinity
      type (type_state_variable_id)      :: id_O2o   ! oxygen
      type (type_state_variable_id)      :: id_O3c   ! dissolved inorganic carbon
      type (type_state_variable_id)      :: id_O5c   ! particulate inorganic carbon
      type (type_state_variable_id)      :: id_O3h   ! alkalinity
      type (type_state_variable_id)      :: id_N4n   ! Ammonium
      type (type_state_variable_id)      :: id_N1p   ! Phosphate
      type (type_state_variable_id)      :: id_R6c,id_R6s,id_R6n,id_R6p   ! particulate organic carbon, silicon, nitrogen, phosphorous
      type (type_state_variable_id)      :: id_X1c   ! colored dissolved organic carbon
      type (type_state_variable_id)      :: id_R1c,id_R1n,id_R1p   ! dissolved organic carbon, nitrogen, phosphorous (R1: labile)

      ! Environmental dependencies
      type (type_dependency_id)          :: id_ETW   ! temperature

      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id) :: id_ETWd   ! temperature Celsius
      type (type_diagnostic_variable_id) :: id_et     ! temperature q10 factor
      type (type_diagnostic_variable_id) :: id_eO2    ! Oxygen limitation
      type (type_diagnostic_variable_id) :: id_rumc   ! total potential food
      type (type_diagnostic_variable_id) :: id_rugc   ! total food uptake rate (eq 38 Vichi et al. 2007)
      type (type_diagnostic_variable_id) :: id_sut    ! specific uptake rate considering potentially available food
      type (type_diagnostic_variable_id) :: id_rugn   ! tbd
      type (type_diagnostic_variable_id) :: id_rugp   ! tbd
      type (type_diagnostic_variable_id) :: id_rrtc   ! tbd
      type (type_diagnostic_variable_id) :: id_rrsc   ! tbd
      type (type_diagnostic_variable_id) :: id_rrac   ! tbd
      type (type_diagnostic_variable_id) :: id_rric   ! tbd
      type (type_diagnostic_variable_id) :: id_reac   ! tbd
      type (type_diagnostic_variable_id) :: id_rdc    ! tbd
      type (type_diagnostic_variable_id) :: id_rr1c   ! exudation flux to DOC
      type (type_diagnostic_variable_id) :: id_rr6c   ! tbd
      type (type_diagnostic_variable_id) :: id_rrin   ! tbd
      type (type_diagnostic_variable_id) :: id_rr1n   ! tbd
      type (type_diagnostic_variable_id) :: id_rr6n   ! tbd
      type (type_diagnostic_variable_id) :: id_rrip   ! tbd
      type (type_diagnostic_variable_id) :: id_rr1p   ! tbd
      type (type_diagnostic_variable_id) :: id_rr6p   ! tbd
      type (type_diagnostic_variable_id) :: id_runc   ! tbd
      type (type_diagnostic_variable_id) :: id_runn   ! tbd
      type (type_diagnostic_variable_id) :: id_runp   ! tbd
      type (type_diagnostic_variable_id) :: id_ren    ! tbd
      type (type_diagnostic_variable_id) :: id_rep    ! tbd
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preycd !prey c
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preynd !prey n
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preypd !prey p
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preysd !prey s
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preyfd !prey f
      ! type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_preydld !prey l

      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_CaCO3precip ! precipitation of PIC
      type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_CaCO3_to_O3h ! consume of alk due to precipitation of PIC
      type (type_diagnostic_variable_id) :: id_varO3h_Nutil ! variation of O3h due to NH44 utilization by Zoo
      type (type_diagnostic_variable_id) :: id_varO3h_Putil ! variation of O3h due to PO4 utilization by Zoo

      ! Parameters (described in subroutine initialize, below)
      integer  :: nprey
      real(rk), allocatable :: p_pa(:)
      integer, allocatable :: p_isP2(:)
      real(rk) :: p_q10, p_srs, p_sum, p_sdo, p_sd
      real(rk) :: p_pu, p_pu_ea, p_chro, p_chuc, p_minfood
      real(rk) :: p_pecaco3, p_qncMIZ, p_qpcMIZ
      real(rk) :: p_pe_R1c, p_pe_R1n, p_pe_R1p
      real(rk) :: p_fX1z
      
      ! Parameters (described in subroutine initialize, below)
  ! integer       :: i
      ! integer       :: ppzooc, ppzoon, ppzoop
      ! integer, save :: first =0
      ! integer(4)    :: local_seed
      ! integer       :: AllocStatus
      ! integer,dimension(NO_BOXES)  :: limit
      ! real(RLEN),allocatable,save,dimension(:) :: sut,et,eO2,rumc,  &
      ! rugc,rugn,rugp,runc,runn,runp, &
      ! rrsc,rrac,reac,rdc,rrtc,ruPBAc,ruPPYc,  &
      ! ruMIZc,rric,rr1c,rr6c,rr1p,rr1n, &
      ! rrip,rr6p,rep,rrin,zooc, tfluxC, tfluxN, tfluxP,MIZ_FLUCT
  ! real(RLEN),allocatable,save,dimension(:)    :: rr6n,ren,pu_ra,r
      ! real(RLEN),allocatable,save,dimension(:)    :: pe_N1p, pe_N4n, pe_R6c
  ! real(RLEN),allocatable,save,dimension(:,:)  :: PBAc,PPYc,MIZc
      ! #ifndef INCLUDE_PELCO2
      ! integer,parameter :: ppO3c = 0
      ! #endif
      contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
      
      end type type_ogs_bfm_microzoo
      
      contains
      
      subroutine initialize(self,configunit)
        !
        ! !DESCRIPTION:
        !
        ! !INPUT PARAMETERS:
        class (type_ogs_bfm_microzoo),intent(inout),target :: self
        integer,                      intent(in)           :: configunit
        !
        ! !REVISION HISTORY:
        !
        ! !LOCAL VARIABLES:
        logical  :: preyisdiat,preyisphyto
        integer           :: iprey
        character(len=16) :: index
        real(rk) :: pippo1
        !EOP
        !-----------------------------------------------------------------------
        !BOC
        ! Obtain the values of all model parameters from FABM.
        ! Specify the long name and units of the parameters, which could be used
        ! by FABM (or its host)
        ! to present parameters to the user for configuration (e.g., through a
        ! GUI)
        call self%get_parameter(self%p_q10,     'p_q10',     '-',         'Characteristic Q10 coefficient')
        call self%get_parameter(self%p_srs,     'p_srs',     '1/d',       'Respiration rate at 10 degrees Celsius')
        call self%get_parameter(self%p_sum,     'p_sum',     '1/d',       'Potential growth rate')
        call self%get_parameter(self%p_sdo,     'p_sdo',     '1/d',       'Mortality rate due to oxygen limitation')
        call self%get_parameter(self%p_sd,      'p_sd',      '1/d',       'Temperature independent mortality rate')
        call self%get_parameter(self%p_pu,      'p_pu',      '-',         'Assimilation efficiency')
        call self%get_parameter(self%p_pu_ea,   'p_pu_ea',   '-',         'Fraction of activity excretion')
        call self%get_parameter(self%p_chro,    'p_chro',    'mmol/m3',   'Half-saturation oxygen concentration')
        call self%get_parameter(self%p_chuc,    'p_chuc',    'mgC/m3',    'Half-saturation Food concentration for Type II')
        call self%get_parameter(self%p_minfood, 'p_minfood', 'mgC/m3',    'Half-saturation food concentration for preference factor')
        call self%get_parameter(self%p_pecaco3, 'p_pecaco3', '-',         'Portion of egested calcified shells during grazing')
        call self%get_parameter(self%p_qncMIZ,  'p_qncMIZ',  'mmolN/mgC', 'Maximum quotum P:C')
        call self%get_parameter(self%p_qpcMIZ,  'p_qpcMIZ',  'mmolN/mgC', 'Maximum quotum N:C')
        call self%get_parameter(self%p_pe_R1c,  'p_pe_R1c',  '-',         'Fractional content of C in cytoplasm')
        call self%get_parameter(self%p_pe_R1n,  'p_pe_R1n',  '-',         'Fractional content of N in cytoplasm')
        call self%get_parameter(self%p_pe_R1p,  'p_pe_R1p',  '-',         'Fractional content of P in cytoplasm')
!              --------- Flux partition CDOM parameters ------------
        call self%get_parameter(self%p_fX1z,    'p_fX1z',    '-',         'colored fraction in labile DOC', default=0.02_rk)
        
        ! Register state variables (handled by type_bfm_pelagic_base)
        call self%initialize_bfm_base()
        call self%add_constituent('c',1.e-4_rk)
        call self%add_constituent('n',1.26e-6_rk)
        call self%add_constituent('p',4.288e-8_rk)
        
        ! Determine number of prey types
        call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
        ! Get prey-specific parameters.
        allocate(self%p_pa(self%nprey)) !Availability of nprey for microzoo group
        allocate(self%p_isP2(self%nprey))   !is P2? [=1 for P2 and 0 otherwise]
        allocate(self%id_prey(self%nprey))
        allocate(self%id_preyc(self%nprey))
        allocate(self%id_preyn(self%nprey))
        allocate(self%id_preyp(self%nprey))
        allocate(self%id_preyl(self%nprey))
        allocate(self%id_preys(self%nprey))


        allocate(self%id_CaCO3precip(self%nprey))
        allocate(self%id_CaCO3_to_O3h(self%nprey))
  
        do iprey=1,self%nprey
          write (index,'(i0)') iprey
          call self%get_parameter(self%p_pa(iprey),'suprey'//trim(index),'-','Availability for prey type '//trim(index))
          call self%get_parameter(self%p_isP2(iprey),'isP2'//trim(index),'-','identify P2 among the preys '//trim(index))
          
          call self%register_state_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mg C/m^3', 'prey '//trim(index)//' carbon') ! mg o mml? anna
          call self%register_state_dependency(self%id_preyn(iprey),'prey'//trim(index)//'n','mmol n/m^3', 'prey '//trim(index)//' nitrogen')
          call self%register_state_dependency(self%id_preyp(iprey),'prey'//trim(index)//'p','mmol p/m^3', 'prey '//trim(index)//' phosphorous')
          call self%register_state_dependency(self%id_preys(iprey),'prey'//trim(index)//'s','mmol Si/m^3', 'prey '//trim(index)//' silica')
          call self%register_state_dependency(self%id_preyl(iprey),'prey'//trim(index)//'Chl','mg Chl/^3', 'prey '//trim(index)//' chlorophyll')
!          call self%register_state_dependency(self%id_preyf(iprey),'prey'//trim(index)//'f','umol Fe/^3', 'prey '//trim(index)//' iron')
          
          call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
          call self%request_coupling_to_model(self%id_preyc(iprey),self%id_prey(iprey),'c')
          call self%request_coupling_to_model(self%id_preyn(iprey),self%id_prey(iprey),'n')
          call self%request_coupling_to_model(self%id_preyp(iprey),self%id_prey(iprey),'p')
          call self%request_coupling_to_model(self%id_preys(iprey),self%id_prey(iprey),standard_variables%total_silicate)
          call self%request_coupling_to_model(self%id_preyl(iprey),self%id_prey(iprey),total_chlorophyll)
!          call self%request_coupling_to_model(self%id_preyf(iprey),self%id_prey(iprey),'f')
        end do
        
        
        ! Register links to external nutrient pools.
        call self%register_state_dependency(self%id_O2o,'O2o','mmol O2/m^3','dissolved oxygen')
        call self%register_state_dependency(self%id_O3c,'O3c','mg C/m^3',   'dissolved organic carbon')
        call self%register_state_dependency(self%id_O5c,'O5c','mg C/m^3',   'particulate inorganic carbon')
        call self%register_state_dependency(self%id_O3h,'O3h','mmol/m^3',   'alkalinity')
        call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3', 'ammonium')
        call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3', 'phosphate')
        call self%register_state_dependency(self%id_R6c,'R6c','mmg C/m^3',  'POC')
        call self%register_state_dependency(self%id_R6s,'R6s','mg Si/m^3',  'POS')
        call self%register_state_dependency(self%id_R6n,'R6n','mmol N/m^3', 'PON')
        call self%register_state_dependency(self%id_R6p,'R6p','mmol P/m^3', 'POP')
        call self%register_state_dependency(self%id_X1c,'X1c','mg C/m^3',   'labile CDOM')
        call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3',   'labile DOC')
        call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3', 'labile DON')
        call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3', 'labile DOP')
        
        ! Register environmental dependencies (temperature, shortwave radiation)
        call self%register_dependency(self%id_ETW,standard_variables%temperature)
        
        ! Register diagnostic variables (i.e., model outputs)
        call self%register_diagnostic_variable(self%id_ETWd, 'ETW',  'C',     'temperature Celsius',output=output_none)
        call self%register_diagnostic_variable(self%id_et,   'et',   '-',     'temperature factor',output=output_none)
        call self%register_diagnostic_variable(self%id_eO2,  'eO2',  '-',     'Oxygen limitation',output=output_none)
        call self%register_diagnostic_variable(self%id_rumc, 'rumc', 'mgC/m3',   'total potential food',output=output_none)
        call self%register_diagnostic_variable(self%id_rugc, 'rugc', 'mgC/m3/d', 'total food uptake rate',output=output_none)
        call self%register_diagnostic_variable(self%id_sut,  'sut',  '1/d',      'specific uptake rate',output=output_none)
        call self%register_diagnostic_variable(self%id_rugn, 'rugn', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rugp, 'rugp', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rrtc, 'rrtc', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rrsc, 'rrsc', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rrac, 'rrac', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rric, 'rric', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_reac, 'reac', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rdc,  'rdc',  'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rr1c, 'rr1c', 'mgC/m3/d', 'exudation flux to DOC')
        call self%register_diagnostic_variable(self%id_rr6c, 'rr6c', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rrin, 'rrin', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rr1n, 'rr1n', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rr6n, 'rr6n', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rrip, 'rrip', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rr1p, 'rr1p', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rr6p, 'rr6p', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_runc, 'runc', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_runn, 'runn', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_runp, 'runp', 'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_ren,  'ren',  'tbd',      'tbd',output=output_none)
        call self%register_diagnostic_variable(self%id_rep,  'rep',  'tbd',      'tbd',output=output_none)

       call self%register_diagnostic_variable(self%id_varO3h_Nutil,'varO3h_Nutil','mmol/m3/d','variaz O3h due to N utiliz',output=output_none)
       call self%register_diagnostic_variable(self%id_varO3h_Putil,'varO3h_Putil','mmol/m3/d','variazO3h due to P utiliz',output=output_none)
      do iprey=1,self%nprey
       write (index,'(i0)') iprey
       if (self%p_isP2(iprey).eq.1) then
         call self%register_diagnostic_variable(self%id_CaCO3precip(iprey),'_'//trim(index)//'_CaCO3precip','mgC/m^3/d','prey '//trim(index)//' CaCO3precip',output=output_none)
         call self%register_diagnostic_variable(self%id_CaCO3_to_O3h(iprey),'_'//trim(index)//'_consumeO3h_for_CaCO3precip','mmol/m^3/d','prey '//trim(index)//' consumeO3h_for_CaCO3precip',output=output_none)
       endif
      end do

     end subroutine
    
    subroutine do(self,_ARGUMENTS_DO_)
      
      class (type_ogs_bfm_microzoo),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      !LOCAL VARIABLES:
      integer  :: iprey,istate
      real(rk), dimension(self%nprey) :: preycP,preypP,preynP,preylP,preysP
      real(rk), dimension(self%nprey) :: PPYc,qpcPPY,qncPPY,qlcPPY,qscPPY
      real(rk) :: preyP
      real(rk) :: zooc, zoop, zoon
      real(rk) :: qncMIZ, qpcMIZ
      real(rk) :: et,ETW,eO2
      real(rk) :: O2o
      real(rk) :: rumc,rugc,sut
      real(rk) :: rugn,rugp,ruPPYc
      real(rk) :: rrtc,rrsc,rrac
      real(rk) :: rric,reac,rdc,rr1c,rr6c
      real(rk) :: rrin,rr1n,rr6n
      real(rk) :: rrip,rr1p,rr6p
      real(rk) :: runc,runn,runp,ren,rep
      
      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Allocate local memory
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! if (first==0) then
      !    ALLOCATE ( PBAc(NO_BOXES,iiPelBacteria),   PPYc(NO_BOXES,iiPhytoPlankton),  &
      !       &       MIZc(NO_BOXES,iiMicroZooPlankton),               &
      !       &       sut(NO_BOXES), et(NO_BOXES), eO2(NO_BOXES),      &
      !       &       rumc(NO_BOXES),  &
      !       &       rugc(NO_BOXES), rugn(NO_BOXES), rugp(NO_BOXES),  &
      !       &       runc(NO_BOXES), runn(NO_BOXES), runp(NO_BOXES),  &
      !       &       rrsc(NO_BOXES), rrac(NO_BOXES), reac(NO_BOXES),  &
      !       &       rdc(NO_BOXES) , rrtc(NO_BOXES), ruPBAc(NO_BOXES), ruPPYc(NO_BOXES), &
      !       &       ruMIZc(NO_BOXES), rric(NO_BOXES), rr1c(NO_BOXES), rr6c(NO_BOXES), &
      !       &       rr1p(NO_BOXES), rr1n(NO_BOXES), zooc(NO_BOXES), rrip(NO_BOXES), &
      !       &       rr6p(NO_BOXES), rep(NO_BOXES), rrin(NO_BOXES), rr6n(NO_BOXES), &
      !       &       ren(NO_BOXES), pu_ra(NO_BOXES), r(NO_BOXES), &
      !       &       tfluxC(NO_BOXES), tfluxN(NO_BOXES), tfluxP(NO_BOXES), &
      !       &       pe_N4n(NO_BOXES), pe_N1p(NO_BOXES), pe_R6c(NO_BOXES), MIZ_FLUCT(NO_BOXES), &
      !       &      STAT = AllocStatus )
      
      !    IF( AllocStatus /= 0 ) call bfm_error('MicroZooDynamics','Error allocating arrays')
      !    first=1
      ! end if
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !  Copy  state var. object in local var
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Retrieve local biomass (carbon, phosphorus, nitrogen, chlorophyll
      ! concentrations).
      
      ! Concentrations excluding background (used in sink terms)
      _GET_(self%id_c,zooc)
      _GET_(self%id_n,zoon)
      _GET_(self%id_p,zoop)
      
      ! Retrieve ambient nutrient concentrations
      _GET_(self%id_O2o,O2o)
      
      ! Retrieve environmental dependencies (water temperature)
      _GET_(self%id_ETW,ETW)
      
      ! Quota collectors
      qncMIZ = zoon/(zooc+p_small) ! add some epsilon (add in shared) to avoid divide by 0
      qpcMIZ = zoop/(zooc+p_small) ! add some epsilon (add in shared) to avoid divide by 0

      ! Get prey concentrations and quotas
      do iprey = 1, self%nprey
        _GET_(self%id_preyc(iprey), preycP(iprey))
        _GET_(self%id_preyn(iprey), preynP(iprey))
        _GET_(self%id_preyp(iprey), preypP(iprey))
        _GET_(self%id_preyl(iprey), preylP(iprey))
        _GET_(self%id_preys(iprey), preysP(iprey))
!#ifdef INCLUDE_PELFE
      ! _GET_(self%id_preyf(iprey), preyfP(iprey))
!#endif

        ! Quota collectors
        qpcPPY(iprey) = preypP(iprey)/(preycP(iprey)+p_small) ! add some epsilon (add in shared) to avoid divide by 0
        qncPPY(iprey) = preynP(iprey)/(preycP(iprey)+p_small) ! add some epsilon (add in shared) to avoid divide by 0
        qlcPPY(iprey) = preylP(iprey)/(preycP(iprey)+p_small) ! add some epsilon (add in shared) to avoid divide by 0
        qscPPY(iprey) = preysP(iprey)/(preycP(iprey)+p_small) ! add some epsilon (add in shared) to avoid divide by 0
      enddo
      ! Prey carbon was returned in mmol (due to units of standard_variables%total_carbon); convert to mg
      ! MAYBE NOT NECESSARY preycP = preycP*CMass
      
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Temperature effect
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      et = eTq(ETW, self%p_q10)
      
      _SET_DIAGNOSTIC_(self%id_ETWd,ETW)
      _SET_DIAGNOSTIC_(self%id_et,et)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Oxygen limitation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      eO2 = min(ONE, MM(O2o, self%p_chro))
      
      _SET_DIAGNOSTIC_(self%id_eO2,eO2)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate total potential food given the non-dim prey availability
      ! and capture efficiency with loops over all LFGs.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rumc   = ZERO
      do iprey = 1, self%nprey
        PPYc(iprey) = self%p_pa(iprey)*preycP(iprey)* &
        MM(preycP(iprey), self%p_minfood)
        rumc = rumc + PPYc(iprey)
      end do
      
      _SET_DIAGNOSTIC_(self%id_rumc,rumc)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate total food uptake rate (eq 38 Vichi et al. 2007) and 
      ! specific uptake rate considering potentially available food (sut)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rugc = et*self%p_sum*MM(rumc, self%p_chuc)*zooc
      sut = rugc / (p_small + rumc)
      
      _SET_DIAGNOSTIC_(self%id_rugc,rugc)
      _SET_DIAGNOSTIC_(self%id_sut,sut)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Total Gross Uptakes from every LFG
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Bacterioplankton
      rugn = ZERO
      rugp = ZERO
      
      do iprey = 1, self%nprey
        ruPPYc = sut*PPYc(iprey)
        ! All the predation flux from the prey are assigned in the istate loop below
        ! call quota_flux(iiPel, ppzooc, ppPelBacteria(i,iiC), ppzooc, ruPBAc, tfluxC)
        _SET_ODE_(self%id_c,            ruPPYc)
        _SET_ODE_(self%id_preyc(iprey),-ruPPYc)
        ! call quota_flux(iiPel, ppzoon, ppPelBacteria(i,iiN), ppzoon, ruPBAc*qncPBA(i,:), tfluxN)
        _SET_ODE_(self%id_n,            ruPPYc*qncPPY(iprey))
        _SET_ODE_(self%id_preyn(iprey),-ruPPYc*qncPPY(iprey))
        ! call quota_flux(iiPel, ppzoop, ppPelBacteria(i,iiP), ppzoop, ruPBAc*qpcPBA(i,:), tfluxP)
        _SET_ODE_(self%id_p,            ruPPYc*qpcPPY(iprey))
        _SET_ODE_(self%id_preyp(iprey),-ruPPYc*qpcPPY(iprey))
        
        rugn = rugn + ruPPYc*qncPPY(iprey)
        rugp = rugp + ruPPYc*qpcPPY(iprey)
        
        ! Chl is transferred to the infinite sink
        ! call flux_vector(iiPel, ppPhytoPlankton(i,iiL), &
        !                  ppPhytoPlankton(i,iiL), -ruPPYc*qlcPPY(i,:))
        _SET_ODE_(self%id_preyl(iprey),-ruPPYc*qlcPPY(iprey))
        
        ! silicon constituent is transferred to biogenic silicate
        ! if ( ppPhytoPlankton(i,iiS) .gt. 0 ) & 
        !   call flux_vector(iiPel, ppPhytoPlankton(i,iiS), ppR6s, ruPPYc*qscPPY(i,:))
        _SET_ODE_(self%id_R6s,          ruPPYc*qscPPY(iprey))
        _SET_ODE_(self%id_preys(iprey),-ruPPYc*qscPPY(iprey))

        ! #ifdef INCLUDE_PELFE
        !     ! Fe constituent is transferred to particulate iron
        !     if ( ppPhytoPlankton(i,iiF) .gt. 0 ) & 
        !        call flux_vector(iiPel, ppPhytoPlankton(i,iiF), ppR6f, ruPPYc*qfcPPY(i,:))
        ! #endif
!#if defined INCLUDE_PELCO2
!    ! PIC (calcite/aragonite) production associated to the grazed biomass
!    ! The idea in PISCES is that the calcite flux exists only when associated
!    ! to a carbon release from phytoplankton (there is no calcite storage in
!    ! phyto)
!    ! Use the realized rain ratio for each phytoplankton species and assume
!    ! that only a portion is egested
!    ! Calcite production is parameterized as a flux between DIC and PIC
!    ! that affects alkalinity
     if (self%p_isP2(iprey).eq.1) then
!    call flux_vector( iiPel, ppO3c,ppO5c, p_pecaco3(zoo)*ruPPYc*qccPPY(i,:))
        _SET_ODE_(self%id_O3c,-self%p_pecaco3*ruPPYc*qccPPY)        ! precipitation of CaCO3 consumes DIC
        _SET_ODE_(self%id_O5c,self%p_pecaco3*ruPPYc*qccPPY)         ! precipitation of CaCO3 produces PIC
!    call flux_vector( iiPel, ppO3h,ppO3h,-C2ALK*p_pecaco3(zoo)*ruPPYc*qccPPY(i,:))
        _SET_ODE_(self%id_O3h,-C2ALK*self%p_pecaco3*ruPPYc*qccPPY)  ! precipitation of CaCO3 consumes 2 alkalinity

        _SET_DIAGNOSTIC_(self%id_CaCO3precip(iprey), self%p_pecaco3*ruPPYc*qccPPY)
        _SET_DIAGNOSTIC_(self%id_CaCO3_to_O3h(iprey),-C2ALK*self%p_pecaco3*ruPPYc*qccPPY)
     endif
!#endif
      
        
      end do
      
      _SET_DIAGNOSTIC_(self%id_rugn,rugn)
      _SET_DIAGNOSTIC_(self%id_rugp,rugp)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Fluxes from microzooplankton
      ! The metabolic balance is the following:
      ! Ingestion = Growth + Excretion + Respiration
      ! Assimilation efficiency p_pu = G/I
      ! Excretion E = I*p_pu_ea
      ! therefore R = (1-p_pu-p_pu_ea)*I
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Rest, activity and total respiration fluxes
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rrsc = self%p_srs*et*zooc
      ! the activity respiration is derived from the other constant parameters
      rrac = rugc*(ONE - self%p_pu - self%p_pu_ea)
      rrtc = rrsc + rrac
      ! call quota_flux(iiPel, ppzooc, ppzooc, ppO3c, rrtc, tfluxC)
      _SET_ODE_(self%id_c, -rrtc)
      _SET_ODE_(self%id_O3c,rrtc)
      ! call flux_vector(iiPel, ppO2o, ppO2o, -rrtc/MW_C)
      _SET_ODE_(self%id_O2o,-(rrtc/MW_C))
      
      _SET_DIAGNOSTIC_(self%id_rrsc,rrsc)
      _SET_DIAGNOSTIC_(self%id_rrac,rrac)
      _SET_DIAGNOSTIC_(self%id_rrtc,rrtc)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Mortality (rdc) + Activity Excretion (reac)
      ! and partitioning between particulate and dissolved
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rdc  = ((ONE - eO2)*self%p_sdo + self%p_sd)*zooc
      reac = rugc*(ONE - self%p_pu)*self%p_pu_ea
      rric = reac + rdc
      rr1c = rric*self%p_pe_R1c
      rr6c = rric*(ONE - self%p_pe_R1c)

      ! call quota_flux(iiPel, ppzooc, ppzooc, ppR1c, 0.98D0*rr1c, tfluxC) ! flux to non CDOM
      _SET_ODE_(self%id_c, -(1.00D0-self%p_fX1z)*rr1c)
      _SET_ODE_(self%id_R1c,(1.00D0-self%p_fX1z)*rr1c)
      ! call quota_flux(iiPel, ppzooc, ppzooc, ppR1l, 0.02D0*rr1c, tfluxC) ! flux to CDOM
      _SET_ODE_(self%id_c, -self%p_fX1z*rr1c)
      _SET_ODE_(self%id_X1c,self%p_fX1z*rr1c)

      !CEA to DOC: (1-fX1) of reac and all rdc
!      _SET_ODE_(self%id_c, -(((1.00D0-self%p_fX1z)*(reac*rr1c/rric))+(rdc*rr1c/rric))  )
!      _SET_ODE_(self%id_R1c,(((1.00D0-self%p_fX1z)*(reac*rr1c/rric))+(rdc*rr1c/rric))  )
      !CEA to CDOM: fX1 of reac
!      _SET_ODE_(self%id_c, -self%p_fX1z*(reac*rr1c/rric))
!      _SET_ODE_(self%id_X1c,self%p_fX1z*(reac*rr1c/rric))

      ! call quota_flux(iiPel, ppzooc, ppzooc, ppR6c, rr6c, tfluxC)
      _SET_ODE_(self%id_c, -rr6c)
      _SET_ODE_(self%id_R6c,rr6c)
      
      _SET_DIAGNOSTIC_(self%id_rdc,rdc)
      _SET_DIAGNOSTIC_(self%id_reac,reac)
      _SET_DIAGNOSTIC_(self%id_rric,rric)
      _SET_DIAGNOSTIC_(self%id_rr1c,rr1c)
      _SET_DIAGNOSTIC_(self%id_rr6c,rr6c)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !     Nutrient dynamics in microzooplankton
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Organic Nitrogen dynamics
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rrin = rugn*self%p_pu_ea + rdc*qncMIZ
      rr1n = rrin*self%p_pe_R1n
      rr6n = rrin - rr1n
      !   call quota_flux(iiPel, ppzoon, ppzoon, ppR1n, rr1n, tfluxN)
      _SET_ODE_(self%id_n, -rr1n)
      _SET_ODE_(self%id_R1n,rr1n)
      !   call quota_flux(iiPel, ppzoon, ppzoon, ppR6n, rr6n, tfluxN)
      _SET_ODE_(self%id_n, -rr6n)
      _SET_ODE_(self%id_R6n,rr6n)
      
      _SET_DIAGNOSTIC_(self%id_rrin,rrin)
      _SET_DIAGNOSTIC_(self%id_rr1n,rr1n)
      _SET_DIAGNOSTIC_(self%id_rr6n,rr6n)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Organic Phosphorus dynamics
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rrip = rugp*self%p_pu_ea + rdc*qpcMIZ
      rr1p = rrip*self%p_pe_R1p
      rr6p = rrip - rr1p
      ! call quota_flux(iiPel, ppzoop, ppzoop, ppR1p, rr1p, tfluxP)
      _SET_ODE_(self%id_p, -rr1p)
      _SET_ODE_(self%id_R1p,rr1p)
      ! call quota_flux(iiPel, ppzoop, ppzoop, ppR6p, rr6p, tfluxP)
      _SET_ODE_(self%id_p, -rr6p)
      _SET_ODE_(self%id_R6p,rr6p)
      
      _SET_DIAGNOSTIC_(self%id_rrip,rrip)
      _SET_DIAGNOSTIC_(self%id_rr1p,rr1p)
      _SET_DIAGNOSTIC_(self%id_rr6p,rr6p)
      
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved nutrient dynamics
      ! Compare the quota of the net growth rates with the optimal quota
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      runc = max(ZERO, rugc*(ONE - self%p_pu_ea) - rrac)
      runn = max(ZERO, rugn*(ONE - self%p_pu_ea) + rrsc*qncMIZ)
      runp = max(ZERO, rugp*(ONE - self%p_pu_ea) + rrsc*qpcMIZ)
      ren  = max(ZERO, runn/(p_small + runc) - self%p_qncMIZ)* runc
      rep  = max(ZERO, runp/(p_small + runc) - self%p_qpcMIZ)* runc
      ! call quota_flux(iiPel, ppzoon, ppzoon, ppN4n, ren, tfluxN)
      _SET_ODE_(self%id_n, -ren)
      _SET_ODE_(self%id_N4n,ren)
      _SET_ODE_(self%id_O3h, ren)      ! release of 1 NH4 increases 1 alkalinity
      ! call quota_flux(iiPel, ppzoop, ppzoop, ppN1p, rep, tfluxP)
      _SET_ODE_(self%id_p, -rep)
      _SET_ODE_(self%id_N1p,rep)
!      _SET_ODE_(self%id_O3h, -rep)      ! release of 1 PO4 decreases 1 alkalinity


      _SET_DIAGNOSTIC_(self%id_runc,runc)
      _SET_DIAGNOSTIC_(self%id_runn,runn)
      _SET_DIAGNOSTIC_(self%id_runp,runp)
      _SET_DIAGNOSTIC_(self%id_ren,ren)
      _SET_DIAGNOSTIC_(self%id_rep,rep)
     
      _SET_DIAGNOSTIC_(self%id_varO3h_Nutil,ren)
      _SET_DIAGNOSTIC_(self%id_varO3h_Putil,-rep)
 
      _LOOP_END_
    end subroutine do
end module
! SKIP FROM HERE (slow)   if ( ppzoon == 0 .or. ppzoop == 0 ) then
!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      ! Eliminate the excess of the non-limiting constituent under fixed quota
!      ! Determine whether C, P or N is limiting (Total Fluxes Formulation)
!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      limit = nutlim(tfluxc,tfluxn,tfluxp,qncMIZ(zoo,:),qpcMIZ(zoo,:),iiC,iiN,iiP)

!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      ! Compute the correction terms depending on the limiting constituent
!      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      WHERE     ( limit == iiC )
!          pe_N1p = max(ZERO,tfluxp  - p_qpcMIZ(zoo)* tfluxc)
!          pe_N4n = max(ZERO,tfluxn  - p_qncMIZ(zoo)* tfluxc)
!          pe_R6c = ZERO
!      ELSEWHERE ( limit == iiP )
!          pe_N1p = ZERO
!          pe_N4n = max(ZERO, tfluxn  - tfluxp/p_qpcMIZ(zoo)*p_qncMIZ(zoo) )
!          pe_R6c = max(ZERO, tfluxc  - tfluxp/p_qpcMIZ(zoo))
!      ELSEWHERE ( limit == iiN )
!          pe_N1p = max(ZERO, tfluxp  - tfluxn/p_qncMIZ(zoo)*p_qpcMIZ(zoo))
!          pe_N4n = ZERO
!          pe_R6c = max(ZERO, tfluxc  - tfluxn/p_qncMIZ(zoo))
!      END WHERE

!      call flux_vector(iiPel, ppzooc, ppR6c, pe_R6c*(ONE-p_pe_R1c) )
!      call flux_vector(iiPel, ppzooc, ppR1c, 0.98D0*pe_R6c*(p_pe_R1c))
!      call flux_vector(iiPel, ppzooc, ppR1l, 0.02D0*pe_R6c*(p_pe_R1c)) ! To CDOM
! !    call flux_vector(iiPel, ppzooc, ppR1c, pe_R6c*(p_pe_R1c))
!      call flux_vector(iiPel, ppzoop, ppN1p, pe_N1p)
!      call flux_vector(iiPel, ppzoon, ppN4n, pe_N4n)

! #ifdef DEBUG
!      write(*,*) '+++++++++++++++'
!      if ( limit(1)==iiC ) then
!      write(*,*) 'tfluxp', tfluxp,'pe_N1p', tfluxp  - p_qpcMIZ(zoo)* tfluxc 
!      write(*,*) 'tfluxn', tfluxn,'pe_N4n', tfluxn  - p_qncMIZ(zoo)* tfluxc
!      write(*,*) 'tfluxc', tfluxc,'pe_R6c', ZERO 
!      write(*,*) 'ooooooooooooooo'
!      endif

!      if ( limit(1)==iiP ) then
!      write(*,*) 'tfluxp', tfluxp,'pe_N1p', ZERO 
!      write(*,*) 'tfluxn', tfluxn,'pe_N4n', tfluxn - tfluxp/p_qpcMIZ(zoo)*p_qncMIZ(zoo)
!      write(*,*) 'tfluxc', tfluxc,'pe_R6c', tfluxc - tfluxp/p_qpcMIZ(zoo)
!      write(*,*) 'ooooooooooooooo'
!      endif

!      if ( limit(1)==iiN ) then
!      write(*,*) 'tfluxp', tfluxp,'pe_N1p', tfluxp  - tfluxn/p_qncMIZ(zoo)*p_qpcMIZ(zoo)
!      write(*,*) 'tfluxn', tfluxn,'pe_N4n', ZERO
!      write(*,*) 'tfluxc', tfluxc,'pe_R6c', tfluxc  - tfluxn/p_qncMIZ(zoo) 
!      endif
!      write(*,*) '+++++++++++++++'
! #endif
! TO HERE
!   endif

!   end subroutine MicroZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
