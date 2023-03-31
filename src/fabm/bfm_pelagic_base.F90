#include "fabm_driver.h"

module ogs_bfm_pelagic_base
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
   use fabm_particle

   use ogs_bfm_shared

   implicit none

   private

   ! type, extends(type_base_model), public :: type_ogs_bfm_pelagic_base
   type,extends(type_particle_model),public :: type_ogs_bfm_pelagic_base
      type (type_state_variable_id)                 :: id_c,id_n,id_p,id_f,id_s,id_chl,id_o,id_r,id_h
      ! Add variable identifiers and parameters here.
      type (type_horizontal_dependency_id)          :: id_bedstress,id_wdepth
      type (type_dependency_id)                     :: id_dens
      type (type_horizontal_diagnostic_variable_id) :: id_w_bot
      type (type_horizontal_diagnostic_variable_id),allocatable,dimension(:) :: id_cdep,id_ndep,id_pdep,id_sdep,id_fdep,id_odep

      ! Target variables for sedimentation
      type (type_bottom_state_variable_id),allocatable,dimension(:) :: id_targetc,id_targetn,id_targetp,id_targets,id_targetf,id_targeto

      real(rk) :: rm = 0.0_rk
      real(rk) :: tdep
      integer :: ndeposition
      logical :: no_river_dilution = .false.

      real(rk),allocatable :: qxc(:),qxn(:),qxp(:),qxs(:),qxf(:),qxo(:)

   contains
      procedure :: initialize
      procedure :: initialize_bfm_base
      procedure :: add_constituent
      ! Reference model procedures here.
   end type

contains

   subroutine initialize(self, configunit)
      class (type_ogs_bfm_pelagic_base), intent(inout), target :: self
      integer,                           intent(in)            :: configunit

      character(len=10) :: composition
      real(rk)          :: c0,s0,rRPmX,EPS,iopABS,iopBBS
      real(rk)          :: qn, qp


      call self%get_parameter(composition, 'composition', '', 'elemental composition')
!     call self%get_parameter(c0, 'c0', 'mg C/m^3', 'background concentration in carbon units', default=0.0_rk, minimum=0.0_rk)

      if (index(composition,'c')/=0) then
         ! Carbon-based light attenuation [optional, off by default]
         call self%get_parameter(EPS,   'EPS',    'm^2/mg C', 'specific shortwave attenuation', default=0.0_rk, minimum=0.0_rk)
         call self%get_parameter(iopABS,'iopABS', 'm^2/mg C', 'specific shortwave absorption',  default=0.0_rk, minimum=0.0_rk)
         call self%get_parameter(iopBBS,'iopBBS', 'm^2/mg C', 'specific shortwave backscatter', default=0.0_rk, minimum=0.0_rk)

         ! Constant N:C stoichiometry [optional, off by default, and only
         ! supported if no explicit nitrogen constituent is present]
         if (index(composition,'n') == 0) then
            call self%get_parameter(qn, 'qn', 'mmol N/mg C', 'nitrogen : carbon ratio', default=0.0_rk, minimum=0.0_rk)
         else
            qn = 0.0_rk
         end if

         ! Constant P:C stoichiometry [optional, off by default, and only
         ! supported if no explicit phosphorus constituent is present]
         if (index(composition,'p') == 0) then
            call self%get_parameter(qp, 'qp', 'mmol P/mg C', 'phosphorus : carbon ratio', default=0.0_rk, minimum=0.0_rk)
         else
            qp = 0.0_rk
         end if
      end if
      call self%get_parameter(self%rm, 'rm', 'm/d', 'sinking velocity', default=0.0_rk)

!     call self%initialize_ersem_base(rm=rRPmX, sedimentation=rRPmX>0._rk)

      if (index(composition,'c')/=0) then
         call self%add_constituent('c', 0.0_rk, c0, qn, qp)

         ! Add contributions to light attenuation, absorption, scattering.
         ! Contributions with a scale_factor of 0.0 will automatically be
         ! ignored.
         call self%add_to_aggregate_variable(particulate_organic_absorption_coefficient, &
            self%id_c,scale_factor=iopABS,include_background=.true.)
         call self%add_to_aggregate_variable(particulate_organic_backscatter_coefficient, &
            self%id_c,scale_factor=iopBBS,include_background=.true.)
         call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
            self%id_c,scale_factor=EPS,include_background=.true.)
      end if
      if (index(composition,'n')/=0) call self%add_constituent('n',0.0_rk)
      if (index(composition,'p')/=0) call self%add_constituent('p',0.0_rk)
      if (index(composition,'s')/=0) then 
         call self%add_constituent('s',0.0_rk,s0)
      end if
      if (index(composition,'f')/=0) call self%add_constituent('f',0.0_rk)
      if (index(composition,'o')/=0) call self%add_constituent('o',0.0_rk)
      if (index(composition,'r')/=0) call self%add_constituent('r',0.0_rk)
      if (index(composition,'h')/=0) call self%add_constituent('h',0.0_rk)
 
      ! Register model parameters and variables here.

   end subroutine initialize

   subroutine initialize_bfm_base(self)
     class (type_ogs_bfm_pelagic_base), intent(inout), target :: self
     

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are
      ! given in d-1.
      self%dt = 86400._rk

   end subroutine initialize_bfm_base

   ! Add model subroutines here.
      subroutine add_constituent(self,name,initial_value,background_value,qn,qp)
      class (type_ogs_bfm_pelagic_base), intent(inout), target :: self
      character(len=*),                intent(in)            :: name
      real(rk),                        intent(in)            :: initial_value
      real(rk),optional,               intent(in)            :: background_value,qn,qp


      select case (name)
      case ('c')
         call register(self%id_c,'c','mg C','carbon',standard_variables%total_carbon,self%qxc,self%id_cdep,self%id_targetc,1._rk/12.011_rk)
         if (present(qn)) call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c,scale_factor=qn)
         if (present(qp)) call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_c,scale_factor=qp)
      case ('n')
         call register(self%id_n,'n','mmol N','nitrogen',standard_variables%total_nitrogen,self%qxn,self%id_ndep,self%id_targetn)
      case ('p')
         call register(self%id_p,'p','mmol P','phosphorus',standard_variables%total_phosphorus,self%qxp,self%id_pdep,self%id_targetp)
      case ('s')
         call register(self%id_s,'s','mmol Si','silicate',standard_variables%total_silicate,self%qxs,self%id_sdep,self%id_targets)
      case ('f')
!        if (use_iron) call register(self%id_f,'f','umol Fe','iron',standard_variables%total_iron,self%qxf,self%id_fdep,self%id_targetf)
      case ('o')
         call register(self%id_o,'o','mmol O2','Oxygen',total_oxygen)
!        call register(self%id_o,'o','mmol O2','Oxygen',standard_variables%total_oxygen,self%qxo,self%id_odep,self%id_targeto)
      case ('chl')
         call register(self%id_chl,'Chl','mg','chlorophyll a',total_chlorophyll)
      case ('r')
         call register(self%id_r,'r','mmol/Eq','mmol Eq',total_reduction_equivalent)
      case ('h')
         call register(self%id_h,'h','mmol eq','mmol Eq alkalinity',alkalinity)
      case default
         call self%fatal_error('add_constituent','Unknown constituent "'//trim(name)//'".')
      end select

   contains

      subroutine register(variable_id,name,base_units,long_name,aggregate_variable,qx,id_xdep,id_dep,scale_factor)
         type (type_state_variable_id),       intent(inout), target     :: variable_id
         character(len=*),                    intent(in)                :: name,base_units,long_name
         type (type_bulk_standard_variable),  intent(in)                :: aggregate_variable
         real(rk),                            intent(inout),allocatable,optional :: qx(:)
         type (type_bottom_state_variable_id),intent(inout),allocatable,optional :: id_dep(:)
         type (type_horizontal_diagnostic_variable_id),intent(inout),allocatable,optional :: id_xdep(:)
         real(rk),                            intent(in), optional      :: scale_factor

         integer :: idep
         character(len=16) :: num

         ! Register state variable
         self%dt = 86400._rk
         call self%register_state_variable(variable_id,name,trim(base_units)//'/m^3',long_name, &
            initial_value,minimum=0._rk,vertical_movement=-self%rm/self%dt,background_value=background_value,no_river_dilution=self%no_river_dilution)
!        write(*,*) "vertical_movement", -self%rm/self%dt

         ! Contribute to aggregate variable
         call self%add_to_aggregate_variable(aggregate_variable,variable_id,scale_factor)

         self%ndeposition=-1.
         if (self%ndeposition>0.and.present(qx)) then
            ! Create array with fractions per depostion target
            allocate(qx(self%ndeposition))
            do idep = 2,self%ndeposition
               write(num,'(i0)') idep
               call self%get_parameter(qx(idep),'qx'//trim(name)//trim(num),'-','fraction of '//trim(long_name)//' sinking into deposition target '//trim(num))
            end do
            qx(1)=1._rk-sum(qx(2:))

            ! Link to pools in which to deposit matter.
            allocate(id_dep(self%ndeposition))
            allocate(id_xdep(self%ndeposition))

            do idep=1,self%ndeposition
               write(num,'(i0)') idep
               if (qx(idep)/=0) then
                  call self%register_state_dependency(id_dep(idep),'deposition_target'//trim(num)//trim(name),trim(base_units)//'/m^2','target pool '//trim(num)//' for '//trim(long_name)//' deposition')
!ERROR-->                 call self%request_coupling_to_model(id_dep(idep),'deposition_target'//trim(num),name)
                  call self%register_diagnostic_variable(id_xdep(idep),'dep'//trim(num)//trim(name),trim(base_units)//'/m^2/d','deposition of '//trim(long_name)//' in target pool '//trim(num),source=source_do_bottom)
               end if
            end do
         end if
      end subroutine register
   end subroutine add_constituent

end module
