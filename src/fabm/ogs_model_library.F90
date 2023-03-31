module ogs_model_library
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

   use fabm_types, only: type_base_model_factory, type_base_model

   use ogs_bfm_shared
   use ogs_bfm_pelagic_base
   use bfm_Phyto 
   use bfm_PelBac
   use bfm_PelChem
   use bfm_PelOxygen
   use bfm_MicroZoo
   use bfm_MesoZoo
   use ogs_bfm_light
!  use ogs_bfm_light_spectral
   use bfm_zenith_angle
   use bfm_CalciteDissolution
   use bfm_PelagicCSYS


   ! Add use statements for new models here

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
!     procedure :: initialize
      procedure :: create
   end type

   type (type_factory), save, target, public :: ogs_model_factory

contains

!     subroutine initialize(self)
!     class (type_factory), intent(inout) :: self
!     call self%register_version('ERSEM',git_commit_id//' ('//git_branch_name//' branch)')
!      end subroutine initialize

   subroutine create(self, name, model)

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         ! Add case statements for new models here
         case ('bfm_pelagic_base'); allocate(type_ogs_bfm_pelagic_base::model)
         case ('Phyto'); allocate(type_ogs_bfm_primary_producer::model)
         case ('PelBac'); allocate(type_ogs_bfm_pelagic_bacteria::model)
         case ('PelChem'); allocate(type_ogs_bfm_PelChem::model)
         case ('PelOxygen'); allocate(type_ogs_bfm_PelOxygen::model)
         case ('MicroZoo'); allocate(type_ogs_bfm_microzoo::model)
         case ('MesoZoo'); allocate(type_ogs_bfm_mesozoo::model)
         case ('light'); allocate(type_ogs_bfm_light::model)
!        case ('light_spectral'); allocate(type_ogs_bfm_light_spectral::model)
         case ('zenith_angle'); allocate(type_ogs_bfm_zenith_angle::model)
         case ('CalciteDissolution'); allocate(type_ogs_bfm_CalciteDissolution::model)
         case ('PelagicCSYS'); allocate(type_ogs_bfm_PelagicCSYS::model)
      end select

   end subroutine create

end module
