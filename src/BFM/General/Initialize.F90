!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: Initialize
!
! DESCRIPTION
!   Initialization of model, with allocation of memory and 
!   reading of data files 
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
!
! INTERFACE
  SUBROUTINE Initialize
!
! USES
  use mem, only: InitializeModel
  use mem_Param
  use mem_PelChem
  use mem_PelBac
  use mem_MesoZoo
  use mem_MicroZoo
  use mem_Phyto
  use mem_PAR
  use mem_PelSinkSet
#if defined BENTHIC_BIO || defined BENTHIC_FULL 
  use mem_BenBac
  use mem_BenOrganisms
  use mem_FilterFeeder
  use mem_Bioturbation
  use mem_BenthicRemin
  use mem_BenOxygen
  use mem_ControlBennutBuffers
#if defined BENTHIC_FULL
  use mem_BenthicNutrient
  use mem_BenAmmonium
  use mem_BenNitrate
  use mem_BenAnoxic
  use mem_BenDenitriDepth
  use mem_BenPhosphate
  use mem_BenSilica
  use mem_BenQ1Transport
#endif
#else
  use mem_BenthicReturn
#endif
#ifdef INCLUDE_PELCO2
  use mem_CO2
#endif
#ifdef INCLUDE_BENCO2
  use mem_BenCO2Transport
  use mem_BenAlkalinity
#endif
#ifdef INCLUDE_SEAICE
  use mem_SeaiceAlgae
  use mem_SeaiceBac
  use mem_SeaiceZoo
  use mem_SeaicetoPel
#endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    InitializeModel=0

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Allocate Memory for All global variables
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    call AllocateMem

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Initialize state variables transport and bottom box porosity
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    call InitTransportStateTypes
    call InitBoxParams

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Read all data files:(namelist files)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    call InitParam
    call InitPelChem
    call InitPelBac
    call InitMesoZoo
    call InitMicroZoo
    call InitPhyto
    call InitPAR
    call InitPelSinkSet

#if defined BENTHIC_BIO
    ! Intermediate benthic return
    call InitBenOrganisms
    call InitFilterFeeder
    call InitBenBac
    call InitBioturbation
    call InitBenthicRemin
    call InitBenOxygen
    call InitControlBennutBuffers
#elif defined BENTHIC_FULL
    ! Full benthic nutrients
    call InitBenOrganisms
    call InitFilterFeeder
    call InitBenBac
    call InitBioturbation
    call InitBenthicNutrient
    call InitBenAmmonium
    call InitBenNitrate
    call InitBenOxygen
    call InitBenAnoxic
    call InitBenDenitriDepth
    call InitBenPhosphate
    call InitBenSilica
    call InitBenQ1Transport
    call InitControlBennutBuffers
#else
    ! Simple benthic return
    call InitBenthicReturn
#endif
#ifdef INCLUDE_BENCO2
    call InitBenCO2Transport
    call InitBenAlkalinity
#endif
#ifdef INCLUDE_PELCO2
    call InitCO2
#endif
#ifdef INCLUDE_SEAICE
    call InitSeaiceAlgae
    call InitSeaiceBac
    call InitSeaiceZoo
    call InitSeaicetoPel
#endif

  END SUBROUTINE Initialize

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
