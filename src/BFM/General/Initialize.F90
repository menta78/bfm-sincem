!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Initialize
!
! DESCRIPTION
!   Initialization of model
!   Allocation of memory for variables, reading of data files 
!
! !INTERFACE
  SUBROUTINE Initialize
!
! USES:
  use mem, only: InitializeModel
  use mem_Param
  use mem_PelGlobal
  use mem_PelChem
  use mem_PelBac
  use mem_MesoZoo
  use mem_MicroZoo
  use mem_Phyto
  use mem_PAR
  use mem_Settling
#ifdef BENTHIC_RETURN
    use mem_BenthicReturn
#endif
#if defined BENTHIC_BIO || defined BENITHIC_FULL 
  use mem_BenBac
  use mem_BenOrganism
  use mem_FilterFeeder
  use mem_Bioturbation
  use mem_BenthicReturn2
  use mem_BenOxygen
  use mem_ControlBennutBuffers
#if defined BENITHIC_FULL
  use mem_BenthicNutrient3
  use mem_BenAmmonium
  use mem_BenNitrate
  use mem_BenAnoxic
  use mem_BenDenitriDepth
  use mem_BenPhosphate
  use mem_BenSilica
  use mem_BenQ1Transport
#endif
#endif
#ifdef INCLUDE_PELCO2
  use mem_CO2
#endif
#ifdef INCLUDE_BENCO2
  use mem_BenCO2Transport
  use mem_BenAlkalinity
#endif
#ifdef INCLUDE_SILT
  use mem_Silt
#endif
#ifdef INCLUDE_SEAICE
  use mem_SeaiceAlgae
  use mem_SeaiceBac
  use mem_SeaiceZoo
  use mem_SeaicetoPel
#endif

!  
!
! !AUTHORS
!   mfstep ERSEM team
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
      call InitSettling
      call InitPelGlobal

      ! Benthic initialization is done only if there is an active model
      ! When INCLUDE_BEN is defined, 
      ! CalcBenthicFlag=0 is used to test the benthic memory only
#ifdef BENTHIC_RETURN
      ! Simple benthic return
      call InitBenthicReturn
#endif
#ifdef BENTHIC_BIO
      ! Intermediate benthic return
      call InitBenOrganism
      call InitFilterFeeder
      call InitBenBac
      call InitBioturbation
      call InitBenthicReturn2
      call InitBenOxygen
      call InitControlBennutBuffers
#endif
#ifdef BENTHIC_FULL
      ! Full benthic nutrients
      call InitBenOrganism
      call InitFilterFeeder
      call InitBenBac
      call InitBioturbation
      call InitBenthicReturn1
      call InitBenthicReturn2
      call InitBenthicNutrient3
      call InitBenAmmonium
      call InitBenNitrate
      call InitBenOxygen
      call InitBenAnoxic
      call InitBenDenitriDepth
      call InitBenPhosphate
      call InitBenSilica
      call InitBenQ1Transport
      call InitControlBennutBuffers
#ifdef INCLUDE_BENCO2
      call InitBenCO2Transport
      call InitBenAlkalinity
#endif
#endif

#ifdef INCLUDE_PELCO2
      call InitCO2
#endif
#ifdef INCLUDE_SILT
      call InitSilt
#endif
#ifdef INCLUDE_SEAICE
      call InitSeaiceAlgae
      call InitSeaiceBac
      call InitSeaiceZoo
      call InitSeaicetoPel
#endif

    END SUBROUTINE Initialize
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
