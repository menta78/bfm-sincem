#include "cppdefs.h"
#include "DEBUG.h"
#include "INCLUDE.h"

#ifdef INCLUDE_PELCO2
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelCO2Dynamics
!
! DESCRIPTION
!   !
!
! !INTERFACE
  subroutine PelagicCSYS()
!
! !USES:

  use global_mem, ONLY: RLEN,ONE,ZERO
  use constants,  ONLY: MW_C, C2ALK
  use mem_Param,  ONLY: AssignAirPelFluxesInBFMFlag
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: iiPel, O3h, O3c, D3STATE, jsurO3c, CO2airflux,    &
                 Depth, flux_vector, DIC, EPCO2air, ALK,           &
                 Source_D3_vector, ppO5c, ppN3n, ppN4n
  use mem, ONLY: ppO3h, ppO3c, NO_BOXES, NO_BOXES_XY, BoxNumber,   &
    N1p,N5s,CO2, HCO3, CO3, pCO2, pH, ETW, ESW, ERHO, EWIND, EICE, &
    OCalc, OArag, EPR, ppO5c, O5c
#endif
  use mem_CO2    
  use mem_CSYS, ONLY : CarbonateSystem
  use AirSeaExchange, ONLY: AirSeaCO2
  use BFM_ERROR_MSG, ONLY: BFM_ERROR
#ifdef BFM_GOTM
  use bio_var, ONLY: SRFindices
#else
  use api_bfm, ONLY: SRFindices,bfm_init
#endif
  IMPLICIT NONE

!  
!
! !AUTHORS
!  T. Lovato 
!
! !REVISION_HISTORY
!
! !LOCAL VARIABLES:
  integer            :: error=0
  integer,save       :: first=0
  integer            :: AllocStatus
  real(RLEN),allocatable,save,dimension(:) :: rateN, excess, rdiss, xflux
!
! COPYING
!   
!   Copyright (C) 2017 BFM System Team (bfm_st@lists.cmcc.it)
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
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     ALLOCATE ( rateN(NO_BOXES), excess(NO_BOXES), rdiss(NO_BOXES),  &
        &      xflux(NO_BOXES_XY),                                   &
        &      STAT = AllocStatus )
     IF( AllocStatus /= 0 ) call bfm_error('PelagicCSYS','Error allocating arrays')
     first=1
  end if
  if (bfm_init == 0 ) pH(:) = -ONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute carbonate system equilibria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! To use the Pressure correction of CSYS here the pr_in=EPS value
    do BoxNumber=1,NO_BOXES
       ! convert DIC and alkalinity from model units to diagnostic output
       ! mg C/m3 --> umol/kg
       ! mmol eq/m3 --> umol/kg
       DIC(BoxNumber) = O3c(BoxNumber)/MW_C/ERHO(BoxNumber)*1000.0_RLEN
       ALK(BoxNumber) = O3h(BoxNumber)/ERHO(BoxNumber)*1000.0_RLEN

       error= CarbonateSystem( ESW(BoxNumber), ETW(BoxNumber),ERHO(BoxNumber), &
               N1p(BoxNumber), N5s(BoxNumber), DIC(BoxNumber), ALK(BoxNumber), &
               CO2(BoxNumber) ,HCO3(BoxNumber), CO3(BoxNumber), pH(BoxNumber), &
               pCO2(BoxNumber), patm=AtmSLP%fnow(BoxNumber), pr_in=EPR(BoxNumber), &
               OmegaC=OCalc(BoxNumber), OmegaA=OArag(BoxNumber))
#ifdef DEBUG
       write(LOGUNIT,*) "in PelagicCSYS:"
       write(LOGUNIT,'(A,'' ='',f12.6)') 'ERHO',ERHO(BoxNumber)
       write(LOGUNIT,'(A,'' ='',f12.6)') 'ESW',ESW(BoxNumber)
       write(LOGUNIT,'(A,'' ='',f12.6)') 'N1p',N1p(BoxNumber)
       write(LOGUNIT,'(A,'' ='',f12.6)') 'N5s',N5s(BoxNumber)
       write(LOGUNIT,'(A,'' ='',f12.6)') 'DIC',DIC(BoxNumber)
       write(LOGUNIT,'(A,'' ='',f12.6)') 'ALK',ALK(BoxNumber)
       write(LOGUNIT,'(A,'' ='',f12.6)') 'OCalc',OCalc(BoxNumber)
       write(LOGUNIT,'(A,'' ='',f12.6)') 'OArag',OArag(BoxNumber)
       write(LOGUNIT,'(''layer:'',I4,'' pH='',f12.6)') BoxNumber,pH(BoxNumber)
#endif
       if ( error > 0 ) then
              call BFM_ERROR("PelagicCSYS","Error in csys computation")
              write(LOGUNIT,'(''layer:'',I4,'' pH='',f12.6)') BoxNumber,pH(BoxNumber)
       endif
    end do
   
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Computes Atmospheric pCO2 value
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Rough approximation: pCO2 is assumed equal to the mixing ratio of CO2
  EPCO2air = AtmCO2%fnow
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Computes air-sea flux (only at surface points)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  CO2airflux = AirSeaCO2(AtmCO2%fnow, AtmSLP%fnow, ETW(SRFindices), ESW(SRFindices), &
          ERHO(SRFindices), EWIND, EICE, CO2(SRFindices) )
  !
  jsurO3c(:) = jsurO3c(:) + CO2airflux * MW_C
  xflux = jsurO3c(:) / Depth(SRFindices) * CO2fluxfac
  if ( AssignAirPelFluxesInBFMFlag)  call flux_vector( iiPel, ppO3c,ppO3c, xflux )
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Changes in alkalinity due to N uptake (see BFM Manual Eq. 2.5.21)
  ! It is computed this way
  ! net_uptakeNO3=dNO3/dt+denit-nit , net_uptakeNH4=dNH4/dt+nit
  ! dTA/dt = -net_uptakeNO3 + net_uptakeNH4 - 2*nit + denit = -dNO3/dt+dNH4/dt
  ! Sulfur reactions associated to reduction equivalents are not
  ! considered as included in the operational TA definition
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if ( CalcBioAlk ) then
     rateN(:) = - Source_D3_vector(ppN3n) + Source_D3_vector(ppN4n)
     call flux_vector( iiPel, ppO3h,ppO3h, rateN)
  endif 

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Dissolution of Particulate Inorganic Carbon (calcite/aragonite) in seawater
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute undersaturation
  excess(:) = max(ZERO,ONE - OCalc(:))
  ! Dissolution rate of C in CaCO3 (mg C/m3/d) from Morse and Berner (1972)
  rdiss(:) = p_kdca * excess(:)**p_nomega * O5c(:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Inorganic carbon and alkalinity flux due to PIC changes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call flux_vector( iiPel, ppO5c, ppO3c, rdiss(:) )
  call flux_vector( iiPel, ppO3h, ppO3h, -C2ALK*rdiss(:) )

  end subroutine PelagicCSYS
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
