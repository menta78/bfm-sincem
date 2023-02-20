!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PelOxygen
!
! DESCRIPTION
!   Compute Oxygen dynamics in the Pelagic environment
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
#include "cppdefs.h"
!
! INTERFACE
  subroutine PelOxygen ()
!
! USES
  use global_mem, ONLY: RLEN,ONE,ZERO
  use constants,  ONLY: HOURS_PER_DAY,ZERO_KELVIN
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: NO_BOXES_XY, NO_BOXES, iiPel, flux_vector
  use mem, ONLY: ppO2o, O2o, jsurO2o, cxoO2, eO2mO2, ETW, ESW, EWIND, EICE, Depth, dpo2
# endif
  use api_bfm,    ONLY: SRFindices
  use mem_param,  ONLY: AssignAirPelFluxesInBFMFlag, CalcPelChemistry, p_small
  use mem_CO2,    ONLY: patm3d
  use AirSeaExchange, ONLY: AirpGas

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)     :: temp, salt, tmpflux, ts, ts2, lnk
  real(RLEN),dimension(NO_BOXES_XY)  :: wind, fice, pschmidt, Kw660, kwgas,   &
                                        tsrf, tsrf2, o2flux, po2air, po2sea
  real(RLEN),dimension(NO_BOXES)  :: o2sat, o2rel, jo2sur
  real(RLEN),parameter    :: kwfac  = 0.01_RLEN * HOURS_PER_DAY ! convert from cm/h to m/d
  real(RLEN),parameter    :: xo2air = 0.20946_RLEN              ! Volume fraction of oxygen in dry air [atm]
  real(RLEN),parameter    :: o2mvol = 22.391_RLEN               ! volume of 1 mole O2 at STP (ICES conversions)
  real(RLEN),parameter    :: o2fac  = 1000._RLEN / o2mvol       ! convert from ml/l to umol/L [=mmol/m3]
  ! O2 solubility parameters [cm3/dm3 = ml/L]
  real(RLEN),parameter :: &
  A(6) = (/ 2.00856, 3.22400,  3.99063, 4.80299, 9.78188E-01, 1.71069 /),  &
  B(4) = (/ -6.24097E-03, -6.93498E-03, -6.90358E-03, -4.29155E-03 /),     &
  C0 = -3.11680e-7
  !! Schmidt coefficients
  real(RLEN),parameter :: Sc = 660.0_RLEN, &
  C(5)  = (/1920.4_RLEN, -135.6_RLEN, 5.2122_RLEN, -0.10939_RLEN, 0.00093777_RLEN/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  ! input arrays
  temp = ETW(:)
  salt = ESW(:)
  wind = EWIND(:)
  fice = EICE(:)

  ! common arrays
  tmpflux(:) = ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Oxygen saturation (Gordon and Garcia, 1992 L&O)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Scaled temperature 
  ts = log ( ( 298.15_RLEN - temp) / ( temp - ZERO_KELVIN ) )
  ts2 = ts * ts
  lnk = A(1) + A(2)*ts + A(3)*ts2 + A(4)*ts*ts2 + A(5)*ts2*ts2 + A(6)*ts2*ts2*ts + &
        salt * (B(1) + B(2)*ts + B(3)*ts2 + B(4)*ts*ts2) + C0*salt*salt
  ! Compute O2sat in (mmol/m3) 
  cxoO2(:)  = exp(lnk) * o2fac
  ! Relative oxygen saturation
  eO2mO2(:) = max( p_small,O2o(:) ) / cxoO2(:)
  
  if ( (SRFindices(1) .EQ. 0 ) .or. ( .NOT. CalcPelChemistry ) )  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Air-Sea Flux of Oxygen
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compute piston velolicty kw660 (at 25 C) from wind speed (Wanninkhof 2014, Limnol. Oceanograph. Methods, 12, 351-362)
  kw660 = 0.251_RLEN * wind**2 * kwfac
  !  Schmidt number (Sc), ratio between the kinematic viscosity and the molecular diffusivity of the gas
  tsrf = temp(SRFindices)
  tsrf2 = tsrf * tsrf 
  pschmidt = C(1) + C(2)*tsrf + C(3)*tsrf2 + C(4)*tsrf*tsrf2 + C(5)*tsrf2*tsrf2
  ! Transfer velocity for gas in m/d (see equation [4] in OCMIP2 design document & OCMIP2 Abiotic HOWTO)
  kwgas = kw660 * (Sc/pschmidt)**0.5
  !
  o2flux = kwgas * ( cxoO2(SRFindices) - O2o(SRFindices) ) * (ONE - Fice)
  !
  ! Update flux in general array
  jsurO2o(:)  = jsurO2o(:) + o2flux
  ! Convert flux to mmol/m2/day
  tmpflux(SRFindices) = jsurO2o(:) / Depth(SRFindices)
  if ( AssignAirPelFluxesInBFMFlag) call flux_vector( iiPel, ppO2o, ppO2o, tmpflux )
  !
  ! Atmospheric partial pressure [atm]
  po2air = AirpGas(xo2air, patm3d(SRFindices), temp(SRFindices), salt(SRFindices)) 
  ! Sea surface partial pressure [atm] (rearranged from Orr et al, 2017)
  po2sea = O2o(SRFindices) * ( po2air / cxoO2(SRFindices) )
  ! difference [atm]
  dpo2(:) = po2air - po2sea

  return

  end subroutine PelOxygen

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
