!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: Air-Sea flux Exchange
!
! DESCRIPTION
!   Compute Air-Sea CO2 exchanges
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
  module AirSeaExchange
!
! USES
  use global_mem, ONLY: RLEN,ONE,ZERO
  use constants,  ONLY: HOURS_PER_DAY,ZERO_KELVIN

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public AirSeaCO2, AirpGas

 contains

  elemental function AirSeaCO2(xco2,patm,temp,salt,rho,wind,Fice,co2ocn)
  ! 
    implicit none
    real(RLEN)                       :: AirSeaCO2 
    !     
    real(RLEN),intent(IN)            :: xco2   ! Atmospheric mixing ratio of Gas [ppmv]
    real(RLEN),intent(IN)            :: patm   ! Atmospheric pressure [mbar=hPa]
    real(RLEN),intent(IN)            :: temp   ! in-situ temperature at surface
    real(RLEN),intent(IN)            :: salt   ! practical salinity at surface
    real(RLEN),intent(IN)            :: rho    ! in-situ density 
    real(RLEN),intent(IN)            :: wind   ! wind speed [m/s]
    real(RLEN),intent(IN)            :: Fice   ! Sea-ice cover fraction [0-1]
    real(RLEN),intent(IN)            :: co2ocn ! gas concentration at surface [umol/kg]
    !
    real(RLEN)           :: xflux  ! Air-Sea Gas Flux (positive downward) [mol/m2/d]
    !---------------------------------------------------------------------------
    real(RLEN),parameter    :: PERMIL=ONE/1000.0_RLEN
    real(RLEN),parameter    :: bar2atm  = ONE/1.01325_RLEN ! Conversion factor from bar to atm
    real(RLEN),parameter    :: Rgas_atm = 82.05736_RLEN    ! (cm3 * atm) / (mol * K)  CODATA (2006)
    real(RLEN),parameter    :: kwfac = 0.01_RLEN * HOURS_PER_DAY ! convert from cm/h to m/d
    real(RLEN)   :: pco2atm, patma, co2starair, co2star
    real(RLEN)   :: pH20, fco2atm, fugcoeff, bb, Del, xc2
    real(RLEN)   :: tk, tk100, tk1002, itk100
    real(RLEN)   :: pschmidt, Kw660, kwgas, K0
    ! CO2 solubility coefficients
    real(RLEN),parameter :: A(3) = (/-60.2409_RLEN, 93.4517_RLEN, 23.3585_RLEN/), &
                            B(3) = (/0.023517_RLEN, -0.023656_RLEN, 0.0047036_RLEN/)
    !! Schmidt coefficients
    real(RLEN),parameter :: Sc = 660.0_RLEN, &
         C(5)  = (/2116.8_RLEN, -136.25_RLEN, 4.7353_RLEN, -0.092307_RLEN, 0.0007555_RLEN/)
    !
    !---------------------------------------------------------------------------
    !
    ! common arrays
    tk = temp - ZERO_KELVIN
    tk100 = tk/100.0_RLEN
    tk1002 = tk100*tk100
    itk100 = 1.0_RLEN/tk100
    patma = patm * PERMIL * bar2atm ! convert patm from mbar to atm
    !
    ! Oceanic concentration [CO2*] from umol/kg to mol/kg
    co2star = co2ocn * 1.0e-6_RLEN 
    !
    ! Compute vapor pressure of seawater [in atm] (Weiss and Price, 1980, Eq. 10)
    pH20 = exp(24.4543_RLEN - 67.4509_RLEN*(100._RLEN/tk) - 4.8489_RLEN*log(tk/100.0_RLEN) - 0.000544_RLEN*salt)
    !
    ! Compute pco2atm [uatm] from xco2 [ppm], atmospheric pressure [atm], & vapor pressure of seawater pH20 [atm]
    pco2atm = (patma - pH20) * xco2
    ! 
    ! Compute fCO2atm [uatm] from pCO2atm [uatm] & fugacity coefficient [unitless]
    bb = -1636.75_RLEN + 12.0408_RLEN*tk - 0.0327957_RLEN*(tk*tk) + 0.0000316528_RLEN*(tk*tk*tk)
    Del = 57.7_RLEN - 0.118_RLEN*tk
    xc2 = (ONE - (pco2atm*1.e-6_RLEN) )**2
    fugcoeff = EXP( patma*(bb + 2.0d0*xc2*Del)/(Rgas_atm*tk) )
    fco2atm = pco2atm * fugcoeff
    !
    ! Surface solubility of gas K0  [(mol/kg) / atm] at T, S of surface water according to Weiss (1974)
    K0 = exp( A(1) + A(2)/tk100 + A(3)*log(tk100) + salt*(B(1) + B(2)*tk100 + B(3)*tk1002) )
    ! Equilbrium [CO2*air] for atm gas at Patm & sfc-water T,S [mol/kg]
    co2starair = K0 * fco2atm * 1.0e-6_RLEN

    ! Compute piston velocity kw660 (at 25 C) from wind speed (Wanninkhof 2014, Limnol. Oceanograph. Methods, 12, 351-362)
    kw660 = 0.251_RLEN * wind**2 * kwfac
    !  Schmidt number (Sc), ratio between the kinematic viscosity and the molecular diffusivity of the gas
    pschmidt = C(1) + C(2)*temp + C(3)*temp**2 + C(4)*temp**3 + C(5)*temp**4
    ! Transfer velocity for gas in m/d (see equation [4] in OCMIP2 design document & OCMIP2 Abiotic HOWTO)
    kwgas = kw660 * (Sc/pschmidt)**0.5

    ! Air-sea gas flux [mmol/(m2 * day)]
    AirSeaCO2 = kwgas * (co2starair - co2star) * (ONE - Fice) * rho * 1000

    return 

  end function AirSeaCO2

! -----------------------------------------------------------------------------

  elemental function AirpGas(xgas,patm,temp,salt)
  !
    implicit none
    real(RLEN)                       :: AirpGas ! Atmospheric partial pressure of gas [uatm]
    !
    real(RLEN),intent(IN)            :: xgas    ! Mixing ratio of Gas in dry air [ppmv]
    real(RLEN),intent(IN)            :: patm    ! Atmospheric pressure [mbar]
    real(RLEN),intent(IN)            :: temp    ! in-situ temperature at surface
    real(RLEN),intent(IN)            :: salt    ! practical salinity at surface
    !---------------------------------------------------------------------------
    real(RLEN)   :: tk
    real(RLEN)   :: patma, pH20
    !
    real(RLEN),parameter    :: PERMIL=ONE/1000.0_RLEN
    real(RLEN),parameter    :: bar2atm  = ONE/1.01325_RLEN ! Conversion factor from bar to atm
    !---------------------------------------------------------------------------
    ! common arrays
    tk = temp - ZERO_KELVIN
    patma = patm * PERMIL * bar2atm ! convert patm from mbar to atm
    !
    ! Compute vapor pressure of seawater [in atm] (Weiss and Price, 1980, Eq. 10)
    pH20 = exp(24.4543_RLEN - 67.4509_RLEN*(100._RLEN/tk) - 4.8489_RLEN*log(tk/100.0_RLEN) - 0.000544_RLEN*salt)
    !
    ! Compute pGasatm [uatm] from xgas [ppm], atmospheric pressure [atm], & vapor pressure of seawater pH20 [atm]
    AirpGas = (patma - pH20) * xgas

    return
  end function AirpGas

end module AirSeaExchange

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
