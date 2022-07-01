#include "cppdefs.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Air-Sea flux Exchange
!
! DESCRIPTION
!   ! compute Air-Sea CO2 Exchange

! !INTERFACE
Module AirSeaExchange
!
! !USES:
  use global_mem, ONLY: RLEN,ONE,ZERO
  use constants,  ONLY: HOURS_PER_DAY,ZERO_KELVIN

  IMPLICIT NONE
!
! !AUTHORS
!   T. Lovato (CMCC) 2017
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2020 BFM System Team (bfm_st@cmcc.it)
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
  public AirSeaCO2, AirpGas

 contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  elemental function AirSeaCO2(xco2,patm,temp,salt,rho,wind,Fice,co2ocn)
  ! 
    implicit none
    real(RLEN)                       :: AirSeaCO2 
    !     
    real(RLEN),intent(IN)            :: xco2   ! Atmospheric mixing ratio of Gas [ppmv]
    real(RLEN),intent(IN)            :: patm   ! Atmospheric pressure [mbar]
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
    real(RLEN),parameter :: A(3) = (/-58.0931_RLEN, 90.5069_RLEN, 22.2940_RLEN/), &
                            B(3) = (/0.027766_RLEN, -0.025888_RLEN, 0.0050578_RLEN/)
    !! Schmidt coefficients
    real(RLEN),parameter :: Sc = 660.0_RLEN, &
         C(5)  = (/2116.8_RLEN, -136.25_RLEN, 4.7353_RLEN, -0.092307_RLEN, 0.0007555_RLEN/)
    ! This is for alternative fit of co2starair (see Orr et al, 2017 GMD)
    !real(RLEN)   :: phi0atm
    !real(RLEN),parameter :: F(7) = (/-160.7333, 215.4152, 89.8920, -1.47759, 0.029941, -0.027455, 0.0053407/)
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

    ! Alternative computation
    !! Solubility function for atmospheric CO2 saturation concentration (Orr et al. GMD 2017, Eq. 15)
    !phi0atm = exp (F(1) + F(2)*itk100 + F(3)*log(tk100) + F(4)*tk1002 + salt*(F(5) + F(6)*tk100 + F(7)*tk1002))
    !! Compute saturation concentration in atmosphere [mol/kg] at total pressure Pa as in Orr et al (GMD 2017, Eq. 16)
    !co2starair = patma * phi0atm * xco2 * 1.0e-6_RLEN
    ! 
    ! Compute piston velolicty kw660 (at 25 C) from wind speed (Wanninkhof 2014, Limnol. Oceanograph. Methods, 12, 351-362)
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
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
