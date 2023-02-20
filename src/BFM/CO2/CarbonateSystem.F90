!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! MODULE: CO2System
!
! DESCRIPTION
! Solve for pH and other carbonate system variables using total alkalinity and
! dissolved inorganic carbon concentrations
! Use SolveSAPHE v1.0.1 routines from Munhoven (2013, GMD) modified to use local Ks instead of its own
! INPUT:
!   temp    = in situ temperature                  [degrees C]
!   salt    = practical salinity                   [unitless]
!   ta      = total alkalinity                     [ueq/kg]
!   tc      = dissolved inorganic carbon           [umol/kg]
!   pt      = total dissolved inorganic phosphorus [mmol/m3]
!   sit     = total dissolved inorganic silicon    [mmol/m3]
!   Bt      = total dissolved inorganic boron      computed 
!   St      = total dissolved inorganic sulfur     computed 
!   Ft      = total dissolved inorganic fluorine   computed 
!   K's     = K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
!   Patm    = atmospheric pressure [mbar=hPa]
!   rho     = in-situ densiaty     [kg/m3]
!   pr_in   = hydrostatic pressure [dbar]
!   If pr_in is given, considers in situ T & total pressure (atm + hydrostatic) to compute fCO2 and pCO2
!   otherwise, in situ T & only atm pressure (hydrostatic=0) for 'zero order' fCO2 and pCO2.
!
!   ---------
!
! OUTPUT:
!   ph   = pH on total scale
!   pco2 = CO2 partial pressure (uatm)
!   fco2 = CO2 fugacity (uatm)
!   co2  = aqueous CO2 concentration in [umol/kg]
!   hco3 = bicarbonate (HCO3-) concentration in [umol/kg]
!   co3  = carbonate (CO3--) concentration in [umol/kg]
!   OmegaA = Omega for aragonite, i.e., the aragonite saturation state
!   OmegaC = Omega for calcite, i.e., the   calcite saturation state
!
! ORIGINAL REFERENCES
!   Orr and Epitaloni, 2015: "Improved routines to model the ocean
!   carbonate system: mocsy 2.0." GMD 8.3 : 485-499.
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

#if defined INCLUDE_PELCO2 || defined INCLUDE_BENCO2

!
! INTERFACE
  module mem_CSYS
!
! USES
  use global_mem, ONLY: RLEN, LOGUNIT, ONE, ZERO
  use mem_CO2, ONLY: MaxIterPHsolver

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Shared Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),public  ::  K0 ! solubility : [Co2]=k0Ac Pco2
  real(RLEN),public  ::  K1 ! carbonate equilibrium I
  real(RLEN),public  ::  K2 ! carbonate equilibrium II
  real(RLEN),public  ::  Kw ! water dissociation
  real(RLEN),public  ::  Kb ! constant for Boron equilibrium
  real(RLEN),public  ::  Ks ! constant for bisulphate equilibrium
  real(RLEN),public  ::  Kf ! constant for hidrogen fluoride equilibirum
  real(RLEN),public  ::  K1p,K2p,K3p ! constants for phosphate equilibirum
  real(RLEN),public  ::  Ksi ! constant for silicic acid equilibrium
  real(RLEN),public  ::  Kspc ! constant for calcite equilibrium
  real(RLEN),public  ::  Kspa ! constant for aragonite equilibrium

  public CarbonateSystem

  contains

  function CarbonateSystem(salt,temp,rho,n1p,n5s,dic,alk,      &
                   co2,hco3,co3,pH,pco2,patm,pr_in,OmegaC,OmegaA,fCO2)
  ! USES
  use constants,  ONLY : ZERO_KELVIN,MW_C,MW_Ca,Rgas
  use mem_Param,  ONLY : p_atm0

  IMPLICIT NONE

  ! INPUT
  integer :: CarbonateSystem
  real(RLEN),intent(IN)            :: salt ! practical salinity
  real(RLEN),intent(IN)            :: temp ! in-situ temperature
  real(RLEN),intent(IN)            :: rho  ! in-situ density
  real(RLEN),intent(IN)            :: n1p
  real(RLEN),intent(IN)            :: n5s
  real(RLEN),intent(IN)            :: dic
  real(RLEN),intent(IN)            :: alk
  real(RLEN),intent(IN),optional   :: patm  ! Atmospheric pressure  [atm]
  real(RLEN),intent(IN),optional   :: pr_in ! Water column Pressure [dbar]
  ! INPUT/OUTPUT
  real(RLEN),intent(INOUT)         :: pH
  ! OUTPUT
  real(RLEN),intent(OUT)           :: co2
  real(RLEN),intent(OUT)           :: hco3
  real(RLEN),intent(OUT)           :: co3
  real(RLEN),intent(OUT)           :: pco2
  real(RLEN),intent(OUT),optional  :: OmegaC
  real(RLEN),intent(OUT),optional  :: OmegaA
  real(RLEN),intent(OUT),optional  :: fCO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),parameter    :: MEG=1.E6_RLEN, &
                             PERMIL=ONE/1000.0_RLEN,&
                             PERMEG=ONE/MEG
  real(RLEN)              :: press    = 0.0_RLEN         ! local hydrostatic pressure in bar
  real(RLEN)              :: pratm    = 0.0_RLEN         ! local atmospheric pressure in bar
  real(RLEN),parameter    :: Rgas_atm = 82.05736_RLEN    ! (cm3 * atm) / (mol * K)  CODATA (2006)
  real(RLEN),parameter    :: vbarCO2  = 32.3_RLEN        ! partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
  real(RLEN),parameter    :: CaRelCon = 0.02128_RLEN     ! Calcium ion relative concentration See Dickson (2007) and Munhoven (2013)
  real(RLEN),parameter    :: bar2atm  = ONE/1.01325_RLEN ! Conversion factor from bar to atm

  integer     :: i, l, error
  real(RLEN)  :: H, Hi, ta, tc, sit, pt, tfco2
  real(RLEN)  :: Bt, Ft, St
  real(RLEN)  :: lnk, Ks_0p, Kf_0p, deltav, deltak
  real(RLEN)  :: tk, tk100, tk1002, temp2, invtk, dlogtk, is, is2, sqrtis
  real(RLEN)  :: s, s2, sqrts, s15, scl 
  real(RLEN)  :: ptot, pr, pr2
  real(RLEN)  :: total2SWS, SWS2total, free2SWS, free2SWS_0p,total2SWS_0p, &
                  total2free, total2free_0p
  ! Carbonate Alkalinity 
  real(RLEN)  :: HSO4, HF, HSI, HPO4, ab, aw, ac, cu, cb, cc, Ca
  ! Fugacity
  real(RLEN)  :: B, Del, xc2, xCO2approx, fugcoeff
  ! Pressure correction
  real(RLEN), DIMENSION(12) :: a0, a1, a2, b0, b1, b2, lnkpok0
  ! 
  ! Constants in formulation for Pressure effect on K's (Millero, 95)
  ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
  !            9) K1P, 10) K2P, 11) K3P, 12) Ksi
  DATA a0 / 25.50_RLEN, 15.82_RLEN, 29.48_RLEN, 25.60_RLEN, 18.03_RLEN,  9.78_RLEN, &
            48.76_RLEN, 46.00_RLEN, 14.51_RLEN, 23.12_RLEN, 26.57_RLEN, 29.48_RLEN /
  DATA a1 / 0.1271_RLEN, -0.0219_RLEN, 0.1622_RLEN, 0.2324_RLEN, 0.0466_RLEN, -0.0090_RLEN, &
            0.5304_RLEN,  0.5304_RLEN, 0.1211_RLEN, 0.1758_RLEN, 0.2020_RLEN,  0.1622_RLEN /
  DATA a2 / 0.0_RLEN, 0.0_RLEN, 2.608_RLEN, -1.409_RLEN, 0.316_RLEN, -0.942_RLEN, &
            0.0_RLEN, 0.0_RLEN,-0.321_RLEN, -2.647_RLEN, -3.042_RLEN, -2.608_RLEN /
  DATA b0 /  3.08_RLEN, -1.13_RLEN, 2.84_RLEN, 5.13_RLEN, 4.53_RLEN, 3.91_RLEN, &
            11.76_RLEN, 11.76_RLEN, 2.67_RLEN, 5.15_RLEN, 4.08_RLEN, 2.84_RLEN /
  DATA b1 / 0.0877_RLEN, -0.1475_RLEN, 0.0_RLEN,  0.0794_RLEN, 0.09_RLEN,  0.054_RLEN, &
            0.3692_RLEN,  0.3692_RLEN, 0.0427_RLEN, 0.09_RLEN, 0.0714_RLEN, 0.0_RLEN  /
  !DATA b2 / 12*0.0_RLEN / ! not used as it is zero
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  CarbonateSystem = 0

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! 1. PREPARE INPUT FIELDS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Convert input fields
  ! from umol/kg to mol/kg
  ta = alk * PERMEG
  tc = dic * PERMEG
  ! from mmol/m^3 -> mol/kg
  pt  = n1p / rho * PERMIL
  sit = n5s / rho * PERMIL
  
  ! Absolute temperature (Kelvin) and related values
  tk     = temp - ZERO_KELVIN  !ZERO_KELVIN=-273.16; ETW is in degC; tk is in degK
  temp2  = temp * temp
  tk100  = tk / 100.0_RLEN
  tk1002 = tk100 * tk100
  invtk  = ONE / tk
  dlogtk = log(tk)

  ! Salinity and simply related values
  s  = salt
  s2 = salt*salt
  sqrts = sqrt(salt)
  s15 = salt**1.5_RLEN
 
  ! Hydrostatic Pressure [dbar]
  if (present(pr_in)) then
    press = pr_in * 0.1_RLEN  ! convert from dbar to bar
    pr2   = press * press / Rgas
    pr    = press / Rgas
  endif
  
  ! Atmospheric pressure
  pratm = p_atm0 * PERMIL     ! convert reference atmospheric pressure from mbar to bar
  if ( present(patm) ) pratm = patm * PERMIL  ! convert from mbar to bar
  Ptot  = (pratm + press) * bar2atm   !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
  ! 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! 2. COMPUTE CONSTANTS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! chlorinity
  scl = s/1.80655_RLEN

  ! ionic strength
  is = 19.924_RLEN*s/ (1000._RLEN-1.005_RLEN*s)
  is2 = is*is
  sqrtis = sqrt(is)
  
  ! Total concentrations for sulfate, fluoride, and boron
  ! Sulfate: Morris & Riley (1966)
  St = 0.14_RLEN * scl/96.062_RLEN
  ! Fluoride: Riley (1965)
  Ft = 0.000067_RLEN * scl/18.9984_RLEN
  ! Boron: Lee et al (2010)
  Bt = 0.0002414_RLEN * scl/10.811_RLEN

  ! -----------------------------------------------------------------------
  ! K0, solubility of co2 in the water (K Henry) from Weiss 1974
  ! K0 = [CO2]/ fCO2 [mol/kg/atm]
  ! -----------------------------------------------------------------------
  lnK = 93.4517_RLEN/tk100 - 60.2409_RLEN + 23.3585_RLEN * log(tk100) +   &
       s * (0.023517_RLEN - 0.023656_RLEN * tk100 + 0.0047036_RLEN * tk1002)
  K0 = exp ( lnK )

  ! -----------------------------------------------------------------------
  ! Choice of Acidity constants
  ! K1 = [H][HCO3]/[H2CO3]   ,   K2 = [H][CO3]/[HCO3]
  ! -----------------------------------------------------------------------
  ! Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total scale)
  K1 = 10.0_RLEN**(-ONE*(3633.86_RLEN*invtk - 61.2172_RLEN + 9.6777_RLEN*dlogtk - &
       0.011555_RLEN * s + 0.0001152_RLEN * s2))
  K2 = 10.0_RLEN**(-ONE*(471.78_RLEN*invtk + 25.9290_RLEN - &
       3.16967_RLEN*dlogtk - 0.01781_RLEN * s + 0.0001122_RLEN * s2))
  ! ! Millero (2010, Mar. Fresh Wat. Res.) (seawater scale)
  ! K1 = 10.0_RLEN**(-ONE*( (6320.813_RLEN*invtk + 19.568224_RLEN*dlogtk -126.34048_RLEN + &
  !      13.4038_RLEN*sqrts + 0.03206_RLEN*s - (5.242e-5_RLEN)*s2) + &
  !      (-530.659_RLEN*sqrts - 5.8210_RLEN*s)*invtk -2.0664_RLEN*sqrts*dlogtk) )
  ! K2 = 10.0_RLEN**(-ONE*( (5143.692_RLEN*invtk + 14.613358_RLEN*dlogtk -90.18333_RLEN + &
  !      21.3728_RLEN*sqrts + 0.1218_RLEN*s - (3.688e-4_RLEN)*s2 ) + &
  !      (-788.289_RLEN*sqrts - 19.189_RLEN*s)*invtk -3.374_RLEN*sqrts*dlogtk) )
  ! ! Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014) (seawater scale)
  ! K1 = 10.0_RLEN**(-ONE*( (6320.813_RLEN*invtk + 19.568224_RLEN*dlogtk -126.34048_RLEN + &
  !      13.409160_RLEN*sqrts + 0.031646_RLEN*s - (5.1895e-5_RLEN)*s2 ) + &
  !      (-531.3642_RLEN*sqrts - 5.713_RLEN*s)*invtk -2.0669166_RLEN*sqrts*dlogtk) )
  ! K2 = 10.0_RLEN**(-ONE*( (5143.692_RLEN*invtk + 14.613358_RLEN*dlogtk -90.18333_RLEN + & 
  !      21.225890_RLEN*sqrts + 0.12450870_RLEN*s - (3.7243e-4_RLEN)*s2 ) + &
  !      (-779.3444_RLEN*sqrts - 19.91739_RLEN*s)*invtk -3.3534679_RLEN*sqrts*dlogtk) )
 
  !-----------------------------------------------------------------------
  ! Kb = [H][BO2]/[HBO2] 
  ! Millero p.669 (1995) using data from Dickson (1990)    (total scale)
  !-----------------------------------------------------------------------
  lnK = (-8966.90_RLEN - 2890.53_RLEN*sqrts - 77.942_RLEN*s +    &
       1.728_RLEN*s15 - 0.0996_RLEN*s2)*invtk +                  &
       (148.0248_RLEN + 137.1942_RLEN*sqrts + 1.62142_RLEN*s) +  &
       (-24.4344_RLEN - 25.085_RLEN*sqrts - 0.2474_RLEN*s) *     &
       dlogtk + 0.053105_RLEN*sqrts*tk
  Kb = exp(lnK)

  ! -----------------------------------------------------------------------
  ! K1p = [H][H2PO4]/[H3PO4]
  ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
  ! Millero (1995), p.670, eq. 65                        (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = -4576.752_RLEN*invtk + 115.540_RLEN - 18.453_RLEN * dlogtk + &
       (-106.736_RLEN*invtk + 0.69171_RLEN) * sqrts +                &
       (-0.65643_RLEN*invtk - 0.01844_RLEN) * s
  K1p = exp(lnK)

  ! -----------------------------------------------------------------------
  ! K2p = [H][HPO4]/[H2PO4]
  ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
  ! Millero (1995), p.670, eq. 66                        (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = -8814.715_RLEN*invtk + 172.1033_RLEN - 27.927_RLEN * dlogtk +  &
       (-160.340_RLEN*invtk + 1.3566_RLEN) * sqrts +                   &
       (0.37335_RLEN*invtk - 0.05778_RLEN) * s
  K2p = exp(lnK)

  !------------------------------------------------------------------------
  ! K3p = [H][PO4]/[HPO4]
  ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  ! Millero (1995), p.670, eq. 67                        (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = -3070.75_RLEN*invtk - 18.126_RLEN +     &
       (17.27039_RLEN*invtk + 2.81197_RLEN) *   &
       sqrts + (-44.99486_RLEN*invtk - 0.09984_RLEN) * s
  K3p = exp(lnK)

  !------------------------------------------------------------------------
  ! ksi = [H][SiO(OH)3]/[Si(OH)4]
  ! Millero (1995), p.671, eq. 72                        (seawater scale)
  ! -----------------------------------------------------------------------
  ksi = -8904.2_RLEN*invtk + 117.400_RLEN - 19.334_RLEN * dlogtk + &
       (-458.79_RLEN*invtk + 3.5913_RLEN) * sqrtis +               &
       (188.74_RLEN*invtk - 1.5998_RLEN) * is +                    &
       (-12.1652_RLEN*invtk + 0.07871_RLEN) * is2 +                &
       log(ONE-0.001005_RLEN*s)
  Ksi = exp(lnK)

  ! -----------------------------------------------------------------------
  ! Kw = [H][OH]
  ! Millero (1995) p.670, eq. 63 from composite data     (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = 148.9802_RLEN -13847.26_RLEN*invtk - 23.6521_RLEN * dlogtk +      &
       (118.67_RLEN*invtk - 5.977_RLEN + 1.0495_RLEN * dlogtk) *          &
       sqrts - 0.01615_RLEN * s
  Kw = exp(lnK)

  !------------------------------------------------------------------------
  ! ks = [H][SO4]/[HSO4]
  ! Dickson (1990, J. chem. Thermodynamics 22, 113)          (free scale)
  !------------------------------------------------------------------------
  lnK = -4276.1_RLEN*invtk + 141.328_RLEN - 23.093_RLEN*dlogtk +          &
       (-13856._RLEN*invtk + 324.57_RLEN - 47.986_RLEN*dlogtk) * sqrtis + &
       (35474._RLEN*invtk - 771.54_RLEN + 114.723_RLEN*dlogtk) * is -     &
       2698._RLEN*invtk*is**1.5_RLEN + 1776._RLEN*invtk*is2 +             &
       log(ONE - 0.001005_RLEN*s)
  Ks_0p = exp(lnK)

  !------------------------------------------------------------------------
  ! kf = [H][F]/[HF]
  ! Perez & Fraga (1987) recom. by Dickson et al., (2007)   (total scale)
  !------------------------------------------------------------------------
  lnK = 874.0_RLEN*invtk - 9.68_RLEN + 0.111_RLEN*sqrts
  Kf_0p = exp(lnK)

  !------------------------------------------------------------------------
  ! Kspc = [Ca2+] [CO32-] - apparent solubility product of Calcite
  ! Mucci (1983)  [mol/kg-soln]
  !------------------------------------------------------------------------
  kspc = 10.0_RLEN**(ONE*( -171.9065_RLEN - 0.077993_RLEN * tk + 2839.319_RLEN * invtk &
          + 71.595_RLEN * log10(tk) + sqrts * (-0.77712_RLEN +  &
          0.0028426_RLEN * tk + 178.34_RLEN * invtk)            &
          - 0.07711_RLEN * s + 0.0041249_RLEN *s15 ) )
  !------------------------------------------------------------------------
  ! Kspa = [Ca2+] [CO32-] - apparent solubility product of Aragonite
  ! Mucci (1983)  [mol/kg-soln]
  !------------------------------------------------------------------------
  kspa = 10.0_RLEN**(ONE*(  -171.945_RLEN - 0.077993_RLEN * tk + 2903.293_RLEN * invtk   &
          + 71.595_RLEN * log10(tk) + sqrts * (-0.068393_RLEN +   &
          0.0017276_RLEN * tk + 88.135_RLEN * invtk)   &
          - 0.10018_RLEN * s + 0.0059415_RLEN * s15 ) ) 

  ! Pressure effect on K0 based on Weiss (1974, equation 5)
  K0 = K0 * exp( ((1-Ptot)*vbarCO2)/(Rgas_atm*tk) )

  lnkpok0 = 0.0_RLEN
  if (present(pr_in) .AND. press .GT. 0) then
  ! Pressure effect on all other K's (based on Millero, (1995)
  ! Index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
  !           K1p(9), K2p(10), K3p(11), Ksi(12)
     do l = 1, 12
        deltav  =  -a0(l) + a1(l)*temp + 1.0e-3_RLEN*a2(l)*temp2
        deltak  = 1.0e-3_RLEN*(-b0(l) + b1(l)*temp)
        lnkpok0(l) = -deltav*invtk*pr + 0.5_RLEN*deltak*invtk*pr2
     end do
  endif 
  
  ! Pressure correction on Ks (Free scale)
  Ks = Ks_0p*EXP(lnkpok0(5))
  ! Conversion factor total -> free scale
  total2free     = 1.0_RLEN/(1.0_RLEN + St/Ks)   ! Kfree = Ktotal*total2free
  ! Conversion factor total -> free scale at pressure zero
  total2free_0p  = 1.0_RLEN/(1.0_RLEN + St/Ks_0p)   ! Kfree = Ktotal*total2free

  ! Pressure correction on Kf
  ! Kf must be on FREE scale before correction
  Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
  Kf    = Kf_0p * EXP(lnkpok0(6)) !Pressure correction (on Free scale)
  Kf    = Kf / total2free         !Convert back from Free to Total scale

  ! Convert between seawater and total hydrogen (pH) scales
  free2SWS  = 1.0_RLEN + St/Ks + Ft/(Kf*total2free)  ! using Kf on free scale
  total2SWS = total2free * free2SWS                  ! KSWS = Ktotal*total2SWS
  SWS2total = 1.0_RLEN / total2SWS
  ! Conversion at pressure zero
  free2SWS_0p  = 1.d0 + St/Ks_0p + Ft/(Kf_0p)  ! using Kf on free scale
  total2SWS_0p = total2free_0p * free2SWS_0p         ! KSWS = Ktotal*total2SWS

  ! Convert from Total to Seawater scale before pressure correction
  ! Must change to SEAWATER scale: K1, K2, Kb
  K1 = K1 * total2SWS_0p
  K2 = K2 * total2SWS_0p
  Kb = Kb * total2SWS_0p
  ! Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw
  ! Other contants (keep on another scale):
  !    - K0         (independent of pH scale, already pressure corrected)
  !    - Ks         (already on Free scale;   already pressure corrected)
  !    - Kf         (already on Total scale;  already pressure corrected)
  !    - Kspc, Kspa (independent of pH scale; pressure-corrected below)
  
  ! Perform actual pressure correction (on seawater scale)
  K1   = K1   * EXP(lnkpok0(1))
  K2   = K2   * EXP(lnkpok0(2))
  Kb   = Kb   * EXP(lnkpok0(3))
  Kw   = Kw   * EXP(lnkpok0(4))
  Kspc = Kspc * EXP(lnkpok0(7))
  Kspa = Kspa * EXP(lnkpok0(8))
  K1p  = K1p  * EXP(lnkpok0(9))
  K2p  = K2p  * EXP(lnkpok0(10))
  K3p  = K3p  * EXP(lnkpok0(11))
  Ksi  = Ksi  * EXP(lnkpok0(12))

  ! Convert back to original total scale:
  K1  = K1  * SWS2total
  K2  = K2  * SWS2total
  K1p = K1p * SWS2total
  K2p = K2p * SWS2total
  K3p = K3p * SWS2total
  Kb  = Kb  * SWS2total
  Ksi = Ksi * SWS2total
  Kw  = Kw  * SWS2total

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! 3. COMPUTE CARBONATE SYSTEM
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! H+ concentration (mol/kg) at previous step
  if ( pH <= 0._RLEN ) then
     Hi = Hini_for_at(ta,tc,bt,K1,K2,Kb)
  else
     Hi = 10.0_RLEN**(-pH)
  endif

  ! Solve for H+ using above result as the initial H+ value (mol/kg)
  H = solve_at_general(ta, tc, Bt, pt, sit, St, Ft,            &
               K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi, Hi )
  if (H < ZERO) CarbonateSystem = 1

  ! Calculate pH from from H+ concentration (mol/kg)
  pH = -ONE * LOG10( H )

  ! Compute carbonate Alk (Ac) by difference: from total Alk and other Alk components
  HSO4 = St/(1.0_RLEN + Ks/(H/(1.0_RLEN + St/Ks)))
  HF   = 1.0_RLEN / (1.0_RLEN + Kf/H)
  HSI  = 1.0_RLEN / (1.0_RLEN + H/Ksi)
  HPO4 = (K1p*K2p*(H + 2.*K3p) - H**3)    /         &
         (H**3 + K1p*H**2 + K1p*K2p*H + K1p*K2p*K3p)
  ab = Bt/(1.0_RLEN + H/Kb)
  aw = Kw/H - H/(1.0_RLEN + St/Ks)
  ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4

  ! Calculate CO2*, HCO3-, & CO32- (in mol/kg soln) from Ct, Ac, H+, K1, & K2
  cu = (2.0_RLEN * tc - ac) / (2.0_RLEN + K1 / H)
  cb = K1 * cu / H
  cc = K2 * cb / H 
  
  ! Determine Omega Calcite and Aragonite (see Munhoven (2013, GMD))
  Ca = (CaRelCon / MW_Ca) * s/1.80655_RLEN
  if (present(OmegaA)) OmegaA = (Ca*cc) / Kspa
  if (present(OmegaC)) OmegaC = (Ca*cc) / Kspc

  ! Determine CO2 fugacity [uatm]
  tfco2 = cu * 1.e6_RLEN/K0
  if (present(fCO2)) fCO2 = tfco2

  ! Determine CO2 partial pressure from CO2 fugacity [uatm]
  ! compute fugacity coefficient terms : B, Del, xc2
  B = -1636.75_RLEN + 12.0408_RLEN*tk - 0.0327957_RLEN*(tk*tk) + 0.0000316528_RLEN*(tk*tk*tk)
  Del = 57.7_RLEN - 0.118_RLEN*tk
  xCO2approx = tfco2 * 1.e-6_RLEN

  if (present(pr_in) .AND. press .GT. 0) &
     xCO2approx = xCO2approx * exp( ((ONE-Ptot)*32.3_RLEN)/(82.05736_RLEN*tk) )   ! of K0 press. correction, see Weiss (1974, equation 5)
  xc2 = (ONE - xCO2approx)**2
  fugcoeff = exp( Ptot*(B + 2.0_RLEN*xc2*Del)/(Rgas_atm*tk) )
  pco2 = tfco2 / fugcoeff

  ! scale from mol/kg -----> umol/kg
  co2  = cu * MEG 
  hco3 = cb * MEG
  co3  = cc * MEG

  return
  end function CarbonateSystem

! ----------------

  elemental FUNCTION Hini_for_at(p_alkcb, p_dictot, p_bortot, K1, K2, Kb)
  ! Function returns the root for the 2nd order approximation of the
  ! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
  ! around the local minimum, if it exists.
  
  ! Returns * 1E-03_RLEN if p_alkcb <= 0
  !         * 1E-10_RLEN if p_alkcb >= 2*p_dictot + p_bortot
  !         * 1E-07_RLEN if 0 < p_alkcb < 2*p_dictot + p_bortot
  !                    and the 2nd order approximation does not have a solution 
    IMPLICIT NONE
    REAL(KIND=RLEN) :: Hini_for_at   ! [H+] in mol/kg
    
    ! Argument variables
    !--------------------
    REAL(KIND=RLEN), INTENT(IN)   ::  p_alkcb, p_dictot, p_bortot ! [mol/kg]
    REAL(KIND=RLEN), INTENT(IN)   ::  K1, K2, Kb
    
    ! Local variables
    !-----------------
    REAL(KIND=RLEN)  ::  p_hini
    REAL(KIND=RLEN)  ::  zca, zba
    REAL(KIND=RLEN)  ::  zd, zsqrtd, zhmin
    REAL(KIND=RLEN)  ::  za2, za1, za0
    
    IF (p_alkcb <= 0._RLEN) THEN
      p_hini = 1.e-3_RLEN
    ELSEIF (p_alkcb >= (2._RLEN*p_dictot + p_bortot)) THEN
      p_hini = 1.e-10_RLEN
    ELSE
      zca = p_dictot/p_alkcb
      zba = p_bortot/p_alkcb
    
      ! Coefficients of the cubic polynomial
      za2 = Kb*(1._RLEN - zba) + K1*(1._RLEN-zca)
      za1 = K1*Kb*(1._RLEN - zba - zca) + K1*K2*(1._RLEN - (zca+zca))
      za0 = K1*K2*Kb*(1._RLEN - zba - (zca+zca))
                                            ! Taylor expansion around the minimum
      zd = za2*za2 - 3._RLEN*za1              ! Discriminant of the quadratic equation
                                            ! for the minimum close to the root
    
      IF(zd > 0._RLEN) THEN                   ! If the discriminant is positive
        zsqrtd = SQRT(zd)
        IF(za2 < 0._RLEN) THEN
          zhmin = (-za2 + zsqrtd)/3._RLEN
        ELSE
          zhmin = -za1/(za2 + zsqrtd)
        ENDIF
        p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
      ELSE
        p_hini = 1.e-7_RLEN
      ENDIF
    
    Hini_for_at = p_hini 

    ENDIF
    RETURN
  END FUNCTION hini_for_at
  
! ----------------

  elemental FUNCTION solve_at_general(p_alktot, p_dictot, p_bortot,               &
                              p_po4tot, p_siltot,                                 &
                              p_so4tot, p_flutot,                                 &
                              K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,     &
                              p_h)
    
    ! Purpose: Compute [H°] ion concentration from sea-water ion concentrations,
    !          alkalinity, DIC, and equilibrium constants
    ! Universal pH solver that converges from any given initial value,
    ! determines upper an lower bounds for the solution if required
    
    IMPLICIT NONE
    REAL(KIND=RLEN) :: SOLVE_AT_GENERAL    ! [H+] in mol/kg
    
    ! Argument variables
    !--------------------
    REAL(KIND=RLEN), INTENT(IN)            :: p_alktot  ! [mol/kg] 
    REAL(KIND=RLEN), INTENT(IN)            :: p_dictot  ! [mol/kg]
    REAL(KIND=RLEN), INTENT(IN)            :: p_bortot  ! [mol/kg]
    REAL(KIND=RLEN), INTENT(IN)            :: p_po4tot  ! [mol/kg]
    REAL(KIND=RLEN), INTENT(IN)            :: p_siltot  ! [mol/kg]
    REAL(KIND=RLEN), INTENT(IN)            :: p_so4tot  ! [mol/kg]
    REAL(KIND=RLEN), INTENT(IN)            :: p_flutot  ! [mol/kg]
    REAL(KIND=RLEN), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
    REAL(KIND=RLEN), INTENT(IN)            :: K1p, K2p, K3p, Ksi
    REAL(KIND=RLEN), INTENT(IN)            :: p_h       ! [H+] in mol/kg
    
    ! Local variables
    !-----------------
    REAL(KIND=RLEN)  ::  zh, zh_prev, zh_lnfactor
    REAL(KIND=RLEN)  ::  zalknw_inf, zalknw_sup
    REAL(KIND=RLEN)  ::  zh_min, zh_max
    REAL(KIND=RLEN)  ::  zdelta, zh_delta
    REAL(KIND=RLEN)  ::  zeqn, zdeqndh, zeqn_absmin
    REAL(KIND=RLEN)  ::  aphscale
    LOGICAL          ::  l_exitnow
    INTEGER          ::  niter_atgen
    REAL(KIND=RLEN)  ::  znumer, zdenom, zdnumer
    REAL(KIND=RLEN)  ::  zalk_dic, zalk_bor, zalk_po4, zalk_sil, zalk_so4, zalk_flu, zalk_wat
    REAL(KIND=RLEN)  ::  zdalk_dic, zdalk_bor, zdalk_po4, zdalk_sil, zdalk_so4, zdalk_flu, zdalk_wat

    ! General parameters
    !-----------------
    REAL(KIND=RLEN), PARAMETER :: pz_exp_threshold = 1.0_RLEN
    REAL(KIND=RLEN), PARAMETER :: pp_rdel_ah_target = 1.E-8_RLEN
    
    
    ! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
    aphscale = 1._RLEN + p_so4tot/Ks

    ! lower and upper bounds of "non-water-selfionization"
    ! contributions to total alkalinity (the infimum and the supremum), i.e
    ! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])
    zalknw_inf = -p_po4tot - p_so4tot - p_flutot
    zalknw_sup = p_dictot + p_dictot + p_bortot + p_po4tot + p_po4tot + p_siltot

    zdelta = (p_alktot-zalknw_inf)**2 + 4._RLEN*Kw/aphscale
    
    IF(p_alktot >= zalknw_inf) THEN
       zh_min = 2._RLEN*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
    ELSE
       zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._RLEN
    ENDIF
    
    zdelta = (p_alktot-zalknw_sup)**2 + 4._RLEN*Kw/aphscale
    
    IF(p_alktot <= zalknw_sup) THEN
       zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._RLEN
    ELSE
       zh_max = 2._RLEN*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
    ENDIF
    
    zh = MAX(MIN(zh_max, p_h), zh_min)
    niter_atgen        = 0                 ! Reset counters of iterations
    zeqn_absmin        = HUGE(1._RLEN)

    DO
       IF(niter_atgen >= MaxIterPHsolver ) THEN
          zh = -1._RLEN
          EXIT
       ENDIF

       zh_prev = zh

       ! Compute total alkalinity from ion concentrations and equilibrium constants
       ! H2CO3 - HCO3 - CO3 : n=2, m=0
       znumer     = 2._RLEN*K1*K2 + zh*       K1
       zdenom     =         K1*K2 + zh*(      K1 + zh)
       zalk_dic   = p_dictot * (znumer/zdenom)
       zdnumer    = K1*K1*K2 + zh*(4._RLEN*K1*K2 + zh* K1 ) 
       zdalk_dic  = -p_dictot * (zdnumer/zdenom**2)
       ! B(OH)3 - B(OH)4 : n=1, m=0
       znumer     = Kb
       zdenom     = Kb + zh
       zalk_bor   = p_bortot * (znumer/zdenom)
       zdnumer    = Kb
       zdalk_bor  = -p_bortot * (zdnumer/zdenom**2)
       ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
       znumer     = 3._RLEN*K1p*K2p*K3p + zh*(2._RLEN*K1p*K2p + zh* K1p)
       zdenom     =         K1p*K2p*K3p + zh*(        K1p*K2p + zh*(K1p + zh))
       zalk_po4   = p_po4tot * (znumer/zdenom - 1._RLEN) ! Zero level of H3PO4 = 1
       zdnumer    = K1p*K2p*K1p*K2p*K3p + zh*(4._RLEN*K1p*K1p*K2p*K3p            &
                                        + zh*(9._RLEN*K1p*K2p*K3p + K1p*K1p*K2p  &
                                        + zh*(4._RLEN*K1p*K2p                    &
                                        + zh*       K1p)))
       zdalk_po4  = -p_po4tot * (zdnumer/zdenom**2)
       ! H4SiO4 - H3SiO4 : n=1, m=0
       znumer     = Ksi
       zdenom     = Ksi + zh
       zalk_sil   = p_siltot * (znumer/zdenom)
       zdnumer    = Ksi
       zdalk_sil  = -p_siltot * (zdnumer/zdenom**2)
       ! HSO4 - SO4 : n=1, m=1
       znumer     = Ks
       zdenom     = Ks + zh
       zalk_so4   = p_so4tot * (znumer/zdenom - 1._RLEN)
       zdnumer    = Ks
       zdalk_so4  = -p_so4tot * (zdnumer/zdenom**2)
       ! HF - F : n=1, m=1
       znumer     = Kf
       zdenom     = Kf + zh
       zalk_flu   = p_flutot * (znumer/zdenom - 1._RLEN)
       zdnumer    = Kf
       zdalk_flu  = -p_flutot * (zdnumer/zdenom**2)
       ! H2O - OH
       zalk_wat   = Kw/zh - zh/aphscale
       zdalk_wat  = -Kw/zh**2 - 1._RLEN/aphscale
        
       zeqn = zalk_dic + zalk_bor + zalk_po4 + zalk_sil   &
               + zalk_so4 + zalk_flu + zalk_wat - p_alktot 
      
       zdeqndh = zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
               + zdalk_so4 + zdalk_flu + zdalk_wat

    
       ! Adapt bracketing interval
       IF(zeqn > 0._RLEN) THEN
          zh_min = zh_prev
       ELSEIF(zeqn < 0._RLEN) THEN
          zh_max = zh_prev
       ELSE
          ! zh is the root; unlikely but, one never knows
          EXIT
       ENDIF
       ! Now determine the next iterate zh
       niter_atgen = niter_atgen + 1

       IF(ABS(zeqn) >= 0.5_RLEN*zeqn_absmin) THEN
          zh = SQRT(zh_max * zh_min)
          zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
       ELSE
          zh_lnfactor = -zeqn/(zdeqndh*zh_prev)
          IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
             zh          = zh_prev*EXP(zh_lnfactor)
          ELSE
             zh_delta    = zh_lnfactor*zh_prev
             zh          = zh_prev + zh_delta
          ENDIF
          IF( zh < zh_min ) THEN ! if [H]_new < [H]_min
             zh                = SQRT(zh_prev * zh_min)
             zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
          ENDIF
          IF( zh > zh_max ) THEN ! if [H]_new > [H]_min
             zh                = SQRT(zh_prev * zh_max)
             zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
          ENDIF
       ENDIF
       zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)
       l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)
    
       IF(l_exitnow) EXIT
    ENDDO
    
    solve_at_general = zh
    
    
  END FUNCTION solve_at_general

end module mem_CSYS

#endif

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
