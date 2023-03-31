#include "fabm_driver.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelCO2Dynamics
!
! DESCRIPTION
!   !
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
module bfm_PelagicCSYS

   use fabm_types
   use fabm_builtin_models
   use ogs_bfm_shared

   implicit none
   
   private

   type,extends(type_base_model),public :: type_ogs_bfm_PelagicCSYS
!     Variable identifiers
      type (type_state_variable_id)     :: id_O3c,id_O3h,id_N1p,id_N5s
      type (type_dependency_id)         :: id_ETW, id_ESW, id_ERHO, id_EPR,id_depth

      type (type_horizontal_dependency_id) :: id_EWIND,id_PCO2A, id_EPRatm

      type (type_diagnostic_variable_id) :: id_pH,id_pco2,id_CarbA, id_BiCarb, id_Carb
      type (type_diagnostic_variable_id) :: id_OCalc,id_OArag,id_ALK,id_DIC,id_EPRdiag !,id_depth_diag

      type (type_horizontal_diagnostic_variable_id) :: id_fair,id_wnd_diag

   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type type_ogs_bfm_PelagicCSYS

   public :: CarbonateSystem  

contains

   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_PelagicCSYS), intent(inout), target :: self
      integer,                      intent(in)            :: configunit
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%register_state_variable(self%id_O3c,'c','mgC/m^3','total dissolved inorganic carbon', 26400._rk,minimum=10000._rk)
!      call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_O3c) !not used yet and attention, scale_factor must be used (standard_variable for C is mmol/m3)
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      call self%register_state_dependency(self%id_N5s,'N5s','mmol Si/m^3','Silicate')
      call self%register_state_dependency(self%id_O3h,'O3h','mmol /m^3','alkalinity')

! dependency for carb sys solver
      call self%register_dependency(self%id_ETW, standard_variables%temperature)
      call self%register_dependency(self%id_ESW,standard_variables%practical_salinity)
      call self%register_dependency(self%id_ERHO,standard_variables%density)
      call self%register_dependency(self%id_depth,standard_variables%depth)  
      call self%register_dependency(self%id_EPR,standard_variables%pressure)  
      call self%register_dependency(self%id_EPRatm,standard_variables%surface_air_pressure)  
! dependecy for previous values of carb syst variables
!      call self%register_dependency(self%id_pco2_in,'pCO2','1e-6','previous pCO2')
!      call self%register_dependency(self%id_Carb_in,'Carb','mmol/m^3','previous carbonate concentration')
!      call self%register_dependency(self%id_CarbA_in,'CarbA','mmol/m^3','previous carbonic acid concentration')
!      call self%register_dependency(self%id_BiCarb_in,'BiCarb','mmol/m^3','previous bicarbonate concentration')
!      call self%register_dependency(self%id_pH_in,'pH','-','previous pH')
! dependency for computing CO2 air-sea exchange
      call self%register_dependency(self%id_EWIND,standard_variables%wind_speed)
      call self%register_dependency(self%id_PCO2A,standard_variables%mole_fraction_of_carbon_dioxide_in_air) ! in [ppm]
!diagnostics for carbonate system
      call self%register_diagnostic_variable(self%id_ALK,  'ALK',  'umol/kg','alkalinity [umol/kg]',missing_value=0._rk)
      call self%register_diagnostic_variable(self%id_DIC,  'DIC',  'umol/kg','DIC [umol/kg]',missing_value=0._rk)
      call self%register_diagnostic_variable(self%id_ph,'pH',    '-','pH on total scale',standard_variable=standard_variables%ph_reported_on_total_scale,missing_value=0._rk)

      call self%register_diagnostic_variable(self%id_pco2,  'pCO2',  'ppm','partial pressure of CO2',missing_value=0._rk)
      call self%register_diagnostic_variable(self%id_CarbA, 'CarbA', 'mmol/m^3','carbonic acid concentration',missing_value=0._rk)
      call self%register_diagnostic_variable(self%id_BiCarb,'BiCarb','mmol/m^3','bicarbonate concentration',missing_value=0._rk)
      call self%register_diagnostic_variable(self%id_Carb,  'Carb','mmol/m^3','carbonate concentration',standard_variable=standard_variables%mole_concentration_of_carbonate_expressed_as_carbon,missing_value=0._rk)

      call self%register_diagnostic_variable(self%id_OCalc,'OCalc','-','calcite saturation',missing_value=4._rk)
      call self%register_diagnostic_variable(self%id_OArag,'OArag','-','aragonite saturation',missing_value=3._rk)
      call self%register_diagnostic_variable(self%id_EPRdiag,'EPR','dbar','pressure water',missing_value=0._rk)
!      call self%register_diagnostic_variable(self%id_depth_diag,'depth','m','depth water column cell',missing_value=0._rk)

! diagnostics for air-sea CO2 flux
      call self%register_diagnostic_variable(self%id_fair,'fair','mmolC/m^2/d','air-sea flux of CO2',source=source_do_surface)
      call self%register_diagnostic_variable(self%id_wnd_diag,'wind','m/s','surface windspeed',source=source_do_surface)


      
      self%dt = 3600._rk*24._rk

   end subroutine
!==================================================================
  subroutine do(self,_ARGUMENTS_DO_)
      class (type_ogs_bfm_PelagicCSYS), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: O3c,ETW,ESW,ERHO,EPR,EPRatm,Depth
      real(rk) :: O3h,N1p,N5s 
      real(rk) :: DIC,ALK,PHOS,SILIC

      real(rk) :: pH,PCO2,H2CO3,HCO3,CO3,CO2
      real(rk) :: OCalc,OArag, ffCO2
      integer  :: error

      _LOOP_BEGIN_
         _GET_(self%id_O3c,O3c)  ! dissolved inorganic carbon in mg/m3
         _GET_(self%id_O3h,O3h)  ! alkalinity in mmol/m3
         _GET_(self%id_N1p,N1p)  ! 
         _GET_(self%id_N5s,N5s)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_ESW,ESW)
         _GET_(self%id_ERHO,ERHO)! density kg/m3
         _GET_(self%id_EPR,EPR)  ! pressure water in dbar
         _GET_(self%id_depth,Depth)  ! depth water column cell in m
         _GET_HORIZONTAL_(self%id_EPRatm,EPRatm)  ! pressure water in dbar

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Compute carbonate system equilibria
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! To use the Pressure correction of CSYS here the pr_in=EPS value
! convert DIC and alkalinity from model units to diagnostic output
! mg C/m3 --> umol/kg
! mmol eq/m3 --> umol/kg
         DIC = O3c/MW_C/ERHO*1000.0_rk
         ALK = O3h/ERHO*1000.0_rk  
         PHOS=N1p  !/ERHO*1000.0_rk
         SILIC=N5s ! /ERHO*1000.0_rk
         pH=8.1_rk ! first guess to be passed to CarbonateSystem
       
!         EPR=5.04097033087945_rk !! used to test bfm0d 

       error= CarbonateSystem( ESW, ETW,ERHO, &
               PHOS, SILIC, DIC, ALK, &
               CO2, HCO3, CO3, pH, &
               pCO2, patm=p_atm0, pr_in=EPR, &
               OmegaC=OCalc,OmegaA=OArag,fCO2=ffCO2)
! CO2,HCO3,CO3 in umol/kg
! pH on total scale
! pCO2 in uatm which is equivalent to ppm

         if (error > 0 ) then
              stop 'carbonate solver failed to converge' 
! if error copy the previous values
!            _GET_(self%id_CarbA_in,CO2)
!            _GET_(self%id_BiCarb_in,HCO3)
!            _GET_(self%id_Carb_in,CO3)
!            _GET_(self%id_pCO2_in,pCO2)
!            _GET_(self%id_pH_in,pH)
         endif
!            CO3 = CO3/1.e3_rk/density  ! from mmol/m3 to mol/kg
!            HCO3 = HCO3/1.e3_rk/density  ! from mmol/m3 to mol/kg
!            H2CO3 = H2CO3/1.e3_rk/density  ! from mmol/m3 to mol/kg
!            pCO2 = pCO2/1.e6_rk  ! from uatm to atm
     
         _SET_DIAGNOSTIC_(self%id_ph,pH)
         _SET_DIAGNOSTIC_(self%id_pco2,PCO2)  ! pCO2 in ppm
         _SET_DIAGNOSTIC_(self%id_CarbA, CO2)
         _SET_DIAGNOSTIC_(self%id_Bicarb,HCO3)
         _SET_DIAGNOSTIC_(self%id_Carb,  CO3)
         _SET_DIAGNOSTIC_(self%id_DIC,  DIC)
         _SET_DIAGNOSTIC_(self%id_ALK,  ALK)
         _SET_DIAGNOSTIC_(self%id_OCalc,OCalc)
         _SET_DIAGNOSTIC_(self%id_OArag,OArag)
         _SET_DIAGNOSTIC_(self%id_EPRdiag,EPR)
!         _SET_DIAGNOSTIC_(self%id_depth_diag,Depth)
 
     _LOOP_END_
   end subroutine



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Computes air-sea flux (only at surface points)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ogs_bfm_PelagicCSYS), intent(in) :: self

      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: O3c,ETW,ESW,ERHO,EWIND, EPR, EPRatm, Depth
      real(rk) :: O3h,N1p,N5s  
      real(rk) :: DIC,ALK,PHOS,SILIC

      real(rk) :: pH,PCO2,H2CO3,HCO3,CO3,CO2
      real(rk) :: OCalc,OArag,ffCO2

      real(rk) :: wind,PCO2A
      real(rk) :: fwind,UPTAKE,FAIRCO2
      integer  :: error
    ! local variable
    real(rk),parameter    :: PERMIL=ONE/1000.0_rk
    real(rk),parameter    :: bar2atm  = ONE/1.01325_rk ! Conversion factor from bar to atm
    real(rk),parameter    :: Rgas_atm = 82.05736_rk    ! (cm3 * atm) / (mol * K) CODATA (2006)
    real(rk),parameter    :: kwfac = 0.01_rk * HOURS_PER_DAY ! convert from cm/h to m/d
    real(rk)   :: pco2atm, patma, co2starair, co2star
    real(rk)   :: pH20, fco2atm, fugcoeff, bb, Del, xc2
    real(rk)   :: tk, tk100, tk1002, itk100
    real(rk)   :: pschmidt, Kw660, kwgas, K0
    real(rk)   :: rho,Fice
    real(rk)   :: co2ocn,xco2,temp,salt, AirSeaCO2
    ! CO2 solubility coefficients
    real(rk),parameter :: A(3) = (/-58.0931_rk, 90.5069_rk,22.2940_rk/)
    real(rk),parameter :: B(3) = (/0.027766_rk,-0.025888_rk,0.0050578_rk/)
    !! Schmidt coefficients
    real(rk),parameter :: Sc = 660.0_rk
    real(rk),parameter :: C(5)  = (/2116.8_rk, -136.25_rk, 4.7353_rk,-0.092307_rk, 0.0007555_rk/)
    ! This is for alternative fit of co2starair (see Orr et al, 2017 GMD)
    !real(RLEN)   :: phi0atm
    !real(RLEN),parameter :: F(7) = (/-160.7333, 215.4152, 89.8920, -1.47759,
    !0.029941, -0.027455, 0.0053407/)
    !---------------------------------------------------------------------------

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_O3c,O3c)
         _GET_(self%id_O3h,O3h)
         _GET_(self%id_N1p,N1p)
         _GET_(self%id_N5s,N5s)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_ESW,ESW)
         _GET_(self%id_ERHO,ERHO)
         _GET_(self%id_EPR,EPR)
         _GET_(self%id_depth,Depth)  ! depth water column cell in m
         _GET_HORIZONTAL_(self%id_EPRatm,EPRatm)
         _GET_HORIZONTAL_(self%id_EWIND,EWIND)
         _GET_HORIZONTAL_(self%id_PCO2A,PCO2A)
        
         wind = max(EWIND,0.0_rk)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Compute carbonate system equilibria
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! To use the Pressure correction of CSYS here the pr_in=EPS value
! convert DIC and alkalinity from model units to diagnostic output
! mg C/m3 --> umol/kg
! mmol eq/m3 --> umol/kg
         DIC = O3c/MW_C/ERHO*1000.0_rk
         ALK = O3h/ERHO*1000.0_rk
         PHOS=N1p  !/ERHO*1000.0_rk
         SILIC=N5s ! /ERHO*1000.0_rk
         pH=8.1_rk ! first guess to be passed to CarbonateSystem

!         EPR=5.04097033087945_rk !! used to test bfm0d 
       error= CarbonateSystem( ESW, ETW,ERHO, &
               PHOS, SILIC, DIC, ALK, &
               CO2,HCO3, CO3, pH, &
               pCO2, patm=p_atm0, pr_in=EPR, &
               OmegaC=OCalc,OmegaA=OArag,fCO2=ffCO2)
  
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   ! compute Air-Sea CO2 Exchange
! from AirSeaExchange.F90 routine of BFM
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  xco2   ! Atmospheric mixing ratio of Gas [ppmv]
! patm    ! Atmospheric pressure [mbar]
! temp    ! in-situ temperature at surface
! salt    ! practical salinity at surface
! rho     ! in-situ density 
! wind    ! wind speed [m/s]
! Fice    ! Sea-ice cover fraction [0-1]
! co2ocn  ! gas concentration at surface [umol/kg]

    co2ocn=CO2
    xco2=PCO2A ! PCO2A is a standard variable in ppm
    rho= ERHO
    Fice=0._rk ! no ice
    temp = ETW
    salt = ESW
    tk = temp - ZERO_KELVIN
    tk100 = tk/100.0_rk
    tk1002 = tk100*tk100
    itk100 = 1.0_rk/tk100
    patma = p_atm0 * PERMIL * bar2atm ! convert patm from mbar to atm
    !
    ! Oceanic concentration [CO2*] from umol/kg to mol/kg
    co2star = co2ocn * 1.0e-6_rk
    !
    ! Compute vapor pressure of seawater [in atm] (Weiss and Price, 1980, Eq. 10)
    pH20 = exp(24.4543_rk - 67.4509_rk*(100._rk/tk) - 4.8489_rk*log(tk/100.0_rk) - 0.000544_rk*salt)
    !
    ! Compute pco2atm [uatm] from xco2 [ppm], atmospheric pressure [atm], &
    ! vapor pressure of seawater pH20 [atm]
    pco2atm = (patma - pH20) * xco2
    ! 
    ! Compute fCO2atm [uatm] from pCO2atm [uatm] & fugacity coefficient
    ! [unitless]
    bb = -1636.75_rk + 12.0408_rk*tk - 0.0327957_rk*(tk*tk) +0.0000316528_rk*(tk*tk*tk)
    Del = 57.7_rk - 0.118_rk*tk
    xc2 = (ONE - (pco2atm*1.e-6_rk) )**2
    fugcoeff = EXP( patma*(bb + 2.0d0*xc2*Del)/(Rgas_atm*tk) )
    fco2atm = pco2atm * fugcoeff
    !
    ! Surface solubility of gas K0  [(mol/kg) / atm] at T, S of surface water
    ! according to Weiss (1974)
    K0 = exp( A(1) + A(2)/tk100 + A(3)*log(tk100) + salt*(B(1) + B(2)*tk100 + B(3)*tk1002) )
    ! Equilbrium [CO2*air] for atm gas at Patm & sfc-water T,S [mol/kg]
    co2starair = K0 * fco2atm * 1.0e-6_rk

    ! Alternative computation
    !! Solubility function for atmospheric CO2 saturation concentration (Orr et
    !al. GMD 2017, Eq. 15)
    !phi0atm = exp (F(1) + F(2)*itk100 + F(3)*log(tk100) + F(4)*tk1002 +
    !salt*(F(5) + F(6)*tk100 + F(7)*tk1002))
    !! Compute saturation concentration in atmosphere [mol/kg] at total pressure
    !Pa as in Orr et al (GMD 2017, Eq. 16)
    !co2starair = patma * phi0atm * xco2 * 1.0e-6_RLEN
    ! 
    ! Compute piston velolicty kw660 (at 25 C) from wind speed (Wanninkhof 2014,
    ! Limnol. Oceanograph. Methods, 12, 351-362)
    kw660 = 0.251_rk * wind**2 * kwfac
    !  Schmidt number (Sc), ratio between the kinematic viscosity and the
    !  molecular diffusivity of the gas
    pschmidt = C(1) + C(2)*temp + C(3)*temp**2 + C(4)*temp**3 + C(5)*temp**4
    ! Transfer velocity for gas in m/d (see equation [4] in OCMIP2 design
    ! document & OCMIP2 Abiotic HOWTO)
    kwgas = kw660 * (Sc/pschmidt)**0.5

    ! Air-sea gas flux [mmol/(m2 * day)]
    AirSeaCO2 = kwgas * (co2starair - co2star) * (ONE - Fice) * rho * 1000

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wnd_diag,wind) ! diagnostic in m/s
         _SET_SURFACE_EXCHANGE_(self%id_O3c,AirSeaCO2*12._rk) ! co2 flux in mgC/m2/d because O3c is in mgC/m3
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fair,AirSeaCO2)
      _HORIZONTAL_LOOP_END_
   end subroutine

! -------------------------------------------------------------------
! function CarbonateSystem
! -------------------------------------------------------------------
  function CarbonateSystem(salt,temp,rho,n1p,n5s,dic,alk,      &
                   co2,hco3,co3,pH,pco2,patm,pr_in,OmegaC,OmegaA,fCO2)
! !DESCRIPTION
! See module preamble
!
! USES:
!  use constants,  ONLY : ZERO_KELVIN,MW_C,MW_Ca,Rgas
!  use mem_Param,  ONLY : p_atm0

! !INPUT PARAMETERS:
  IMPLICIT NONE
  integer :: CarbonateSystem
  real(rk),intent(IN)            :: salt ! practical salinity
  real(rk),intent(IN)            :: temp ! in-situ temperature
  real(rk),intent(IN)            :: rho  ! in-situ density
  real(rk),intent(IN)            :: n1p
  real(rk),intent(IN)            :: n5s
  real(rk),intent(IN)            :: dic
  real(rk),intent(IN)            :: alk
  real(rk),intent(IN),optional   :: patm  ! Atmospheric pressure  [atm]
  real(rk),intent(IN),optional   :: pr_in ! Water column Pressure [dbar]
!
! !INPUT/OUTPUT PARAMETERS:
  real(rk),intent(INOUT)         :: pH
!
! !OUTPUT PARAMETERS:
!
  real(rk),intent(OUT)           :: co2
  real(rk),intent(OUT)           :: hco3
  real(rk),intent(OUT)           :: co3
  real(rk),intent(OUT)           :: pco2
  real(rk),intent(OUT),optional  :: OmegaC
  real(rk),intent(OUT),optional  :: OmegaA
  real(rk),intent(OUT),optional  :: fCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(rk),parameter    :: MEG=1.E6_rk, &
                             PERMIL=ONE/1000.0_rk,&
                             PERMEG=ONE/MEG
  real(rk)              :: press    = 0.0_rk         ! local hydrostatic pressure in bar
  real(rk)              :: pratm    = 0.0_rk         ! local atmospheric pressure in bar
  real(rk),parameter    :: Rgas_atm = 82.05736_rk    ! (cm3 * atm) / (mol * K)  CODATA (2006)
  real(rk),parameter    :: Rgas     = 83.14472_rk    ! Gas constant: L·mbar·K-1·mol-1
  real(rk),parameter    :: vbarCO2  = 32.3_rk        ! partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
  real(rk),parameter    :: CaRelCon = 0.02128_rk     ! Calcium ion relative concentration See Dickson (2007) and Munhoven (2013)
  real(rk),parameter    :: bar2atm  = ONE/1.01325_rk ! Conversion factor from bar to atm
  real(rk),parameter    :: MW_Ca=40.078_rk           ! Molecular weight Calcium

  real(rk)  ::  K0 ! solubility : [Co2]=k0Ac Pco2
  real(rk)  ::  K1 ! carbonate equilibrium I
  real(rk)  ::  K2 ! carbonate equilibrium II
  real(rk)  ::  Kw ! water dissociation
  real(rk)  ::  Kb ! constant for Boron equilibrium
  real(rk)  ::  Ks ! constant for bisulphate equilibrium
  real(rk)  ::  Kf ! constant for hidrogen fluoride equilibirum
  real(rk)  ::  K1p,K2p,K3p ! constants for phosphate equilibirum
  real(rk)  ::  Ksi ! constant for silicic acid equilibrium
  real(rk)  ::  Kspc ! constant for calcite equilibrium
  real(rk)  ::  Kspa ! constant for aragonite equilibrium

  integer     :: i, l, error
  real(rk)  :: H, Hi, ta, tc, sit, pt, tfco2
  real(rk)  :: Bt, Ft, St
  real(rk)  :: lnk, Ks_0p, Kf_0p, deltav, deltak
  real(rk)  :: tk, tk100, tk1002, temp2, invtk, dlogtk, is, is2, sqrtis
  real(rk)  :: s, s2, sqrts, s15, scl
  real(rk)  :: ptot, pr, pr2
  real(rk)  :: total2SWS, SWS2total, free2SWS, free2SWS_0p,total2SWS_0p, &
                  total2free, total2free_0p
  ! Carbonate Alkalinity 
  real(rk)  :: HSO4, HF, HSI, HPO4, ab, aw, ac, cu, cb, cc, Ca
  ! Fugacity
  real(rk)  :: B, Del, xc2, xCO2approx, fugcoeff
  ! Pressure correction
  real(rk), DIMENSION(12) :: a0, a1, a2, b0, b1, b2, lnkpok0
  ! 
  ! Constants in formulation for Pressure effect on K's (Millero, 95)
  ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
  !            9) K1P, 10) K2P, 11) K3P, 12) Ksi
  DATA a0 / 25.50_rk, 15.82_rk, 29.48_rk, 25.60_rk, 18.03_rk, 9.78_rk, &
            48.76_rk, 46.00_rk, 14.51_rk, 23.12_rk, 26.57_rk, 29.48_rk /
  DATA a1 / 0.1271_rk, -0.0219_rk, 0.1622_rk, 0.2324_rk, 0.0466_rk, -0.0090_rk, &
            0.5304_rk,  0.5304_rk, 0.1211_rk, 0.1758_rk, 0.2020_rk, 0.1622_rk /
  DATA a2 / 0.0_rk, 0.0_rk, 2.608_rk, -1.409_rk, 0.316_rk, -0.942_rk, &
            0.0_rk, 0.0_rk,-0.321_rk, -2.647_rk, -3.042_rk, -2.608_rk /
  DATA b0 /  3.08_rk, -1.13_rk, 2.84_rk, 5.13_rk, 4.53_rk, 3.91_rk, &
            11.76_rk, 11.76_rk, 2.67_rk, 5.15_rk, 4.08_rk, 2.84_rk /
  DATA b1 / 0.0877_rk, -0.1475_rk, 0.0_rk,  0.0794_rk, 0.09_rk, 0.054_rk, &
            0.3692_rk,  0.3692_rk, 0.0427_rk, 0.09_rk, 0.0714_rk, 0.0_rk  /
  !DATA b2 / 12*0.0_rk / ! not used as it is zero
  !
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
  tk     = temp - ZERO_KELVIN  !ZERO_KELVIN=-273.15; ETW is in degC; tk is in degK
  temp2  = temp * temp
  tk100  = tk / 100.0_rk
  tk1002 = tk100 * tk100
  invtk  = ONE / tk
  dlogtk = log(tk)
  
  ! Salinity and simply related values
  s  = salt
  s2 = salt*salt
  sqrts = sqrt(salt)
  s15 = salt**1.5_rk

  ! Hydrostatic Pressure [dbar]
  if (present(pr_in)) then
    press = pr_in * 0.1_rk  ! convert from dbar to bar
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
  scl = s/1.80655_rk

  ! ionic strength
  is = 19.924_rk*s/ (1000._rk-1.005_rk*s)
  is2 = is*is
  sqrtis = sqrt(is)

  ! Total concentrations for sulfate, fluoride, and boron
  ! Sulfate: Morris & Riley (1966)
  St = 0.14_rk * scl/96.062_rk
  ! Fluoride: Riley (1965)
  Ft = 0.000067_rk * scl/18.9984_rk
  ! Boron: Lee et al (2010)
  Bt = 0.0002414_rk * scl/10.811_rk

  ! -----------------------------------------------------------------------
  ! K0, solubility of co2 in the water (K Henry) from Weiss 1974
  ! K0 = [CO2]/ fCO2 [mol/kg/atm]
  ! -----------------------------------------------------------------------
  lnK = 93.4517_rk/tk100 - 60.2409_rk + 23.3585_rk * log(tk100) +   &
       s * (0.023517_rk - 0.023656_rk * tk100 + 0.0047036_rk * tk1002)
  K0 = exp ( lnK )

  ! -----------------------------------------------------------------------
  ! Choice of Acidity constants
  ! K1 = [H][HCO3]/[H2CO3]   ,   K2 = [H][CO3]/[HCO3]
  ! -----------------------------------------------------------------------
  ! Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total scale)
  K1 = 10.0_rk**(-ONE*(3633.86_rk*invtk - 61.2172_rk + 9.6777_rk*dlogtk - &
       0.011555_rk * s + 0.0001152_rk * s2))
  K2 = 10.0_rk**(-ONE*(471.78_rk*invtk + 25.9290_rk - &
       3.16967_rk*dlogtk - 0.01781_rk * s + 0.0001122_rk * s2))
  ! ! Millero (2010, Mar. Fresh Wat. Res.) (seawater scale)
  ! K1 = 10.0_rk**(-ONE*( (6320.813_rk*invtk + 19.568224_rk*dlogtk
  ! -126.34048_rk + &
  !      13.4038_rk*sqrts + 0.03206_rk*s - (5.242e-5_rk)*s2) + &
  !      (-530.659_rk*sqrts - 5.8210_rk*s)*invtk -2.0664_rk*sqrts*dlogtk)
  !      )
  ! K2 = 10.0_rk**(-ONE*( (5143.692_rk*invtk + 14.613358_rk*dlogtk
  ! -90.18333_rk + &
  !      21.3728_rk*sqrts + 0.1218_rk*s - (3.688e-4_rk)*s2 ) + &
  !      (-788.289_rk*sqrts - 19.189_rk*s)*invtk -3.374_rk*sqrts*dlogtk) )
  ! ! Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014) (seawater scale)
  ! K1 = 10.0_rk**(-ONE*( (6320.813_rk*invtk + 19.568224_rk*dlogtk
  ! -126.34048_rk + &
  !      13.409160_rk*sqrts + 0.031646_rk*s - (5.1895e-5_rk)*s2 ) + &
  !      (-531.3642_rk*sqrts - 5.713_rk*s)*invtk
  !      -2.0669166_rk*sqrts*dlogtk) )
  ! K2 = 10.0_rk**(-ONE*( (5143.692_rk*invtk + 14.613358_rk*dlogtk
  ! -90.18333_rk + & 
  !      21.225890_rk*sqrts + 0.12450870_rk*s - (3.7243e-4_rk)*s2 ) + &
  !      (-779.3444_rk*sqrts - 19.91739_rk*s)*invtk
  !      -3.3534679_rk*sqrts*dlogtk) )
  !-----------------------------------------------------------------------
  ! Kb = [H][BO2]/[HBO2] 
  ! Millero p.669 (1995) using data from Dickson (1990)    (total scale)
  !-----------------------------------------------------------------------
  lnK = (-8966.90_rk - 2890.53_rk*sqrts - 77.942_rk*s +    &
       1.728_rk*s15 - 0.0996_rk*s2)*invtk +                  &
       (148.0248_rk + 137.1942_rk*sqrts + 1.62142_rk*s) +  &
       (-24.4344_rk - 25.085_rk*sqrts - 0.2474_rk*s) *     &
       dlogtk + 0.053105_rk*sqrts*tk
  Kb = exp(lnK)

  ! -----------------------------------------------------------------------
  ! K1p = [H][H2PO4]/[H3PO4]
  ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
  ! Millero (1995), p.670, eq. 65                        (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = -4576.752_rk*invtk + 115.540_rk - 18.453_rk * dlogtk + &
       (-106.736_rk*invtk + 0.69171_rk) * sqrts +                &
       (-0.65643_rk*invtk - 0.01844_rk) * s
  K1p = exp(lnK)

  ! -----------------------------------------------------------------------
  ! K2p = [H][HPO4]/[H2PO4]
  ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
  ! Millero (1995), p.670, eq. 66                        (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = -8814.715_rk*invtk + 172.1033_rk - 27.927_rk * dlogtk +  &
       (-160.340_rk*invtk + 1.3566_rk) * sqrts +                   &
       (0.37335_rk*invtk - 0.05778_rk) * s
  K2p = exp(lnK)

  !------------------------------------------------------------------------
  ! K3p = [H][PO4]/[HPO4]
  ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  ! Millero (1995), p.670, eq. 67                        (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = -3070.75_rk*invtk - 18.126_rk +     &
       (17.27039_rk*invtk + 2.81197_rk) *   &
       sqrts + (-44.99486_rk*invtk - 0.09984_rk) * s
  K3p = exp(lnK)

  !------------------------------------------------------------------------
  ! ksi = [H][SiO(OH)3]/[Si(OH)4]
  ! Millero (1995), p.671, eq. 72                        (seawater scale)
  ! -----------------------------------------------------------------------
  ksi = -8904.2_rk*invtk + 117.400_rk - 19.334_rk * dlogtk + &
       (-458.79_rk*invtk + 3.5913_rk) * sqrtis +               &
       (188.74_rk*invtk - 1.5998_rk) * is +                    &
       (-12.1652_rk*invtk + 0.07871_rk) * is2 +                &
       log(ONE-0.001005_rk*s)
  Ksi = exp(lnK)

  ! -----------------------------------------------------------------------
  ! Kw = [H][OH]
  ! Millero (1995) p.670, eq. 63 from composite data     (seawater scale)
  ! -----------------------------------------------------------------------
  lnK = 148.9802_rk -13847.26_rk*invtk - 23.6521_rk * dlogtk +      &
       (118.67_rk*invtk - 5.977_rk + 1.0495_rk * dlogtk) *          &
       sqrts - 0.01615_rk * s
  Kw = exp(lnK)

  !------------------------------------------------------------------------
  ! ks = [H][SO4]/[HSO4]
  ! Dickson (1990, J. chem. Thermodynamics 22, 113)          (free scale)
  !------------------------------------------------------------------------
  lnK = -4276.1_rk*invtk + 141.328_rk - 23.093_rk*dlogtk +          &
       (-13856._rk*invtk + 324.57_rk - 47.986_rk*dlogtk) * sqrtis + &
       (35474._rk*invtk - 771.54_rk + 114.723_rk*dlogtk) * is -     &
       2698._rk*invtk*is**1.5_rk + 1776._rk*invtk*is2 +             &
       log(ONE - 0.001005_rk*s)
  Ks_0p = exp(lnK)

  !------------------------------------------------------------------------
  ! kf = [H][F]/[HF]
  ! Perez & Fraga (1987) recom. by Dickson et al., (2007)   (total scale)
  !------------------------------------------------------------------------
  lnK = 874.0_rk*invtk - 9.68_rk + 0.111_rk*sqrts
  Kf_0p = exp(lnK)

  !------------------------------------------------------------------------
  ! Kspc = [Ca2+] [CO32-] - apparent solubility product of Calcite
  ! Mucci (1983)  [mol/kg-soln]
  !------------------------------------------------------------------------
  kspc = 10.0_rk**(ONE*( -171.9065_rk - 0.077993_rk * tk + 2839.319_rk *invtk &
          + 71.595_rk * log10(tk) + sqrts * (-0.77712_rk +  &
          0.0028426_rk * tk + 178.34_rk * invtk)            &
          - 0.07711_rk * s + 0.0041249_rk *s15 ) )
  !------------------------------------------------------------------------
  ! Kspa = [Ca2+] [CO32-] - apparent solubility product of Aragonite
  ! Mucci (1983)  [mol/kg-soln]
  !------------------------------------------------------------------------
  kspa = 10.0_rk**(ONE*(  -171.945_rk - 0.077993_rk * tk + 2903.293_rk * invtk   &
          + 71.595_rk * log10(tk) + sqrts * (-0.068393_rk +   &
          0.0017276_rk * tk + 88.135_rk * invtk)   &
          - 0.10018_rk * s + 0.0059415_rk * s15 ) )

  ! Pressure effect on K0 based on Weiss (1974, equation 5)
  K0 = K0 * exp( ((1-Ptot)*vbarCO2)/(Rgas_atm*tk) )

  lnkpok0 = 0.0_rk
  if (present(pr_in) .AND. press .GT. 0) then
  ! Pressure effect on all other K's (based on Millero, (1995)
  ! Index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
  !           K1p(9), K2p(10), K3p(11), Ksi(12)
     do l = 1, 12
        deltav  =  -a0(l) + a1(l)*temp + 1.0e-3_rk*a2(l)*temp2
        deltak  = 1.0e-3_rk*(-b0(l) + b1(l)*temp)
        lnkpok0(l) = -deltav*invtk*pr + 0.5_rk*deltak*invtk*pr2
     end do
  endif

  ! Pressure correction on Ks (Free scale)
  Ks = Ks_0p*EXP(lnkpok0(5))
  ! Conversion factor total -> free scale
  total2free     = 1.0_rk/(1.0_rk + St/Ks)   ! Kfree = Ktotal*total2free
  ! Conversion factor total -> free scale at pressure zero
  total2free_0p  = 1.0_rk/(1.0_rk + St/Ks_0p)   ! Kfree = Ktotal*total2free

  ! Pressure correction on Kf
  ! Kf must be on FREE scale before correction
  Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
  Kf    = Kf_0p * EXP(lnkpok0(6)) !Pressure correction (on Free scale)
  Kf    = Kf / total2free         !Convert back from Free to Total scale

  ! Convert between seawater and total hydrogen (pH) scales
  free2SWS  = 1.0_rk + St/Ks + Ft/(Kf*total2free)  ! using Kf on free scale
  total2SWS = total2free * free2SWS                  ! KSWS = Ktotal*total2SWS
  SWS2total = 1.0_rk / total2SWS
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
  if ( pH < 0._rk ) then
     Hi = Hini_for_at(ta,tc,bt,K1,K2,Kb)
  else
     Hi = 10.0_rk**(-pH)
  endif
 
  ! Solve for H+ using above result as the initial H+ value (mol/kg)
  H = solve_at_general(ta, tc, Bt, pt, sit, St, Ft,            &
               K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi, Hi )
  if (H < ZERO) CarbonateSystem = 1

  ! Calculate pH from from H+ concentration (mol/kg)
  pH = -ONE * LOG10( H )

  ! Compute carbonate Alk (Ac) by difference: from total Alk and other Alk
  ! components
  HSO4 = St/(1.0_rk + Ks/(H/(1.0_rk + St/Ks)))
  HF   = 1.0_rk / (1.0_rk + Kf/H)
  HSI  = 1.0_rk / (1.0_rk + H/Ksi)
  HPO4 = (K1p*K2p*(H + 2.*K3p) - H**3)    /         &
         (H**3 + K1p*H**2 + K1p*K2p*H + K1p*K2p*K3p)
  ab = Bt/(1.0_rk + H/Kb)
  aw = Kw/H - H/(1.0_rk + St/Ks)
  ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4

  ! Calculate CO2*, HCO3-, & CO32- (in mol/kg soln) from Ct, Ac, H+, K1, & K2
  cu = (2.0_rk * tc - ac) / (2.0_rk + K1 / H)
  cb = K1 * cu / H
  cc = K2 * cb / H

  ! Determine Omega Calcite and Aragonite (see Munhoven (2013, GMD))
  Ca = (CaRelCon / MW_Ca) * s/1.80655_rk
  if (present(OmegaA)) OmegaA = (Ca*cc) / Kspa
  if (present(OmegaC)) OmegaC = (Ca*cc) / Kspc

  ! Determine CO2 fugacity [uatm]
  tfco2 = cu * 1.e6_rk/K0
  if (present(fCO2)) fCO2 = tfco2

  ! Determine CO2 partial pressure from CO2 fugacity [uatm]
  ! compute fugacity coefficient terms : B, Del, xc2
  B = -1636.75_rk + 12.0408_rk*tk - 0.0327957_rk*(tk*tk) + 0.0000316528_rk*(tk*tk*tk)
  Del = 57.7_rk - 0.118_rk*tk
  xCO2approx = tfco2 * 1.e-6_rk

  if (present(pr_in) .AND. press .GT. 0) &
     xCO2approx = xCO2approx * exp( ((ONE-Ptot)*32.3_rk)/(82.05736_rk*tk) ) ! of K0 press. correction, see Weiss (1974, equation 5)
  xc2 = (ONE - xCO2approx)**2
  fugcoeff = exp( Ptot*(B + 2.0_rk*xc2*Del)/(Rgas_atm*tk) )
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

  ! Returns * 1E-03_rk if p_alkcb <= 0
  !         * 1E-10_rk if p_alkcb >= 2*p_dictot + p_bortot
  !         * 1E-07_rk if 0 < p_alkcb < 2*p_dictot + p_bortot
  !                    and the 2nd order approximation does not have a solution 
    IMPLICIT NONE
    REAL(KIND=rk) :: Hini_for_at   ! [H+] in mol/kg

    ! Argument variables
    !--------------------
    REAL(KIND=rk), INTENT(IN)   ::  p_alkcb, p_dictot, p_bortot ! [mol/kg]
    REAL(KIND=rk), INTENT(IN)   ::  K1, K2, Kb

    ! Local variables
    !-----------------
    REAL(KIND=rk)  ::  p_hini
    REAL(KIND=rk)  ::  zca, zba
    REAL(KIND=rk)  ::  zd, zsqrtd, zhmin
    REAL(KIND=rk)  ::  za2, za1, za0

    IF (p_alkcb <= 0._rk) THEN
      p_hini = 1.e-3_rk
    ELSEIF (p_alkcb >= (2._rk*p_dictot + p_bortot)) THEN
      p_hini = 1.e-10_rk
    ELSE
      zca = p_dictot/p_alkcb
      zba = p_bortot/p_alkcb

      ! Coefficients of the cubic polynomial
      za2 = Kb*(1._rk - zba) + K1*(1._rk-zca)
      za1 = K1*Kb*(1._rk - zba - zca) + K1*K2*(1._rk - (zca+zca))
      za0 = K1*K2*Kb*(1._rk - zba - (zca+zca))
                                            ! Taylor expansion around the
                                            ! minimum
      zd = za2*za2 - 3._rk*za1              ! Discriminant of the quadratic equation
                                            ! for the minimum close to the root

      IF(zd > 0._rk) THEN                   ! If the discriminant is positive
        zsqrtd = SQRT(zd)
        IF(za2 < 0._rk) THEN
          zhmin = (-za2 + zsqrtd)/3._rk
        ELSE
          zhmin = -za1/(za2 + zsqrtd)
        ENDIF
        p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
      ELSE
        p_hini = 1.e-7_rk
      ENDIF

    Hini_for_at = p_hini

    ENDIF
    RETURN
  END FUNCTION hini_for_at

! ----------------

  elemental FUNCTION solve_at_general(p_alktot, p_dictot, p_bortot, &
                              p_po4tot, p_siltot, &
                              p_so4tot, p_flutot, &
                              K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi, &
                              p_h)

    ! Purpose: Compute [H°] ion concentration from sea-water ion
    ! concentrations,
    !          alkalinity, DIC, and equilibrium constants
    ! Universal pH solver that converges from any given initial value,
    ! determines upper an lower bounds for the solution if required

    IMPLICIT NONE
    REAL(KIND=rk) :: SOLVE_AT_GENERAL    ! [H+] in mol/kg

    ! Argument variables
    !--------------------
    REAL(KIND=rk), INTENT(IN)            :: p_alktot  ! [mol/kg] 
    REAL(KIND=rk), INTENT(IN)            :: p_dictot  ! [mol/kg]
    REAL(KIND=rk), INTENT(IN)            :: p_bortot  ! [mol/kg]
    REAL(KIND=rk), INTENT(IN)            :: p_po4tot  ! [mol/kg]
    REAL(KIND=rk), INTENT(IN)            :: p_siltot  ! [mol/kg]
    REAL(KIND=rk), INTENT(IN)            :: p_so4tot  ! [mol/kg]
    REAL(KIND=rk), INTENT(IN)            :: p_flutot  ! [mol/kg]
    REAL(KIND=rk), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
    REAL(KIND=rk), INTENT(IN)            :: K1p, K2p, K3p, Ksi
    REAL(KIND=rk), INTENT(IN)            :: p_h       ! [H+] in mol/kg

    ! Local variables
    !-----------------
    REAL(KIND=rk)  ::  zh, zh_prev, zh_lnfactor
    REAL(KIND=rk)  ::  zalknw_inf, zalknw_sup
    REAL(KIND=rk)  ::  zh_min, zh_max
    REAL(KIND=rk)  ::  zdelta, zh_delta
    REAL(KIND=rk)  ::  zeqn, zdeqndh, zeqn_absmin
    REAL(KIND=rk)  ::  aphscale
    LOGICAL          ::  l_exitnow
    INTEGER          ::  niter_atgen
    REAL(KIND=rk)  ::  znumer, zdenom, zdnumer
    REAL(KIND=rk)  ::  zalk_dic, zalk_bor, zalk_po4, zalk_sil, zalk_so4,zalk_flu, zalk_wat
    REAL(KIND=rk)  ::  zdalk_dic, zdalk_bor, zdalk_po4, zdalk_sil, zdalk_so4,zdalk_flu, zdalk_wat

    ! General parameters
  !-----------------
    REAL(KIND=rk), PARAMETER :: pz_exp_threshold = 1.0_rk
    REAL(KIND=rk), PARAMETER :: pp_rdel_ah_target = 1.E-8_rk

  ! MaxIterPHsolver integer          Maximum number of iterations (default 50)
    integer,parameter              :: MaxIterPHsolver = 50 ! in BFM it was in NAMELIST PelCO2_parameters
    
    ! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
    aphscale = 1._rk + p_so4tot/Ks

    ! lower and upper bounds of "non-water-selfionization"
    ! contributions to total alkalinity (the infimum and the supremum), i.e
    ! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])
    zalknw_inf = -p_po4tot - p_so4tot - p_flutot
    zalknw_sup = p_dictot + p_dictot + p_bortot + p_po4tot + p_po4tot + p_siltot

    zdelta = (p_alktot-zalknw_inf)**2 + 4._rk*Kw/aphscale

    IF(p_alktot >= zalknw_inf) THEN
       zh_min = 2._rk*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
    ELSE
       zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._rk
    ENDIF

    zdelta = (p_alktot-zalknw_sup)**2 + 4._rk*Kw/aphscale

    IF(p_alktot <= zalknw_sup) THEN
       zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._rk
    ELSE
       zh_max = 2._rk*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
    ENDIF

    zh = MAX(MIN(zh_max, p_h), zh_min)
    niter_atgen        = 0                 ! Reset counters of iterations
    zeqn_absmin        = HUGE(1._rk)

    DO
       IF(niter_atgen >= MaxIterPHsolver ) THEN
          zh = -1._rk
          EXIT
       ENDIF

       zh_prev = zh

       ! Compute total alkalinity from ion concentrations and equilibrium
       ! constants
       ! H2CO3 - HCO3 - CO3 : n=2, m=0
       znumer     = 2._rk*K1*K2 + zh*       K1
       zdenom     =         K1*K2 + zh*(      K1 + zh)
       zalk_dic   = p_dictot * (znumer/zdenom)
       zdnumer    = K1*K1*K2 + zh*(4._rk*K1*K2 + zh* K1 )
       zdalk_dic  = -p_dictot * (zdnumer/zdenom**2)
       ! B(OH)3 - B(OH)4 : n=1, m=0
       znumer     = Kb
       zdenom     = Kb + zh
       zalk_bor   = p_bortot * (znumer/zdenom)
       zdnumer    = Kb
       zdalk_bor  = -p_bortot * (zdnumer/zdenom**2)
       ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
       znumer     = 3._rk*K1p*K2p*K3p + zh*(2._rk*K1p*K2p + zh* K1p)
       zdenom     =         K1p*K2p*K3p + zh*(        K1p*K2p + zh*(K1p + zh))
       zalk_po4   = p_po4tot * (znumer/zdenom - 1._rk) ! Zero level of H3PO4 = 1
       zdnumer    = K1p*K2p*K1p*K2p*K3p + zh*(4._rk*K1p*K1p*K2p*K3p &
                                        + zh*(9._rk*K1p*K2p*K3p + K1p*K1p*K2p &
                                        + zh*(4._rk*K1p*K2p &
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
       zalk_so4   = p_so4tot * (znumer/zdenom - 1._rk)
       zdnumer    = Ks
       zdalk_so4  = -p_so4tot * (zdnumer/zdenom**2)
       ! HF - F : n=1, m=1
       znumer     = Kf
       zdenom     = Kf + zh
       zalk_flu   = p_flutot * (znumer/zdenom - 1._rk)
       zdnumer    = Kf
       zdalk_flu  = -p_flutot * (zdnumer/zdenom**2)
       ! H2O - OH
       zalk_wat   = Kw/zh - zh/aphscale
       zdalk_wat  = -Kw/zh**2 - 1._rk/aphscale

       zeqn = zalk_dic + zalk_bor + zalk_po4 + zalk_sil   &
               + zalk_so4 + zalk_flu + zalk_wat - p_alktot

       zdeqndh = zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
               + zdalk_so4 + zdalk_flu + zdalk_wat

       ! Adapt bracketing interval
       IF(zeqn > 0._rk) THEN
          zh_min = zh_prev
       ELSEIF(zeqn < 0._rk) THEN
          zh_max = zh_prev
       ELSE
          ! zh is the root; unlikely but, one never knows
          EXIT
       ENDIF
       ! Now determine the next iterate zh
       niter_atgen = niter_atgen + 1

       IF(ABS(zeqn) >= 0.5_rk*zeqn_absmin) THEN
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



end module
 
!GP !EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



































