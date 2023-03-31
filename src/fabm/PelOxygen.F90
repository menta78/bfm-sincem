#include "fabm_driver.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
  module bfm_PelOxygen 
! !ROUTINE: PelOxygen
!
! DESCRIPTION
!   ! compute Oxygen dynamics in the Pelagic environment
!
! !AUTHORS
!   T. Lovato (CMCC) 2017
!
! !REVISION_HISTORY
!
! COPYING
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
!EOP
!-------------------------------------------------------------------------!
!BOC

   use fabm_types
   use ogs_bfm_shared
!   use ogs_bfm_pelagic_base

   implicit none

   private
   
   type,extends(type_base_model),public :: type_ogs_bfm_PelOxygen
!     Variable identifiers
      type (type_state_variable_id)     :: id_O2o
      type (type_dependency_id)         :: id_ETW, id_ESW
      type (type_horizontal_dependency_id) :: id_EWIND

      type (type_diagnostic_variable_id) :: id_eO2mO2,id_osat,id_aou
      type (type_horizontal_diagnostic_variable_id) ::  id_o2flux
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type type_ogs_bfm_PelOxygen
   
   
contains

   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_PelOxygen), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

! Initialize pelagic base model (this also sets the time unit to per day,instead
! of the default per second)
!      call self%initialize_bfm_base()

 
      call self%register_state_variable(self%id_O2o,'o','mmol O_2/m^3','oxygen',300._rk)

      call self%register_diagnostic_variable(self%id_eO2mO2,'eO2mO2','1','relative saturation', &  
         standard_variable=standard_variables%fractional_saturation_of_oxygen,output=output_none)
      call self%register_diagnostic_variable(self%id_osat,'osat','mmol O_2/m^3','saturation concentration',output=output_none)
      call self%register_diagnostic_variable(self%id_aou,'AOU','mmol O_2/m^3','apparent utilisation',output=output_none)
      call self%register_diagnostic_variable(self%id_o2flux,'o2flux','mmolO_2/m^2/d','air-sea flux', source=source_do_surface,output=output_none)



      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_ESW,standard_variables%practical_salinity)
      call self%register_dependency(self%id_EWIND,standard_variables%wind_speed)
!      call self%register_dependency(self%id_EFICE,standard_variables%fraction_ice)

      ! Set time unit to d-1
      ! This implies that all rates (sink/source terms, vertical velocities) are
      ! given in d-1.
      self%dt = 86400._rk

   end subroutine



   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ogs_bfm_PelOxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: O2o,ETW,ESW,cxoO2
      
      _LOOP_BEGIN_
         _GET_(self%id_O2o,O2o)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_ESW,ESW)
         cxoO2 = oxygen_saturation_concentration(self,ETW,ESW)
         _SET_DIAGNOSTIC_(self%id_osat,cxoO2)
         _SET_DIAGNOSTIC_(self%id_eO2mO2,max(0.0_rk,O2o/cxoO2))
         _SET_DIAGNOSTIC_(self%id_aou,cxoO2-O2o)
      _LOOP_END_
   end subroutine



   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ogs_bfm_PelOxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: O2o,ETW,ESW,EWIND
      real(rk) :: wind, temp, Fice
      real(rk) :: kw660, tsrf, tsrf2, pschmidt, kwgas 
      real(rk) :: cxoO2,o2flux
      
    real(rk),parameter    :: kwfac  = 0.01_rk * 24.0_rk ! convert from cm/h to m/d
      
    !! Schmidt coefficients
      real(rk),parameter :: Sc = 660.0_rk
      real(rk),parameter :: C1=1920.4_rk
      real(rk),parameter :: C2=-135.6_rk
      real(rk),parameter :: C3=5.2122_rk
      real(rk),parameter :: C4=-0.10939_rk
      real(rk),parameter :: C5=0.00093777_rk

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_O2o,O2o)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_ESW,ESW)
         _GET_HORIZONTAL_(self%id_EWIND,EWIND)

         EWIND = max(EWIND, 0.0_rk)

         cxoO2 = oxygen_saturation_concentration(self,ETW,ESW)
         
    wind = EWIND
    temp = ETW
    Fice = 0.0
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Air-Sea Flux of Oxygen
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Compute piston velolicty kw660 (at 25 C) from wind speed (Wanninkhof 2014, Limnol. Oceanograph. Methods, 12, 351-362)
    kw660 = 0.251_rk * wind**2 * kwfac
    !  Schmidt number (Sc), ratio between the kinematic viscosity and the molecular diffusivity of the gas
    tsrf = temp
    tsrf2 = tsrf * tsrf 
    pschmidt = C1 + C2*tsrf + C3*tsrf2 + C4*tsrf*tsrf2 + C5*tsrf2*tsrf2
    ! Transfer velocity for gas in m/d (see equation [4] in OCMIP2 design document & OCMIP2 Abiotic HOWTO)
    kwgas = kw660 * (Sc/pschmidt)**0.5
    !
    o2flux = kwgas * ( cxoO2 - O2o ) * (1.0_rk - Fice)
  
        _SET_SURFACE_EXCHANGE_(self%id_O2o,o2flux)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_o2flux,o2flux)
      _HORIZONTAL_LOOP_END_
   end subroutine

   function oxygen_saturation_concentration(self,ETW,ESW) result(cxoO2)
      class (type_ogs_bfm_PelOxygen), intent(in) :: self
      real(rk),                      intent(in) :: ETW,ESW
      real(rk)                                  :: cxoO2
    
      real(rk),parameter :: A1= 2.00856_rk
      real(rk),parameter :: A2= 3.22400_rk
      real(rk),parameter :: A3= 3.99063_rk
      real(rk),parameter :: A4= 4.80299_rk
      real(rk),parameter :: A5= 9.78188E-01_rk
      real(rk),parameter :: A6= 1.71069_rk
      real(rk),parameter :: B1=-6.24097E-03_rk
      real(rk),parameter :: B2=-6.93498E-03_rk
      real(rk),parameter :: B3=-6.90358E-03_rk
      real(rk),parameter :: B4=-4.29155E-03_rk
      
      real(rk),parameter :: C0 = -3.11680e-7_rk
             
             
    real(rk),parameter    :: o2mvol = 22.391_rk               ! volume of 1 mole O2 at STP (ICES conversions)         
    real(rk),parameter    :: o2fac  = 1000._rk / o2mvol       ! convert from ml/l to umol/L [=mmol/m3]
     
    real(rk) :: temp, salt, tmpflux, ts, ts2, lnk
    real(rk) :: wind, fice, pschmidt, Kw660, kwgas, tsrf, tsrf2, o2flux, po2air, po2sea           
    !---------------------------------------------------------------------------
    ! input arrays
    temp = ETW
    salt = ESW

    ! common arrays
    tmpflux = 0.0_rk

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Oxygen saturation (Gordon and Garcia, 1992 L&O)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Scaled temperature 
    ts = log ( ( 298.15_rk - temp) / ( temp - (-273.15_rk) ) )
    ts2 = ts * ts
    lnk = A1 + A2*ts + A3*ts2 + A4*ts*ts2 + A5*ts2*ts2 + A6*ts2*ts2*ts + &
          salt * (B1 + B2*ts + B3*ts2 + B4*ts*ts2) + C0*salt*salt
    ! Compute O2sat in (mmol/m3) 
    cxoO2  = exp(lnk) * o2fac
    ! Relative oxygen saturation
    ! eO2mO2 = max( p_small,O2o ) / cxoO2
  end function

end module
 
!GP !EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
