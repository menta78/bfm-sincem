!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! MODULE: sw_tool
!
! DESCRIPTION
!   SeaWater tools for conversion of ocean physical quantities
!   Uses functions from UNESCO 1983 (EOS80) and TEOS10
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
   module sw_tool
!
! USES
   use global_mem, only: RLEN,ZERO,bfm_lwp,LOGUNIT,NMLUNIT,bfm_file_FirstUnit
   use mem,        only: NO_BOXES
   use global_mem, only: PI

   implicit none

   ! Multiplicative factor to convert the potential temperature from ITS90
   ! to IPTS-68 scale (TEOS-10, t90_from_t68)
   real (RLEN), parameter :: tem9068 = 1.00024_RLEN

!  TEOS-10 constants
   real (RLEN), parameter :: db2pa = 1.0e4_RLEN
   real (RLEN), parameter :: deg2rad = pi/180.0_RLEN
   real (RLEN), parameter :: gamma = 2.26e-7_RLEN

! !PUBLIC MEMBER FUNCTIONS:
   public sw_t_from_pt, sw_rho_t
   public gsw_p_from_z

   contains
   
!  EOS-80 Functions

   elemental function sw_t_from_pt (sp, pt, p, p_ref )
   !==========================================================================
   ! Calculates the in-situ temperature from potential temperature
   ! sp           : practical salinity                       [psu  PSS-78]
   ! pt           : potential temperature                    [degC ITS-90]
   ! p            : pressure                                 [dbar]
   ! p_ref        : reference pressure                       [dbar]
   ! sw_t_from_pt :                                          [degC ITS-90]
   ! Carry out inverse calculation by swapping P_ref (p_ref) and Pressure (p)
   ! in routine that is used to compute potential temp from in-situ temp
   !==========================================================================
   implicit none 

   real (RLEN), intent(in) :: sp, pt, p, p_ref
   real (RLEN) :: sw_t_from_pt

   sw_t_from_pt = sw_pt_from_t(sp, pt, p_ref, p)

   return
   end function

   !-------------------------------------------------------------------------- 

   elemental function sw_pt_from_t (sp, t, p, p_ref )
   !==========================================================================
   ! Calculates the potential temperature from in-situ temperature
   ! sp           : practical salinity                       [psu  PSS-78]
   ! t            : in-situ temperature                      [degC ITS-90]
   ! p            : pressure                                 [dbar]
   ! p_ref        : reference pressure                       [dbar]
   ! sw_pt_from_t :                                          [degC ITS-90]
   ! Potential temperature from the in-situ one as per UNESCO 1983 routines.
   ! Original function from Seawater 3.3.1 (sw_ptmp.m)
   !==========================================================================
   implicit none

   real (RLEN), intent(in) :: sp, t, p, p_ref
   real (RLEN) :: sw_pt_from_t
   real (RLEN) :: t68, del_P ,del_th, th, q
   real (RLEN), parameter :: sqrttwo  = sqrt(2.0_RLEN)
   real (RLEN), parameter :: item9068 = 1._RLEN / tem9068

   ! theta1
   del_P  = p_ref - p
   del_th = del_P * sw_adtg(sp,t,p)
   th     = (t * tem9068) + 0.5_RLEN * del_th
   q      = del_th
   ! theta2
   del_th = del_P * sw_adtg(sp,th*item9068,p+0.5_RLEN*del_P)
   th     = th + (1.0_RLEN - 1.0_RLEN/ sqrttwo ) * (del_th - q)
   q      = (2.0_RLEN-sqrttwo)*del_th + (-2.0_RLEN+3.0_RLEN/sqrttwo)*q
   ! theta3
   del_th = del_P * sw_adtg(sp,th*item9068,p+0.5_RLEN*del_P)
   th     = th + (1.0_RLEN + 1.0_RLEN/sqrttwo) * (del_th - q)
   q      = (2.0_RLEN + sqrttwo)*del_th + (-2.0_RLEN -3.0_RLEN/sqrttwo)*q
   ! theta4
   del_th = del_P * sw_adtg(sp,th*item9068,p+del_P)
   sw_pt_from_t = (th + (del_th - 2.0_RLEN*q) / (6.0_RLEN) ) * item9068

   return
   end function

   !--------------------------------------------------------------------------
   
   elemental function sw_adtg (sp, t, p)
   !==========================================================================
   ! Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
   ! sp            : practical salinity                       [psu  PSS-78]
   ! t             : in-situ temperature                      [degC ITS-90]
   ! p             : pressure                                 [dbar]
   ! sw_adtg       :                                          [degC / dbar]
   ! Original function from Seawater 3.3.1 (sw_adtg.m)
   !==========================================================================
   implicit none

   real (RLEN), intent(in) :: sp, t, p
   real (RLEN) :: sw_adtg
   real (RLEN) :: t68, aa, bb, cc, dd, ee

   real (RLEN), parameter :: sref = 35.00000e+00_RLEN
   real (RLEN), parameter :: a0   =  3.58030e-05_RLEN
   real (RLEN), parameter :: a1   =  8.52580e-06_RLEN
   real (RLEN), parameter :: a2   = -6.83600e-08_RLEN
   real (RLEN), parameter :: a3   =  6.62280e-10_RLEN
   real (RLEN), parameter :: b0   =  1.89320e-06_RLEN
   real (RLEN), parameter :: b1   = -4.23930e-08_RLEN
   real (RLEN), parameter :: c0   =  1.87410e-08_RLEN
   real (RLEN), parameter :: c1   = -6.77950e-10_RLEN
   real (RLEN), parameter :: c2   =  8.73300e-12_RLEN
   real (RLEN), parameter :: c3   = -5.44810e-14_RLEN
   real (RLEN), parameter :: d0   = -1.13510e-10_RLEN
   real (RLEN), parameter :: d1   =  2.77590e-12_RLEN
   real (RLEN), parameter :: e0   = -4.62060e-13_RLEN
   real (RLEN), parameter :: e1   =  1.86760e-14_RLEN
   real (RLEN), parameter :: e2   = -2.16870e-16_RLEN
   
   t68 = t * tem9068

   aa = a0 + ( a1 + (a2 + a3*t68)*t68 ) *t68
   bb = ( b0 + b1*t68 ) * ( sp - sref )
   cc = c0 + ( c1 + (c2 + c3*t68)*t68 ) *t68
   dd = (d0 + d1*t68) * ( sp - sref )
   ee = ( e0 + (e1 + e2*t68)*t68 ) *p *p

   sw_adtg = aa + bb + ( cc + dd )*p + ee

   return
   end function

   !--------------------------------------------------------------------------

   elemental function sw_rho_t (sp, pt, p)
   !==========================================================================
   ! Calculates in-situ density of seawater from PT and SP
   ! sp           : practical salinity                       [psu  PSS-78]
   ! pt           : potential temperature                    [degC ITS-90]
   ! p            : pressure                                 [dbar]
   ! sw_rho_t     :                                          [kg/m^3]
   ! This uses EOS-80 formula, where PT is IPTS-68 scale and pressue in bar
   ! Conversions are applied at the beginning
   !==========================================================================
   implicit none

   real (RLEN), intent(in) :: sp, pt, p
   real (RLEN) :: sw_rho_t
   real (RLEN) :: X, pr, rhow, rho0, drho
   real (RLEN) :: a, b, c
   real (RLEN) :: Ksbmw, Ksbm0, Ksbm 
     
   ! convert potential temperature to IPTS-68 scale (TEOS-10, t90_from_t68)
   X = 1.00024_RLEN * pt

   ! convert pressure from dbar to bar
   pr = p * 0.1_RLEN
   
   ! Density of pure water
   rhow = 999.842594e+00_RLEN + 6.793952e-02_RLEN*X          &
        &  -9.095290e-03_RLEN*X*X + 1.001685e-04_RLEN*X**3   &
        &  -1.120083e-06_RLEN*X**4 + 6.536332e-09_RLEN*X**5
   
   ! Density of seawater at 1 atm, P=0
   A = 8.24493e-01_RLEN - 4.0899e-03_RLEN*X                  &
     &  + 7.6438e-05_RLEN*X*X - 8.2467e-07_RLEN*X**3 + 5.3875e-09_RLEN*X**4
   B = -5.72466e-03_RLEN + 1.0227e-04_RLEN*X - 1.6546e-06_RLEN*X*X
   C =  4.8314e-04_RLEN
   rho0 = rhow + A*sp + B*sp*SQRT(sp) + C*sp**2.0e+00_RLEN
   
   ! Secant bulk modulus of pure water
   ! The secant bulk modulus is the average change in pressure
   ! divided by the total change in volume per unit of initial volume.
   Ksbmw = 19652.21e+00_RLEN + 148.4206e+00_RLEN*X - 2.327105e+00_RLEN*X*X    &
         &  + 1.360477e-02_RLEN*X**3 - 5.155288e-05_RLEN*X**4
   
   ! Secant bulk modulus of seawater at 1 atm
   Ksbm0 = Ksbmw + sp*( 54.6746e+00_RLEN - 0.603459e+00_RLEN*X + 1.09987e-02_RLEN*X**2   &
         &  - 6.1670e-05_RLEN*X**3)                                                      &
         &  + sp*SQRT(sp)*( 7.944e-02_RLEN + 1.6483e-02_RLEN*X - 5.3009e-04_RLEN*X**2)
   
   ! Secant bulk modulus of seawater at S,T,P
   Ksbm = Ksbm0 &
        &  + pr*(3.239908e+00_RLEN+ 1.43713e-03_RLEN*X + 1.16092e-04_RLEN*X**2 - 5.77905e-07_RLEN*X**3) &
        &  + pr*sp*(2.2838e-03_RLEN - 1.0981e-05_RLEN*X - 1.60780e-06_RLEN*X**2)                        &
        &  + pr*sp*SQRT(sp)*1.91075e-04_RLEN                                                            &
        &  + pr*pr*(8.50935e-05_RLEN - 6.12293e-06_RLEN*X + 5.2787e-08_RLEN*X**2)                       &
        &  + pr*pr*sp*(-9.9348e-07_RLEN + 2.0816e-08_RLEN*X + 9.1697e-10_RLEN*X**2)
   
   ! Density of seawater at S,T,P
     sw_rho_t = rho0/(1.0e+00_RLEN - pr/Ksbm)

   return
   end function

   !--------------------------------------------------------------------------


   ! TEOS-10 Functions gsw V3.05

   elemental function gsw_p_from_z (z, lat)
   !==========================================================================
   ! Calculates the pressure p from the height z
   ! z      : height                                   [m]
   ! lat    : latitude                                 [deg]
   ! gsw_p_from_z : pressure                           [dbar]
   !
   ! NOTE: If the graviational acceleration were to be regarded as being depth-independent,
   ! which is often the case in ocean models, then gamma would be set to be zero here,
   ! and the code below works perfectly well
   ! geo_strf_dyn_height and sea_surface_geopotental are here neglected
   !==========================================================================
   implicit none

   real (RLEN), intent(in) :: z, lat
   real (RLEN) :: gsw_p_from_z
   real (RLEN) :: sin2, gs, c1, p, df_dp, f, p_mid
   real (RLEN) :: gamma

   gamma = 0._RLEN
   sin2 = ( sin(lat*deg2rad) ) ** 2
   gs = 9.780327_RLEN*(1.0_RLEN + (5.2792e-3_RLEN + (2.32e-5_RLEN*sin2)) * sin2)

   ! get the first estimate of p from Saunders (1981)
   c1 =  5.25e-3_RLEN * sin2 + 5.92e-3_RLEN
   p  = -2._RLEN*z/((1._RLEN-c1) + sqrt((1._RLEN-c1)*(1._RLEN-c1) + 8.84e-6_RLEN*z))

   ! initial value of the derivative of f
   df_dp = db2Pa * gsw_specvol_sso_0(p)

   f = gsw_enthalpy_sso_0(p) + gs * (z - 0.5_RLEN * gamma * (z**2))

   p_mid = 0.5 * ( p + (p - (f/df_dp)) );
   df_dp = db2Pa * gsw_specvol_sso_0(p_mid);

   gsw_p_from_z = p - f/df_dp;

   return
   end function

   !--------------------------------------------------------------------------

   elemental function gsw_specvol_sso_0 (p)
   !==========================================================================
   !  This function calculates specifc volume at the Standard Ocean Salinity,
   !  SSO, and at a Conservative Temperature of zero degrees C, as a function
   !  of pressure, p, in dbar, using a streamlined version of the CT version
   !  of specific volume, that is, a streamlined version of the code
   !  "gsw_specvol(SA,CT,p)".
   !==========================================================================

   implicit none

   real (RLEN), intent(in) :: p
   real (RLEN) :: gsw_specvol_sso_0
   real (RLEN) :: z
   ! gsw specvol coefficients
   real (RLEN), parameter :: v005 = -1.2647261286e-8_RLEN
   real (RLEN), parameter :: v006 =  1.9613503930e-9_RLEN

   z = p*1e-4_RLEN

   gsw_specvol_sso_0 = &
            9.726613854843870e-04_RLEN + z*(-4.505913211160929e-05_RLEN &
       + z*(7.130728965927127e-06_RLEN + z*(-6.657179479768312e-07_RLEN &
       + z*(-2.994054447232880e-08_RLEN + z*(v005 + v006*z)))))

   return
   end function

   !--------------------------------------------------------------------------

   elemental function gsw_enthalpy_sso_0 (p)
   !==========================================================================
   !  This function calculates enthalpy at the Standard Ocean Salinity, SSO,
   !  and at a Conservative Temperature of zero degrees C, as a function of
   !  pressure, p, in dbar, using a streamlined version of the
   !  computationally-efficient expression for specific volume, that is, a
   !  streamlined version of the code "gsw_enthalpy(SA,CT,p)".
   !==========================================================================
   implicit none

   real (RLEN), intent(in) :: p
   real (RLEN) :: gsw_enthalpy_sso_0
   real (RLEN) :: dynamic_enthalpy_sso_0_p, z
   ! gsw specvol coefficients
   real (RLEN), parameter :: h006 = -2.10787688100e-9_RLEN
   real (RLEN), parameter :: h007 =  2.80192913290e-10_RLEN

   z = p * 1e-4_RLEN

   dynamic_enthalpy_sso_0_p = &
         z*( 9.726613854843870e-4_RLEN + z*(-2.252956605630465e-5_RLEN &
       + z*( 2.376909655387404e-6_RLEN + z*(-1.664294869986011e-7_RLEN &
       + z*(-5.988108894465758e-9_RLEN + z*(h006 + h007*z))))))

   gsw_enthalpy_sso_0 = dynamic_enthalpy_sso_0_p * db2pa * 1e4_RLEN

   return
   end function

   !--------------------------------------------------------------------------

   end module sw_tool

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
