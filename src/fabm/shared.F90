module ogs_bfm_shared
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

   use fabm_types
   use fabm_standard_variables

   implicit none

   public

   integer,parameter :: LightPeriodFlag=1

   real(rk),parameter :: CMass         = 12.011_rk
   real(rk),parameter :: ONE           = 1._rk
   real(rk),parameter :: ZERO          = 0._rk
   real(rk),parameter :: BASETEMP      = 10._rk
   real(rk),parameter :: p_small       = 1.0E-20_rk
   real(rk),parameter :: qnRPIcX       = 1.26E-02_rk
   real(rk),parameter :: qpRPIcX       = 7.86E-04_rk
   real(rk),parameter :: qsRPIcX       = 15._rk/106._rk/CMass
   real(rk),parameter :: ZeroX         = 1e-8_rk
   real(rk),parameter :: pi            = acos(-1._rk)
   real(rk),parameter :: deg2rad       = pi/180._rk
   real(rk),parameter :: WtoQuanta     = 4.57_rk
   real(rk),parameter :: SEC_PER_DAY   = 86400.0_rk
   real(rk),parameter :: SUNQ          = 24.0_rk
   real(rk),parameter :: HOURS_PER_DAY = 24.0_rk
   real(rk),parameter :: MW_C          = 12.0_rk
   real(rk),parameter :: C2ALK         = 2.0_rk/MW_C   ! Conversion factor between inorganic carbon and alkalinity
   real(rk),parameter :: p_atm0         = 1013.25_rk    !reference sea level pressure
   real(rk),parameter :: ZERO_KELVIN   = -273.15_rk;
   real(rk),parameter :: h_planck      = 6.6256E-34   !Plancks constant J sec
   real(rk),parameter :: c_light       = 2.998E8      !speed of light m/sec
   real(rk),parameter :: oavo          = 1.0D0/6.023E23   ! 1/Avogadros number

   real(rk)           :: flPTN6r    ! total rate of formation of reduction equivalent [mmolHS/m3/d] computed in PelBac and used in PelChem   
   real(rk)           :: qccPPY     ! PIC:POC ration in P2: compputed in Phyto and used in MicroZoo and MesoZoo only for prey P2 (crapy solution)

#ifdef IRON
   logical,parameter :: use_iron = .true.
#else
   logical,parameter :: use_iron = .false.
#endif

   ! Aggregate diagnostics for e.g., carbon budgets.
   type (type_bulk_standard_variable),parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_oxygen = type_bulk_standard_variable(name='total_oxygen',units='mmolO2/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_reduction_equivalent = type_bulk_standard_variable(name='total_reduction_equivalent',units='mmolEq/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: alkalinity = type_bulk_standard_variable(name='alkalinity',units='mmolEq/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_calcite_in_biota = type_bulk_standard_variable(name='total_calcite_in_biota',units='mg C/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_silicate = type_bulk_standard_variable(name='total_silicate',units='mmolSi/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: secchi_depth = type_bulk_standard_variable(name='secchi_depth',units='m')

   ! Aggregate variables for benthic bioturbation and bioirrigation (summed over all fauna).
   type (type_bulk_standard_variable),parameter :: total_bioturbation_activity = type_bulk_standard_variable(name='total_bioturbation_activity',units='mg C/m^2/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_bioirrigation_activity = type_bulk_standard_variable(name='total_bioirrigation_activity',units='mg C/m^2/d',aggregate_variable=.true.)

   ! Spectral light variables
   real(rk), allocatable, dimension(:)                 :: lam,lam1,lam2,aw,bw,bbw,apoc,bpoc,bbpoc,WtoQ,acdom_min
   real(rk), allocatable, dimension(:)                 :: Ed_0,Es_0
   real(rk), allocatable, dimension(:,:)               :: ac,ac_ps,bc,bbc,acdom
!  real(rk), allocatable, dimension(:,:)               :: a_array, b_array, bb_array
!  type (type_surface_standard_variable),parameter     :: surf_direct_downward_irradiance_250_nm = type_surface_standard_variable(name='surf_direct_downward_irradiance_250_nm',units='W/m2')
!  type (type_surface_standard_variable),parameter     :: surf_diffuse_downward_irradiance_250_nm = type_surface_standard_variable(name='surf_diffuse_downward_irradiance_250_nm',units='W/m2')
!  type (type_surface_standard_variable),parameter     :: surf_diffuse_direct_irradiance_325_nm = type_surface_standard_variable(name='surf_direct_downward_irradiance_325_nm',units='W/m2')
!  type (type_surface_standard_variable),parameter     :: surf_diffuse_downward_irradiance_325_nm = type_surface_standard_variable(name='surf_diffuse_downward_irradiance_325_nm',units='W/m2')
    type (type_bulk_standard_variable),parameter       :: PAR_tot = type_bulk_standard_variable(name='PAR_tot',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_dia = type_bulk_standard_variable(name='PAR_dia',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_flag = type_bulk_standard_variable(name='PAR_flag',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_pico = type_bulk_standard_variable(name='PAR_pico',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_dino = type_bulk_standard_variable(name='PAR_dino',units='<UNITS>',aggregate_variable=.true.)

   ! Standard benthic variables used to make implicit based on matching standard names coupling possible.
   type (type_horizontal_standard_variable),parameter :: depth_of_sediment_column = type_horizontal_standard_variable(name='depth_of_sediment_column',units='m')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_1 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_1',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_2 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_2',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_3 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_3',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: particulate_diffusivity_due_to_bioturbation = type_horizontal_standard_variable(name='particulate_diffusivity_due_to_bioturbation',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: bioturbation_depth = type_horizontal_standard_variable(name='bioturbation_depth',units='m')
   type (type_horizontal_standard_variable),parameter :: sediment_porosity = type_horizontal_standard_variable(name='sediment_porosity',units='-')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_1 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_1',units='m')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_2 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_2',units='m')
   type (type_horizontal_standard_variable),parameter :: pelagic_benthic_transfer_constant = type_horizontal_standard_variable(name='pelagic_benthic_transfer_constant',units='d/m')
   type (type_horizontal_standard_variable),parameter :: sediment_erosion = type_horizontal_standard_variable(name='sediment_erosion',units='m/d')

   ! Aggregate absorption and backscatter.
   type (type_bulk_standard_variable),parameter :: particulate_organic_absorption_coefficient = type_bulk_standard_variable(name='particulate_organic_absorption_coefficient',units='1/m',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: particulate_organic_backscatter_coefficient = type_bulk_standard_variable(name='particulate_organic_backscatter_coefficient',units='1/m',aggregate_variable=.true.)

   ! Gelbstoff absorption.
   type (type_horizontal_standard_variable),parameter :: gelbstoff_absorption_from_satellite = type_horizontal_standard_variable(name='gelbstoff_absorption_from_satellite',units='1/m')

   ! Zenith angle.
   type (type_horizontal_standard_variable),parameter :: zenith_angle = type_horizontal_standard_variable(name='zenith_angle',units='degrees')

   contains

    ! temperature dependency for Q10 function
    elemental function eTq(t, q10)

        IMPLICIT NONE
        real(rk),intent(IN) :: t, q10
        real(rk)            :: eTq

        eTq = exp( log(q10) * (t-BASETEMP) / BASETEMP)

    end function eTq

    ! Michaelis-Menten saturation curve
    elemental FUNCTION MM(x, m)

        IMPLICIT NONE
        real(rk),intent(IN) :: x, m
        real(rk)            :: MM

        MM = x / (x + m)

    end function MM

    ! Convert values in 0 or 1 according to input field
    elemental FUNCTION INSW(x)

        IMPLICIT NONE
        real(rk),intent(IN) :: x
        real(rk)            :: INSW
 
        INSW = ZERO
        if (x > ZERO ) INSW=ONE 

    end function INSW

    ! Michaelis-Menten saturation curve at power
    elemental FUNCTION MM_POWER(x, m, p)

        IMPLICIT NONE
        real(rk),intent(IN) :: x, m
        integer   ,intent(IN) :: p
        real(rk)            :: MM_POWER

        MM_POWER = x**p / ( x**p+ m**p)

    end function MM_POWER

    subroutine linear_regression(x, y, n, a, b)
      
        IMPLICIT NONE
        real(rk), intent(IN) :: x(:), y(:)
        integer,  intent(IN) :: n
        real(rk), intent(OUT) :: a, b
        real(rk) :: s1,s2,s3,s4
        integer :: i
        do i=1,n
           s1=s1+x(i)
           s2=s2+x(i)**2
           s3=s3+y(i)
           s4=s4+x(i)*y(i)
        enddo
        b=((n*s4)-(s1*s3))/((n*s2)-(s1**2))
        a=(s3-(s1*b))/n
        
    end subroutine linear_regression


end module
