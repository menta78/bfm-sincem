#include "fabm_driver.h"

module ogs_bfm_light_spectral
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
   use ogs_bfm_shared
   use adj_3stream, only: solve_direct

   implicit none

   private

   type,extends(type_base_model),public :: type_ogs_bfm_light_spectral
      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_par_dia, id_par_flag, id_par_pico, id_par_dino
      type (type_diagnostic_variable_id)   :: id_PAR_tot
      type (type_diagnostic_variable_id)   :: id_anap450, id_aph450
      type (type_diagnostic_variable_id)   :: id_acdom250, id_acdom325, id_acdom400, id_acdom425, id_acdom450
      type (type_diagnostic_variable_id)   :: id_Scdom350_500, id_Scdom250_325       
      type (type_diagnostic_variable_id)   :: id_bbp450, id_bbp550, id_bbp700
      type (type_horizontal_diagnostic_variable_id) :: id_Rrs400, id_Rrs425, id_Rrs450, id_Rrs475
      type (type_horizontal_diagnostic_variable_id) :: id_Rrs500, id_Rrs525, id_Rrs550, id_Rrs575, id_Rrs675
      type (type_horizontal_diagnostic_variable_id) :: id_kd375, id_kd400, id_kd425, id_kd475, id_kd500
      
      type (type_dependency_id)            :: id_dz
      type (type_state_variable_id)        :: id_P1c, id_P2c, id_P3c, id_P4c
      type (type_state_variable_id)        :: id_P1chl, id_P2chl, id_P3chl, id_P4chl
      type (type_state_variable_id)        :: id_R6c, id_X1c, id_X2c, id_X3c
      type (type_horizontal_dependency_id) :: id_zenithA

! BLOCK 1 python generated code see AUX_SCRIPTS/python_light_spectral.py
      type (type_horizontal_dependency_id) ::  id_Ed_0_0250, id_Ed_0_0325, id_Ed_0_0350, id_Ed_0_0375, id_Ed_0_0400
      type (type_horizontal_dependency_id) ::  id_Ed_0_0425, id_Ed_0_0450, id_Ed_0_0475, id_Ed_0_0500, id_Ed_0_0525
      type (type_horizontal_dependency_id) ::  id_Ed_0_0550, id_Ed_0_0575, id_Ed_0_0600, id_Ed_0_0625, id_Ed_0_0650
      type (type_horizontal_dependency_id) ::  id_Ed_0_0675, id_Ed_0_0700, id_Ed_0_0725, id_Ed_0_0775, id_Ed_0_0850
      type (type_horizontal_dependency_id) ::  id_Ed_0_0950, id_Ed_0_1050, id_Ed_0_1150, id_Ed_0_1250, id_Ed_0_1350
      type (type_horizontal_dependency_id) ::  id_Ed_0_1450, id_Ed_0_1550, id_Ed_0_1650, id_Ed_0_1750, id_Ed_0_1900
      type (type_horizontal_dependency_id) ::  id_Ed_0_2200, id_Ed_0_2900, id_Ed_0_3700
      type (type_horizontal_dependency_id) ::  id_Es_0_0250, id_Es_0_0325, id_Es_0_0350, id_Es_0_0375, id_Es_0_0400
      type (type_horizontal_dependency_id) ::  id_Es_0_0425, id_Es_0_0450, id_Es_0_0475, id_Es_0_0500, id_Es_0_0525
      type (type_horizontal_dependency_id) ::  id_Es_0_0550, id_Es_0_0575, id_Es_0_0600, id_Es_0_0625, id_Es_0_0650
      type (type_horizontal_dependency_id) ::  id_Es_0_0675, id_Es_0_0700, id_Es_0_0725, id_Es_0_0775, id_Es_0_0850
      type (type_horizontal_dependency_id) ::  id_Es_0_0950, id_Es_0_1050, id_Es_0_1150, id_Es_0_1250, id_Es_0_1350
      type (type_horizontal_dependency_id) ::  id_Es_0_1450, id_Es_0_1550, id_Es_0_1650, id_Es_0_1750, id_Es_0_1900
      type (type_horizontal_dependency_id) ::  id_Es_0_2200, id_Es_0_2900, id_Es_0_3700

! END BLOCK 1 python generated  code

      ! Parameters
      integer  :: nlt,npft
      real(rk) :: rd, rs, ru, vs, vu
      real(rk) :: SdomX1, X1coeff, SdomX2, X2coeff, SdomX3, X3coeff, lambda_aCDOM, Xmincoeff
      real(rk) :: Sapar, lambda_aPart, aparcoeff
      real(rk) :: Sbpar, lambda_bPart, bparcoeff, bb_to_b
      logical :: compute_acdom
      logical :: compute_anap
      real(rk) :: p_epsP1, p_epsP2, p_epsP3, p_epsP4
      real(rk) :: p_bpsP1, p_bpsP2, p_bpsP3, p_bpsP4
      real(rk) :: p_bbrP1, p_bbrP2, p_bbrP3, p_bbrP4
      logical :: compute_aph, compute_bph, compute_bbc
      
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do_column
   end type type_ogs_bfm_light_spectral

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ogs_bfm_light_spectral),intent(inout),target :: self
      integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      real(rk) :: hc, hcoavo, rlamm, rlamm1, rlamm2, nl
      real(rk) :: n, dummy_p, cu_area, aph_mean, bph_mean

!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%nlt,    'nlt',  '-',   'number of wavelenghts', default=-1)
      call self%get_parameter(self%npft,   'npft', '-',   'number of PFT', default=-1)
      call self%get_parameter(self%rd,     'rd',   '-',   ' ', default=1.0_rk)
      call self%get_parameter(self%rs,     'rs',   '-',   ' ', default=1.5_rk)
      call self%get_parameter(self%ru,     'ru',   '-',   ' ', default=3.0_rk)
      call self%get_parameter(self%vs,     'vs',   '',    'avg cosine diffuse down', default=0.83_rk)
      call self%get_parameter(self%vu,     'vu',   '-',   'avg cosine diffuse up', default=0.4_rk)
      call self%get_parameter(self%SdomX1,        'SdomX1',        'nm-1',     'slope for aCDOM [X1c] wavelength dependence')
      call self%get_parameter(self%X1coeff,       'X1coeff',       'm2 mgC-1', 'specific absorption of X1c at lambda_aCDOM ')
      call self%get_parameter(self%SdomX2,        'SdomX2',        'nm-1',     'slope for aCDOM [X2c] wavelength dependence')
      call self%get_parameter(self%X2coeff,       'X2coeff',       'm2 mgC-1', 'specific absorption of X2c at lambda_aCDOM ')
      call self%get_parameter(self%SdomX3,        'SdomX3',        'nm-1',     'slope for aCDOM [X3c] wavelength dependence')
      call self%get_parameter(self%X3coeff,       'X3coeff',       'm2 mgC-1', 'specific absorption of X3c at lambda_aCDOM ')
      call self%get_parameter(self%lambda_aCDOM,  'lambda_aCDOM',  'nm',       'wavelength where reference aCDOM is given')
      call self%get_parameter(self%Xmincoeff,     'Xmincoeff',     'm-1',      'minimum aCDOM at 450nm')      
      call self%get_parameter(self%Sapar,         'Sapar',         'nm-1',     'slope parameter for aNAP wavelength dependence')
      call self%get_parameter(self%lambda_aPart,  'lambda_aPart',  'nm',       'wavelength where reference aNAP is given')
      call self%get_parameter(self%aparcoeff,     'aparcoeff',     'm2 mgC-1', 'specific absorption at lambda_aPart ')
      call self%get_parameter(self%Sbpar,         'Sbpar',         '-',        'exponent for bNAP wavelength dependence')
      call self%get_parameter(self%lambda_bPart,  'lambda_bPart',  'nm',       'wavelength where reference bNAP is given')
      call self%get_parameter(self%bparcoeff,     'bparcoeff',     'm2 mgC-1', 'specific scatter at lambda_bPart')
      call self%get_parameter(self%bb_to_b,       'bb_to_b',       '-',        'backscatter to total scatter ratio') 
      call self%get_parameter(self%compute_acdom, 'compute_acdom', '[T or F]', 'logical flag to compute acdom') 
      call self%get_parameter(self%compute_anap,  'compute_anap',  '[T or F]', 'logical flag to compute anap')
      call self%get_parameter(self%p_epsP1,       'p_epsP1',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P1', default=0.03_rk)
      call self%get_parameter(self%p_epsP2,       'p_epsP2',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P2', default=0.03_rk)
      call self%get_parameter(self%p_epsP3,       'p_epsP3',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P3', default=0.03_rk)
      call self%get_parameter(self%p_epsP4,       'p_epsP4',       'm2 mgChl-1',  'mean absorption coefficient from 400-700nm for P4', default=0.03_rk)
      call self%get_parameter(self%compute_aph,   'compute_aph',   '[T or F]',    'logical flag to scale aph')
      call self%get_parameter(self%p_bpsP1,       'p_bpsP1',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P1')
      call self%get_parameter(self%p_bpsP2,       'p_bpsP2',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P2')
      call self%get_parameter(self%p_bpsP3,       'p_bpsP3',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P3')
      call self%get_parameter(self%p_bpsP4,       'p_bpsP4',       'm2 mgC-1',  'mean scattering coefficient from 400-700nm for P4')
      call self%get_parameter(self%compute_bph,   'compute_bph',   '[T or F]',    'logical flag to scale bph')
      call self%get_parameter(self%p_bbrP1,       'p_bbrP1',       '-',           'backscattering to total scattering ratio for P1')
      call self%get_parameter(self%p_bbrP2,       'p_bbrP2',       '-',           'backscattering to total scattering ratio for P2')
      call self%get_parameter(self%p_bbrP3,       'p_bbrP3',       '-',           'backscattering to total scattering ratio for P3')
      call self%get_parameter(self%p_bbrP4,       'p_bbrP4',       '-',           'backscattering to total scattering ratio for P4')
      call self%get_parameter(self%compute_bbc,   'compute_bbc',   '[T or F]',    'logical flag to compute bbc from bc*bbr')
      
      if (self%nlt>0) then
          allocate(lam(self%nlt));             lam(:)=huge(lam(1))
          allocate(lam1(self%nlt));            lam1(:)=huge(lam1(1))
          allocate(lam2(self%nlt));            lam2(:)=huge(lam2(1))
          allocate(aw(self%nlt));              aw(:)=huge(aw(1))
          allocate(bw(self%nlt));              bw(:)=huge(bw(1))
          allocate(bbw(self%nlt));             bbw(:)=huge(bbw(1))
          allocate(ac(self%npft,self%nlt));    ac(:,:)=huge(ac(1,1))
          allocate(ac_ps(self%npft,self%nlt)); ac_ps(:,:)=huge(ac_ps(1,1))
          allocate(bc(self%npft,self%nlt));    bc(:,:)=huge(bc(1,1))
          allocate(bbc(self%npft,self%nlt));   bbc(:,:)=huge(bbc(1,1))
          allocate(apoc(self%nlt));            apoc(:)=huge(apoc(1))
          allocate(bpoc(self%nlt));            bpoc(:)=huge(bpoc(1))
          allocate(bbpoc(self%nlt));           bbpoc(:)=huge(bbpoc(1))
          allocate(acdom_min(self%nlt));       acdom_min(:)=huge(acdom_min(1))
          allocate(acdom(3,self%nlt));         acdom(:,:)=huge(acdom(1,1))
          allocate(Ed_0(self%nlt));            Ed_0(:)=huge(Ed_0(1))
          allocate(Es_0(self%nlt));            Es_0(:)=huge(Es_0(1))
          allocate(WtoQ(self%nlt));            WtoQ(:)=huge(WtoQ(1))
!          allocate(equis(7));                  equis(:)=huge(equis(1))
!          allocate(ies(7));                    ies(:)=huge(ies(1))

          
!         load the IOP for the biogeochemical variables considered
          call lidata(self%nlt,self%npft)

      endif 

     hc = 1.0D0/(h_planck*c_light)
     hcoavo = hc*oavo


      do nl = 1,self%nlt
       rlamm = real(lam(nl),8)*1.0E-9      !lambda in m
       WtoQ(nl) = rlamm*hcoavo*1000000.0D0 !Watts to micro mol quanta conversion
       acdom_min(nl)= 0.0_rk
      enddo

!   CDOM minimum absorption (m-1)
!   hardcoded
!      acdom_min(1)= 0.1715_rk
!      acdom_min(2)= 0.0373_rk
!      acdom_min(3)= 0.0239_rk
!      acdom_min(4)= 0.0153_rk
!      acdom_min(5)= 0.0098_rk
!      acdom_min(6)= 0.0063_rk
!      acdom_min(7)= 0.0040_rk
!      acdom_min(8)= 0.0026_rk
!      acdom_min(9)= 0.0016_rk
!      acdom_min(10)= 0.0010_rk
!      acdom_min(11)= 0.00068_rk
!      acdom_min(12)= 0.00043_rk
!      acdom_min(13)= 0.00028_rk
!      acdom_min(14)= 0.00018_rk
!      acdom_min(15)= 0.00011_rk
!      acdom_min(16)= 0.000073_rk
!      acdom_min(17)= 0.000047_rk
!      acdom_min(18)= 0.000027_rk
!!      acdom_min(19)= 0.000013_rk
!!      acdom_min(20)= 0.0000037_rk
!!      acdom_min(21)= 0.00000062_rk
!!      acdom_min(22)= 0.00000010_rk
!!      acdom_min(23)= 0.000000017_rk
!!      acdom_min(24)= 0.0000000030_rk
!!      acdom_min(25)= 0.00000000050_rk
!!      acdom_min(26)= 0.000000000084_rk
!!      acdom_min(27)= 0.000000000014_rk
!!      acdom_min(28)= 0.0000000000024_rk
!!      acdom_min(29)= 0.00000000000040_rk      

!   with parameter in namelist     
      do nl = 1,self%nlt
       rlamm1 = real(lam1(nl),8)
       rlamm2 = real(lam2(nl),8)
       acdom_min(nl) = self%Xmincoeff*(exp(-self%SdomX2*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX2*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX2*(rlamm2-rlamm1))
      enddo

!      do nl = 1,self%nlt
!        write(*,*) real(lam(nl),8), acdom_min(nl)
!      enddo
      
      
 !   CDOM absorption coefficients
      if (self%compute_acdom) then     
      do nl = 1,self%nlt
 !      rlamm = real(lam(nl),8)
       rlamm1 = real(lam1(nl),8)
       rlamm2 = real(lam2(nl),8)
 !      acdom(nl) = self%cdomcoeff * exp(-self%Sdom*(rlamm-self%lambda_aCDOM))
       acdom(1,nl) = self%X1coeff*(exp(-self%SdomX1*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX1*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX1*(rlamm2-rlamm1))
       acdom(2,nl) = self%X2coeff*(exp(-self%SdomX2*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX2*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX2*(rlamm2-rlamm1))
       acdom(3,nl) = self%X3coeff*(exp(-self%SdomX3*(rlamm2-self%lambda_aCDOM))-exp(-self%SdomX3*(rlamm1-self%lambda_aCDOM)))/(-self%SdomX3*(rlamm2-rlamm1))
      enddo
      endif

 !   POC absorption/scatter coefficients
      if (self%compute_anap) then
      do nl = 1,self%nlt
       rlamm = real(lam(nl),8)
       rlamm1 = real(lam1(nl),8)
       rlamm2 = real(lam2(nl),8)
 !      apoc(nl) = self%aparcoeff * exp(-self%Sapar*(rlamm-self%lambda_aPart))
       apoc(nl) = self%aparcoeff*(exp(-self%Sapar*(rlamm2-self%lambda_aPart))-exp(-self%Sapar*(rlamm1-self%lambda_aPart)))/(-self%Sapar*(rlamm2-rlamm1))     
       bpoc(nl) = self%bparcoeff * ((self%lambda_bPart/rlamm)**self%Sbpar)
       bbpoc(nl) = self%bparcoeff * ((self%lambda_bPart/rlamm)**self%Sbpar) * self%bb_to_b
      enddo
      endif

!       write(*,*) "acdom", acdom(1,7), acdom(2,7), acdom(3,7)
!       write(*,*) "apoc", apoc(7)
!       write(*,*) "bpoc", bpoc(7)
!       write(*,*) "bbpoc", bbpoc(7)

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), acdom(1,nl), acdom(2,nl), acdom(3,nl)
      enddo

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), apoc(nl), bpoc(nl), bbpoc(nl)
      enddo
       
 !   PHYTO absorption coefficients
      if (self%compute_aph) then
         do n = 1,self%npft
            if (n == 1) dummy_p = self%p_epsP1
            if (n == 2) dummy_p = self%p_epsP2
            if (n == 3) dummy_p = self%p_epsP3
            if (n == 4) dummy_p = self%p_epsP4
            ! find mean
            cu_area = 0.0_rk
            do nl = 5,17   ! indexes for 400-700nm
!                rlamm = real(lam(nl),8)
                rlamm1 = real(lam1(nl),8)
                rlamm2 = real(lam2(nl),8)
                cu_area = cu_area + ((rlamm2-rlamm1) * ac(n,nl))
            enddo
            aph_mean = cu_area / (real(lam2(17),8)-real(lam1(5),8)) ! total band width 400-700nm
            do nl = 1,self%nlt
                ac(n,nl) = ac(n,nl) * (dummy_p/aph_mean)
            enddo 
         enddo
      endif
       
 !   PHYTO scattering coefficients
      if (self%compute_bph) then
         do n = 1,self%npft
            if (n == 1) dummy_p = self%p_bpsP1
            if (n == 2) dummy_p = self%p_bpsP2
            if (n == 3) dummy_p = self%p_bpsP3
            if (n == 4) dummy_p = self%p_bpsP4
            ! find mean
            cu_area = 0.0_rk
            do nl = 5,17   ! indexes for 400-700nm
                rlamm1 = real(lam1(nl),8)
                rlamm2 = real(lam2(nl),8)
                cu_area = cu_area + ((rlamm2-rlamm1) * bc(n,nl))
            enddo
            bph_mean = cu_area / (real(lam2(17),8)-real(lam1(5),8)) ! total band width 400-700nm
            do nl = 1,self%nlt
                bc(n,nl) = bc(n,nl) * (dummy_p/bph_mean)
            enddo 
         enddo
      endif

 !   PHYTO backscattering coefficients       
      if (self%compute_bbc) then
       do nl = 1,19
          bbc(1,nl) = self%p_bbrP1
          bbc(2,nl) = self%p_bbrP2
          bbc(3,nl) = self%p_bbrP3
          bbc(4,nl) = self%p_bbrP4
       enddo
      endif

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), ac(1,nl), ac_ps(1,nl), bc(1,nl), bbc(1,nl)
      enddo

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), ac(2,nl), ac_ps(2,nl), bc(2,nl), bbc(2,nl)
      enddo

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), ac(3,nl), ac_ps(3,nl), bc(3,nl), bbc(3,nl)
      enddo

      do nl = 1,self%nlt
        write(*,*) real(lam(nl),8), ac(4,nl), ac_ps(4,nl), bc(4,nl), bbc(4,nl)
      enddo
      
     ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_par_dia, 'PAR_dia',  'uE mgChl-1 d-1', 'PAR_diatoms', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_flag,'PAR_flag', 'uE mgChl-1 d-1', 'PAR_flagellates', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_pico,'PAR_pico', 'uE mgChl-1 d-1', 'PAR_picophytoplankton', source=source_do_column)
      call self%register_diagnostic_variable(self%id_par_dino,'PAR_dino', 'uE mgChl-1 d-1', 'PAR_dinoflagellates', source=source_do_column)
      call self%register_diagnostic_variable(self%id_PAR_tot, 'PAR_tot',  'uE m-2 d-1 [400-700]','PAR_total', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Scdom350_500, 'Scdom350_500', 'nm-1','visible spectral slope acdom', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Scdom250_325, 'Scdom250_325', 'nm-1','UV spectral slope acdom', source=source_do_column)
      call self%register_diagnostic_variable(self%id_acdom250, 'acdom250', 'm-1', 'acdom in 250 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_acdom325, 'acdom325', 'm-1', 'acdom in 325 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_acdom400, 'acdom400', 'm-1', 'acdom in 400 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_acdom425, 'acdom425', 'm-1', 'acdom in 425 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_acdom450, 'acdom450', 'm-1', 'acdom in 450 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_anap450, 'anap450',  'm-1', 'anap in 450 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_aph450,  'aph450',   'm-1', 'aph in 450 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_bbp450,  'bbp450',   'm-1', 'particle backscattering in 450 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_bbp550,  'bbp550',   'm-1', 'particle backscattering in 550 nm band', source=source_do_column)      
      call self%register_diagnostic_variable(self%id_bbp700,  'bbp700',   'm-1', 'particle backscattering in 700 nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs400,   'Rrs400',   '-',  'subsurface reflectance in 400nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs425,   'Rrs425',   '-',  'subsurface reflectance in 425nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs450,   'Rrs450',   '-',  'subsurface reflectance in 450nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs475,   'Rrs475',   '-',  'subsurface reflectance in 475nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs500,   'Rrs500',   '-',  'subsurface reflectance in 500nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs525,   'Rrs525',   '-',  'subsurface reflectance in 525nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs550,   'Rrs550',   '-',  'subsurface reflectance in 550nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs575,   'Rrs575',   '-',  'subsurface reflectance in 575nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Rrs675,   'Rrs675',   '-',  'subsurface reflectance in 675nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_kd375,    'kd375',  'm-1',  'extinction coefficient in 375nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_kd400,    'kd400',  'm-1',  'extinction coefficient in 400nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_kd425,    'kd425',  'm-1',  'extinction coefficient in 425nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_kd475,    'kd475',  'm-1',  'extinction coefficient in 475nm band', source=source_do_column)
      call self%register_diagnostic_variable(self%id_kd500,    'kd500',  'm-1',  'extinction coefficient in 500nm band', source=source_do_column)
      
      ! Register biogeochemical dependencies
      call self%register_state_dependency(self%id_P1c,'P1c','mg C/m^3', 'Diatoms carbon')
      call self%register_state_dependency(self%id_P2c,'P2c','mg C/m^3', 'Flagellates carbon')
      call self%register_state_dependency(self%id_P3c,'P3c','mg C/m^3', 'PicoPhytoplankton carbon')
      call self%register_state_dependency(self%id_P4c,'P4c','mg C/m^3', 'DinoFlagellates carbon')

      call self%register_state_dependency(self%id_P1chl,'P1chl','mg chl/m^3', 'Diatoms chlorophyll')
      call self%register_state_dependency(self%id_P2chl,'P2chl','mg chl/m^3', 'Flagellates chlorophyll')
      call self%register_state_dependency(self%id_P3chl,'P3chl','mg chl/m^3', 'PicoPhytoplankton chlorophyll')
      call self%register_state_dependency(self%id_P4chl,'P4chl','mg chl/m^3', 'DinoFlagellates chlorophyll')

      call self%register_state_dependency(self%id_R6c,'R6c','mg C/m^3', 'POC')
      call self%register_state_dependency(self%id_X1c,'X1c','mg C/m^3', 'labile CDOM')
      call self%register_state_dependency(self%id_X2c,'X2c','mg C/m^3', 'semi-labile CDOM')
      call self%register_state_dependency(self%id_X3c,'X3c','mg C/m^3', 'semi-refractory CDOM')

      ! Register environmental dependencies 
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_horizontal_dependency(self%id_zenithA, type_horizontal_standard_variable(name='zenith_angle'))

! BLOCK 2 python generate code see AUX_SCRIPTS/python_light_spectral.py
      call self%register_dependency(self%id_Ed_0_0250,type_surface_standard_variable(name='surf_direct_downward_irradiance_0250_nm'))
      call self%register_dependency(self%id_Ed_0_0325,type_surface_standard_variable(name='surf_direct_downward_irradiance_0325_nm'))
      call self%register_dependency(self%id_Ed_0_0350,type_surface_standard_variable(name='surf_direct_downward_irradiance_0350_nm'))
      call self%register_dependency(self%id_Ed_0_0375,type_surface_standard_variable(name='surf_direct_downward_irradiance_0375_nm'))
      call self%register_dependency(self%id_Ed_0_0400,type_surface_standard_variable(name='surf_direct_downward_irradiance_0400_nm'))
      call self%register_dependency(self%id_Ed_0_0425,type_surface_standard_variable(name='surf_direct_downward_irradiance_0425_nm'))
      call self%register_dependency(self%id_Ed_0_0450,type_surface_standard_variable(name='surf_direct_downward_irradiance_0450_nm'))
      call self%register_dependency(self%id_Ed_0_0475,type_surface_standard_variable(name='surf_direct_downward_irradiance_0475_nm'))
      call self%register_dependency(self%id_Ed_0_0500,type_surface_standard_variable(name='surf_direct_downward_irradiance_0500_nm'))
      call self%register_dependency(self%id_Ed_0_0525,type_surface_standard_variable(name='surf_direct_downward_irradiance_0525_nm'))
      call self%register_dependency(self%id_Ed_0_0550,type_surface_standard_variable(name='surf_direct_downward_irradiance_0550_nm'))
      call self%register_dependency(self%id_Ed_0_0575,type_surface_standard_variable(name='surf_direct_downward_irradiance_0575_nm'))
      call self%register_dependency(self%id_Ed_0_0600,type_surface_standard_variable(name='surf_direct_downward_irradiance_0600_nm'))
      call self%register_dependency(self%id_Ed_0_0625,type_surface_standard_variable(name='surf_direct_downward_irradiance_0625_nm'))
      call self%register_dependency(self%id_Ed_0_0650,type_surface_standard_variable(name='surf_direct_downward_irradiance_0650_nm'))
      call self%register_dependency(self%id_Ed_0_0675,type_surface_standard_variable(name='surf_direct_downward_irradiance_0675_nm'))
      call self%register_dependency(self%id_Ed_0_0700,type_surface_standard_variable(name='surf_direct_downward_irradiance_0700_nm'))
      call self%register_dependency(self%id_Ed_0_0725,type_surface_standard_variable(name='surf_direct_downward_irradiance_0725_nm'))
      call self%register_dependency(self%id_Ed_0_0775,type_surface_standard_variable(name='surf_direct_downward_irradiance_0775_nm'))
      call self%register_dependency(self%id_Ed_0_0850,type_surface_standard_variable(name='surf_direct_downward_irradiance_0850_nm'))
      call self%register_dependency(self%id_Ed_0_0950,type_surface_standard_variable(name='surf_direct_downward_irradiance_0950_nm'))
      call self%register_dependency(self%id_Ed_0_1050,type_surface_standard_variable(name='surf_direct_downward_irradiance_1050_nm'))
      call self%register_dependency(self%id_Ed_0_1150,type_surface_standard_variable(name='surf_direct_downward_irradiance_1150_nm'))
      call self%register_dependency(self%id_Ed_0_1250,type_surface_standard_variable(name='surf_direct_downward_irradiance_1250_nm'))
      call self%register_dependency(self%id_Ed_0_1350,type_surface_standard_variable(name='surf_direct_downward_irradiance_1350_nm'))
      call self%register_dependency(self%id_Ed_0_1450,type_surface_standard_variable(name='surf_direct_downward_irradiance_1450_nm'))
      call self%register_dependency(self%id_Ed_0_1550,type_surface_standard_variable(name='surf_direct_downward_irradiance_1550_nm'))
      call self%register_dependency(self%id_Ed_0_1650,type_surface_standard_variable(name='surf_direct_downward_irradiance_1650_nm'))
      call self%register_dependency(self%id_Ed_0_1750,type_surface_standard_variable(name='surf_direct_downward_irradiance_1750_nm'))
      call self%register_dependency(self%id_Ed_0_1900,type_surface_standard_variable(name='surf_direct_downward_irradiance_1900_nm'))
      call self%register_dependency(self%id_Ed_0_2200,type_surface_standard_variable(name='surf_direct_downward_irradiance_2200_nm'))
      call self%register_dependency(self%id_Ed_0_2900,type_surface_standard_variable(name='surf_direct_downward_irradiance_2900_nm'))
      call self%register_dependency(self%id_Ed_0_3700,type_surface_standard_variable(name='surf_direct_downward_irradiance_3700_nm'))
      call self%register_dependency(self%id_Es_0_0250,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0250_nm'))
      call self%register_dependency(self%id_Es_0_0325,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0325_nm'))
      call self%register_dependency(self%id_Es_0_0350,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0350_nm'))
      call self%register_dependency(self%id_Es_0_0375,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0375_nm'))
      call self%register_dependency(self%id_Es_0_0400,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0400_nm'))
      call self%register_dependency(self%id_Es_0_0425,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0425_nm'))
      call self%register_dependency(self%id_Es_0_0450,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0450_nm'))
      call self%register_dependency(self%id_Es_0_0475,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0475_nm'))
      call self%register_dependency(self%id_Es_0_0500,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0500_nm'))
      call self%register_dependency(self%id_Es_0_0525,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0525_nm'))
      call self%register_dependency(self%id_Es_0_0550,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0550_nm'))
      call self%register_dependency(self%id_Es_0_0575,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0575_nm'))
      call self%register_dependency(self%id_Es_0_0600,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0600_nm'))
      call self%register_dependency(self%id_Es_0_0625,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0625_nm'))
      call self%register_dependency(self%id_Es_0_0650,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0650_nm'))
      call self%register_dependency(self%id_Es_0_0675,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0675_nm'))
      call self%register_dependency(self%id_Es_0_0700,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0700_nm'))
      call self%register_dependency(self%id_Es_0_0725,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0725_nm'))
      call self%register_dependency(self%id_Es_0_0775,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0775_nm'))
      call self%register_dependency(self%id_Es_0_0850,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0850_nm'))
      call self%register_dependency(self%id_Es_0_0950,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_0950_nm'))
      call self%register_dependency(self%id_Es_0_1050,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1050_nm'))
      call self%register_dependency(self%id_Es_0_1150,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1150_nm'))
      call self%register_dependency(self%id_Es_0_1250,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1250_nm'))
      call self%register_dependency(self%id_Es_0_1350,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1350_nm'))
      call self%register_dependency(self%id_Es_0_1450,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1450_nm'))
      call self%register_dependency(self%id_Es_0_1550,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1550_nm'))
      call self%register_dependency(self%id_Es_0_1650,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1650_nm'))
      call self%register_dependency(self%id_Es_0_1750,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1750_nm'))
      call self%register_dependency(self%id_Es_0_1900,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_1900_nm'))
      call self%register_dependency(self%id_Es_0_2200,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_2200_nm'))
      call self%register_dependency(self%id_Es_0_2900,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_2900_nm'))
      call self%register_dependency(self%id_Es_0_3700,type_surface_standard_variable(name='surf_diffuse_downward_irradiance_3700_nm'))
! END BLOCK2 python generated code
   end subroutine initialize

   subroutine do_column(self,_ARGUMENTS_VERTICAL_)
      class (type_ogs_bfm_light_spectral),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      integer  :: kk,nlev,l
      real(rk) :: dz,zenithA,mud
      real(rk) :: phy_a,phy_b,phy_bb
      real(rk) :: cdom_a
      real(rk) :: tot_a,tot_b,tot_bb
      real(rk) :: R6c,X1c,X2c,X3c
      real(rk) :: P1c,P2c,P3c,P4c
      real(rk) :: P1chl, P2chl, P3chl, P4chl
      real(rk) :: aph450, anap450
      real(rk) :: acdom250,acdom325,acdom400,acdom425,acdom450
      real(rk) :: Scdom350_500, Scdom250_325
      real(rk) :: bbp450, bbp550, bbp700
      real(rk) :: equis(self%nlt),ies(self%nlt)
      real(rk) :: rlamm
      real(rk) :: zgrid(cache%n+1)
      real(rk) :: a_array(cache%n, self%nlt)
      real(rk) :: b_array(cache%n, self%nlt)
      real(rk) :: bb_array(cache%n, self%nlt)
      real(rk) :: vd(cache%n, self%nlt)
      real(rk) :: E(3, cache%n+1, self%nlt)
      real(rk) :: E_ave(3, cache%n, self%nlt)
      real(rk) :: rd, rs, ru, vs, vu 
      real(rk) :: E_scalar(cache%n, self%nlt)
      real(rk) :: PAR_diatoms_array(cache%n)
      real(rk) :: PAR_flagellates_array(cache%n)
      real(rk) :: PAR_picophytoplankton_array(cache%n)
      real(rk) :: PAR_dinoflagellates_array(cache%n)
      real(rk) :: PAR_scalar_array(cache%n)
      real(rk) :: PAR_diatoms
      real(rk) :: PAR_flagellates
      real(rk) :: PAR_picophytoplankton
      real(rk) :: PAR_dinoflagellates
      real(rk) :: PAR_scalar

      _GET_HORIZONTAL_(self%id_zenithA,zenithA)   ! Zenith angle
      call getrmud(zenithA,mud) ! average cosine direct component in the water

!START BLOCK3
      _GET_SURFACE_(self%id_Ed_0_0250,Ed_0(1))
      _GET_SURFACE_(self%id_Ed_0_0325,Ed_0(2))
      _GET_SURFACE_(self%id_Ed_0_0350,Ed_0(3))
      _GET_SURFACE_(self%id_Ed_0_0375,Ed_0(4))
      _GET_SURFACE_(self%id_Ed_0_0400,Ed_0(5))
      _GET_SURFACE_(self%id_Ed_0_0425,Ed_0(6))
      _GET_SURFACE_(self%id_Ed_0_0450,Ed_0(7))
      _GET_SURFACE_(self%id_Ed_0_0475,Ed_0(8))
      _GET_SURFACE_(self%id_Ed_0_0500,Ed_0(9))
      _GET_SURFACE_(self%id_Ed_0_0525,Ed_0(10))
      _GET_SURFACE_(self%id_Ed_0_0550,Ed_0(11))
      _GET_SURFACE_(self%id_Ed_0_0575,Ed_0(12))
      _GET_SURFACE_(self%id_Ed_0_0600,Ed_0(13))
      _GET_SURFACE_(self%id_Ed_0_0625,Ed_0(14))
      _GET_SURFACE_(self%id_Ed_0_0650,Ed_0(15))
      _GET_SURFACE_(self%id_Ed_0_0675,Ed_0(16))
      _GET_SURFACE_(self%id_Ed_0_0700,Ed_0(17))
      _GET_SURFACE_(self%id_Ed_0_0725,Ed_0(18))
      _GET_SURFACE_(self%id_Ed_0_0775,Ed_0(19))
      _GET_SURFACE_(self%id_Ed_0_0850,Ed_0(20))
      _GET_SURFACE_(self%id_Ed_0_0950,Ed_0(21))
      _GET_SURFACE_(self%id_Ed_0_1050,Ed_0(22))
      _GET_SURFACE_(self%id_Ed_0_1150,Ed_0(23))
      _GET_SURFACE_(self%id_Ed_0_1250,Ed_0(24))
      _GET_SURFACE_(self%id_Ed_0_1350,Ed_0(25))
      _GET_SURFACE_(self%id_Ed_0_1450,Ed_0(26))
      _GET_SURFACE_(self%id_Ed_0_1550,Ed_0(27))
      _GET_SURFACE_(self%id_Ed_0_1650,Ed_0(28))
      _GET_SURFACE_(self%id_Ed_0_1750,Ed_0(29))
      _GET_SURFACE_(self%id_Ed_0_1900,Ed_0(30))
      _GET_SURFACE_(self%id_Ed_0_2200,Ed_0(31))
      _GET_SURFACE_(self%id_Ed_0_2900,Ed_0(32))
      _GET_SURFACE_(self%id_Ed_0_3700,Ed_0(33))
      _GET_SURFACE_(self%id_Es_0_0250,Es_0(1))
      _GET_SURFACE_(self%id_Es_0_0325,Es_0(2))
      _GET_SURFACE_(self%id_Es_0_0350,Es_0(3))
      _GET_SURFACE_(self%id_Es_0_0375,Es_0(4))
      _GET_SURFACE_(self%id_Es_0_0400,Es_0(5))
      _GET_SURFACE_(self%id_Es_0_0425,Es_0(6))
      _GET_SURFACE_(self%id_Es_0_0450,Es_0(7))
      _GET_SURFACE_(self%id_Es_0_0475,Es_0(8))
      _GET_SURFACE_(self%id_Es_0_0500,Es_0(9))
      _GET_SURFACE_(self%id_Es_0_0525,Es_0(10))
      _GET_SURFACE_(self%id_Es_0_0550,Es_0(11))
      _GET_SURFACE_(self%id_Es_0_0575,Es_0(12))
      _GET_SURFACE_(self%id_Es_0_0600,Es_0(13))
      _GET_SURFACE_(self%id_Es_0_0625,Es_0(14))
      _GET_SURFACE_(self%id_Es_0_0650,Es_0(15))
      _GET_SURFACE_(self%id_Es_0_0675,Es_0(16))
      _GET_SURFACE_(self%id_Es_0_0700,Es_0(17))
      _GET_SURFACE_(self%id_Es_0_0725,Es_0(18))
      _GET_SURFACE_(self%id_Es_0_0775,Es_0(19))
      _GET_SURFACE_(self%id_Es_0_0850,Es_0(20))
      _GET_SURFACE_(self%id_Es_0_0950,Es_0(21))
      _GET_SURFACE_(self%id_Es_0_1050,Es_0(22))
      _GET_SURFACE_(self%id_Es_0_1150,Es_0(23))
      _GET_SURFACE_(self%id_Es_0_1250,Es_0(24))
      _GET_SURFACE_(self%id_Es_0_1350,Es_0(25))
      _GET_SURFACE_(self%id_Es_0_1450,Es_0(26))
      _GET_SURFACE_(self%id_Es_0_1550,Es_0(27))
      _GET_SURFACE_(self%id_Es_0_1650,Es_0(28))
      _GET_SURFACE_(self%id_Es_0_1750,Es_0(29))
      _GET_SURFACE_(self%id_Es_0_1900,Es_0(30))
      _GET_SURFACE_(self%id_Es_0_2200,Es_0(31))
      _GET_SURFACE_(self%id_Es_0_2900,Es_0(32))
      _GET_SURFACE_(self%id_Es_0_3700,Es_0(33))

!END BLOCK3 python generated code
      
      kk=0
      zgrid(1)=0.0_rk

      _DOWNWARD_LOOP_BEGIN_
          kk = kk + 1
       _GET_(self%id_dz,dz)     ! Layer height (m)

         zgrid(kk+1)=zgrid(kk)+dz

       _GET_(self%id_P1chl,P1chl)
       _GET_(self%id_P2chl,P2chl)
       _GET_(self%id_P3chl,P3chl)
       _GET_(self%id_P4chl,P4chl)

       _GET_(self%id_P1c,P1c)
       _GET_(self%id_P2c,P2c)
       _GET_(self%id_P3c,P3c)
       _GET_(self%id_P4c,P4c)       
       
       _GET_(self%id_R6c,R6c)

       _GET_(self%id_X1c,X1c)
       _GET_(self%id_X2c,X2c)
       _GET_(self%id_X3c,X3c)

! Equations determining optical properties in relations to biogeochemical variables
          do l=1,self%nlt
             phy_a  = ac(1,l)*P1chl + ac(2,l)*P2chl + ac(3,l)*P3chl + ac(4,l)*P4chl
             phy_b  = bc(1,l)*P1c + bc(2,l)*P2c + bc(3,l)*P3c + bc(4,l)*P4c
!             phy_b  = bc(1,l)*P1chl + bc(2,l)*P2chl + bc(3,l)*P3chl + bc(4,l)*P4chl
             phy_bb = bc(1,l)*bbc(1,l)*P1c + bc(2,l)*bbc(2,l)*P2c + bc(3,l)*bbc(3,l)*P3c + bc(4,l)*bbc(4,l)*P4c
!             phy_bb = bc(1,l)*bbc(1,l)*P1chl + bc(2,l)*bbc(2,l)*P2chl + bc(3,l)*bbc(3,l)*P3chl + bc(4,l)*bbc(4,l)*P4chl
             cdom_a = acdom(1,l)*X1c  + acdom(2,l)*X2c  + acdom(3,l)*X3c 
             cdom_a = MAX(cdom_a, acdom_min(l))
             
! Need to add also cdom
             tot_a  =  aw(l) + phy_a  + apoc(l) * R6c + cdom_a
             tot_b  =  bw(l) + phy_b  + bpoc(l) * R6c 
             tot_bb = bbw(l) + phy_bb + bbpoc(l)* R6c 

             a_array(kk,l)  = tot_a
             b_array(kk,l)  = tot_b
             bb_array(kk,l) = tot_bb
          enddo

! IOPs observed for diagnostics
       acdom250 = MAX(acdom(1,1)*X1c + acdom(2,1)*X2c + acdom(3,1)*X3c, acdom_min(1))
       acdom325 = MAX(acdom(1,2)*X1c + acdom(2,2)*X2c + acdom(3,2)*X3c, acdom_min(2))
       acdom400 = MAX(acdom(1,5)*X1c + acdom(2,5)*X2c + acdom(3,5)*X3c, acdom_min(5))
       acdom425 = MAX(acdom(1,6)*X1c + acdom(2,6)*X2c + acdom(3,6)*X3c, acdom_min(6))
       acdom450 = MAX(acdom(1,7)*X1c + acdom(2,7)*X2c + acdom(3,7)*X3c, acdom_min(7))
       aph450   = ac(1,7)*P1chl + ac(2,7)*P2chl + ac(3,7)*P3chl + ac(4,7)*P4chl
       anap450  = apoc(7) * R6c
       bbp450 = (bc(1,7)*bbc(1,7)*P1c + bc(2,7)*bbc(2,7)*P2c + bc(3,7)*bbc(3,7)*P3c + bc(4,7)*bbc(4,7)*P4c) + bbpoc(7)*R6c
       bbp550 = (bc(1,11)*bbc(1,11)*P1c + bc(2,11)*bbc(2,11)*P2c + bc(3,11)*bbc(3,11)*P3c + bc(4,11)*bbc(4,11)*P4c) + bbpoc(11)*R6c
       bbp700 = (bc(1,17)*bbc(1,17)*P1c + bc(2,17)*bbc(2,17)*P2c + bc(3,17)*bbc(3,17)*P3c + bc(4,17)*bbc(4,17)*P4c) + bbpoc(17)*R6c

         _SET_DIAGNOSTIC_(self%id_acdom250, max(p_small,acdom250))            
         _SET_DIAGNOSTIC_(self%id_acdom325, max(p_small,acdom325))            
         _SET_DIAGNOSTIC_(self%id_acdom400, max(p_small,acdom400))            
         _SET_DIAGNOSTIC_(self%id_acdom425, max(p_small,acdom425))            
         _SET_DIAGNOSTIC_(self%id_acdom450, max(p_small,acdom450))
         _SET_DIAGNOSTIC_(self%id_aph450, max(p_small,aph450))              
         _SET_DIAGNOSTIC_(self%id_anap450, max(p_small,anap450))
         _SET_DIAGNOSTIC_(self%id_bbp450, max(p_small,bbp450))
         _SET_DIAGNOSTIC_(self%id_bbp550, max(p_small,bbp550))
         _SET_DIAGNOSTIC_(self%id_bbp700, max(p_small,bbp700))
         
       ! linear fit of ln-transformed aCDOM(l) against wavelength:
       ! between 350 and 500 nm (Babin et al 2013, Organelli et al 2014) 
       ! Organelli 2014 uses non-linear least-squares fit!
       ! between 250 and 325 nm (Catalá et al 2018, Galletti et al 2019) 
          do l=1,self%nlt
             rlamm = real(lam(l),8)     
             equis(l) = rlamm-self%lambda_aCDOM
             ies(l) = LOG(acdom(1,l)*X1c+acdom(2,l)*X2c+acdom(3,l)*X3c)
          enddo
          call linear_regression(equis(3:9),ies(3:9),7,acdom450,Scdom350_500) ! indexes for 350-500nm
          call linear_regression(equis(1:2),ies(1:2),7,acdom450,Scdom250_325) ! indexes for 250-325nm
       ! acdom450=EXP(acdom450)
         _SET_DIAGNOSTIC_(self%id_Scdom350_500, -Scdom350_500) 
         _SET_DIAGNOSTIC_(self%id_Scdom250_325, -Scdom250_325) 
       
     _DOWNWARD_LOOP_END_
  
     rd      = self%rd
     rs      = self%rs
     ru      = self%ru
     vd(:,:) = mud
     vs      = self%vs
     vu      = self%vu

     call solve_direct(cache%n+1, zgrid, cache%n, zgrid, self%nlt, a_array, b_array, bb_array, rd, rs, ru, vd, vs, vu, Ed_0, Es_0, E, E_ave)


! Scalar irradiance
     E_scalar(:,:)=E_ave(1,:,:)/vd + E_ave(2,:,:)/vs + E_ave(3,:,:)/vu

     PAR_diatoms_array(:)           = 0.0_rk
     PAR_flagellates_array(:)       = 0.0_rk
     PAR_picophytoplankton_array(:) = 0.0_rk
     PAR_dinoflagellates_array(:)   = 0.0_rk
     PAR_scalar_array(:)            = 0.0_rk

!     do l=1,self%nlt
     do l=5,17     
         PAR_diatoms_array(:)           = PAR_diatoms_array(:)           + (WtoQ(l) * ac_ps(1,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_flagellates_array(:)       = PAR_flagellates_array(:)       + (WtoQ(l) * ac_ps(2,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_picophytoplankton_array(:) = PAR_picophytoplankton_array(:) + (WtoQ(l) * ac_ps(3,l) * E_scalar(:,l)) * SEC_PER_DAY
         PAR_dinoflagellates_array(:)   = PAR_dinoflagellates_array(:)   + (WtoQ(l) * ac_ps(4,l) * E_scalar(:,l)) * SEC_PER_DAY
     enddo

     do l=5,17
         PAR_scalar_array(:)            = PAR_scalar_array(:) + (E_scalar(:,l) * WtoQ(l)) * SEC_PER_DAY
     enddo

! AOPs observed for diagnostics     
!    write(*,*) 'Rrs= ', E(3,1,7)/(Ed_0(7)+Es_0(7))

     _HORIZONTAL_LOOP_BEGIN_
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs400,E(3,1,5)/max(p_small,(E(1,1,5)+E(2,1,5))))
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs425,E(3,1,6)/max(p_small,(E(1,1,6)+E(2,1,6))))
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs450,E(3,1,7)/max(p_small,(E(1,1,7)+E(2,1,7)))) 
!    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs450,E(3,1,7)/(Ed_0(7)+Es_0(7))) 
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs475,E(3,1,8)/max(p_small,(E(1,1,8)+E(2,1,8))))
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs500,E(3,1,9)/max(p_small,(E(1,1,9)+E(2,1,9))))
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs525,E(3,1,10)/max(p_small,(E(1,1,10)+E(2,1,10)))) 
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs550,E(3,1,11)/max(p_small,(E(1,1,11)+E(2,1,11)))) 
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs575,E(3,1,12)/max(p_small,(E(1,1,12)+E(2,1,12)))) 
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Rrs675,E(3,1,16)/max(p_small,(E(1,1,16)+E(2,1,16)))) 

     write(*,*) 'Z9= ', zgrid(26)
     
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd375,-LOG(max(p_small,(E(1,26,4)+E(2,26,4))/(max(p_small,E(1,1,4)+E(2,1,4)))))/9.05_rk)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd400,-LOG(max(p_small,(E(1,26,5)+E(2,26,5))/(max(p_small,E(1,1,5)+E(2,1,5)))))/9.05_rk)     
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd425,-LOG(max(p_small,(E(1,26,6)+E(2,26,6))/(max(p_small,E(1,1,6)+E(2,1,6)))))/9.05_rk)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd475,-LOG(max(p_small,(E(1,26,8)+E(2,26,8))/(max(p_small,E(1,1,8)+E(2,1,8)))))/9.05_rk)
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_kd500,-LOG(max(p_small,(E(1,26,9)+E(2,26,9))/(max(p_small,E(1,1,9)+E(2,1,9)))))/9.05_rk)     
     
      _HORIZONTAL_LOOP_END_

      kk=0

      _DOWNWARD_LOOP_BEGIN_

          kk = kk + 1

         _SET_DIAGNOSTIC_(self%id_par_dia, max(p_small,PAR_diatoms_array(kk)))                  
         _SET_DIAGNOSTIC_(self%id_par_flag,max(p_small,PAR_flagellates_array(kk)))            
         _SET_DIAGNOSTIC_(self%id_par_pico,max(p_small,PAR_picophytoplankton_array(kk)))
         _SET_DIAGNOSTIC_(self%id_par_dino,max(p_small,PAR_dinoflagellates_array(kk)))          
         _SET_DIAGNOSTIC_(self%id_PAR_tot, max(p_small,PAR_scalar_array(kk)))          

     _DOWNWARD_LOOP_END_

   end subroutine do_column

end module
