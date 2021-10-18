MODULE trcsms_my_trc
   !!======================================================================
   !!                         ***  MODULE trcsms_my_trc  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_sms_my_trc       : MY_TRC model main routine
   !! trc_sms_my_trc_alloc : allocate arrays specific to MY_TRC sms
   !!----------------------------------------------------------------------
#include "domzgr_substitute.h90"
#include "do_loop_substitute.h90"
#include"cppdefs.h"
   ! NEMO
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc
   USE trcwri_my_trc,  ONLY: mapping_diags
   ! BFM
   USE global_mem, ONLY: LOGUNIT, bfm_lwp, RLEN, ZERO

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_my_trc       ! called by trcsms.F90 module

   ! Defined HERE the arrays specific to MY_TRC sms and ALLOCATE them in trc_sms_my_trc_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_my_trc.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_my_trc( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_my_trc  ***
      !!
      !! ** Purpose :   main routine of MY_TRC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER ::   ji, jj, jn   ! dummy loop index
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrmyt
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_my_trc')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_my_trc:  BFM ecosystem dynamics'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
     
      IF (kt == nit000) CALL mapping_diags()

      CALL update_bgc_forcings(kt)

      ! add here the call to BGC model
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         !
         ! cycle over land points
         IF (tmask(ji,jj,1) .EQ. 0) CYCLE 
         !
         ! Compose environmental Forcings
         CALL environmental_conditions(kt, Kmm, ji, jj)
         !
         ! Get BGC state variables
         ! call wrapper to execute BFM of functions

      END_2D 

      !
      IF( ln_timing )   CALL timing_stop('trc_sms_my_trc')
      !
   END SUBROUTINE trc_sms_my_trc


   SUBROUTINE environmental_conditions(kt, Kmm, ji, jj)
      !
      ! NEMO
      USE par_my_trc, ONLY: bottom_level
      USE oce_trc,    ONLY: jp_tem, jp_sal, ts, rhd, wndm, fr_i, gdept_0, gphit, tmask
      USE par_oce,    ONLY: jpk
      USE phycst,     ONLY: rho0
      USE eosbn2,     ONLY: neos
      USE sbcapr,     ONLY: apr
      ! BFM
      USE mem,        ONLY: xEPS, ESS, ETW, ESW, EWIND,    &
                            Depth, EIR, ERHO, EICE, EPR, NO_BOXES
      USE api_bfm,    ONLY: SRFindices, BOTindices
      USE sw_tool,    ONLY: sw_t_from_pt, gsw_p_from_z
#ifdef INCLUDE_PELCO2
      USE mem_CO2,    ONLY: AtmCO2, AtmSLP, patm3d
#endif
      !
      INTEGER, INTENT(in) :: kt   ! ocean time-step index
      INTEGER, INTENT(in) :: Kmm  ! time level indices
      INTEGER, INTENT(in) :: ji, jj ! i and j coordinate from 3D arrays
      INTEGER :: jk, bot
      !
      SRFindices(1) = 1
      BOTindices(1) = bottom_level(ji, jj)
      bot = BOTindices(1)
      !
      Depth(:) =  gdept(ji,jj,:,Kmm)
      EPR(:)   = gsw_p_from_z(-Depth(:), gphit(ji,jj))
      ! Environmental conditions
      ETW(:)   =  ts(ji,jj,:,jp_tem,Kmm)
      ESW(:)   =  ts(ji,jj,:,jp_sal,Kmm)
      ERHO(:)  =  (rhd(ji,jj,:) + 1._RLEN) * rho0 * tmask(ji,jj,:)
      EWIND(:) =  wndm(ji,jj)
      EICE(:)  =  fr_i(ji,jj)
      AtmSLP%fnow = apr(ji,jj)
      ! broadcast surface atm pressure over the water column (NO_BOXES)
      patm3d = AtmSLP%fnow
      !LEVEL1 'patm3d', patm3d(:)
#ifdef INCLUDE_PELCO2
      AtmCO2%fnow = atm_co2(ji,jj)
      !LEVEL1 'co2', AtmCO2%fnow, 'slp', AtmSLP%fnow
#endif
      ! convert to in-situ temperature and practical salinity
      ! TEOS-10 conversions not yet available
      if ( neos .eq. 0 ) then
         ETW(1:bot)  = sw_t_from_pt(ESW(1:bot),ETW(1:bot),EPR(1:bot),EPR(1:bot)*0.)
      endif

   END SUBROUTINE environmental_conditions


   SUBROUTINE update_bgc_forcings(kt)
      !
      USE fldread,       ONLY: fld_read
      !
      USE time,          ONLY: bfmtime
      USE sbcapr,        ONLY: apr
      USE SystemForcing, ONLY: FieldRead
      USE mem_CO2,       ONLY: AtmSLP, patm3d
#ifdef INCLUDE_PELCO2
      USE mem,        ONLY: ppO3c, ppO3h, ppN6r
      USE mem_CO2,    ONLY: AtmCO20, AtmCO2
      USE trcbc,      ONLY: sf_trcsbc, n_trc_indsbc
#endif
#ifdef INCLUDE_PELFE
      USE mem_PelChem,   ONLY: p_rDust
#endif
      !
      INTEGER, INTENT(in) :: kt   ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: irondep
      REAL(wp)   :: dustsink, dustinp
      ! 
      ! Update bfm internal time
      bfmtime%stepnow  = kt
      !
#ifdef INCLUDE_PELCO2
      !
      ! CO2 atmospheric mixing ratio
      SELECT CASE ( AtmCO2%init )
        CASE (1) ! Read timeseries in BFM
           call FieldRead(AtmCO2)
           atm_co2(:,:) = AtmCO2%fnow(1)
        !   LEVEL1 'bfmtime', bfmtime%stepnow, ' co2', AtmCO2%fnow(1)
        !CASE (2) ! Read Boundary Conditions using NEMO fldread (namtrc_bc)
        !   n = n_trc_indsbc(ppO3c)
        !   AtmCO2%fnow = pack( sf_trcsbc(n)%fnow(:,:,1),SRFmask(:,:,1) )
        !CASE (3) should not be needed. maybe for CCSMCOUPLED if still issue with kt < 10 in NEMO coupling
      END SELECT
#endif
      !
      ! Atmospheric sea level pressure
      SELECT CASE ( AtmSLP%init )
        CASE (1) ! Read timeseries in BFM
           call FieldRead(AtmSLP)
           apr(:,:) = AtmSLP%fnow(1)
        !CASE (2) ! Read Boundary Conditions using NEMO fldread (namtrc_bc)
        !   n = n_trc_indsbc(ppO3h)
        !   AtmSLP%fnow = pack( sf_trcsbc(n)%fnow(:,:,1),SRFmask(:,:,1) )
        !CASE (3) should not be needed. maybe for CCSMCOUPLED if still issue with kt < 10 in NEMO coupling
      END SELECT

      ! print control on received fluxes
      !IF ( (kt-nit000)<20 .OR. MOD(kt,200)==0 .OR. kt==nitend) THEN
      !   write(LOGUNIT,'(a,i14,a,f10.4,a,f10.3,a,f10.3)') 'envforcing_bfm - Step ', kt,  &
      !      '  Air_xCO2: ',maxval(AtmCO2%fnow), '  SLP:', maxval(AtmSLP%fnow),   &
      !      ' EIR: ',maxval(EIR)
      !ENDIF
      !
#ifdef INCLUDE_PELFE
      ! Iron dust deposition below the surface
      irondep = ZERO
      IF (p_rDust > 0. ) THEN
         CALL fld_read( kt, 1, sf_dust )
         ! dust sink length scale (1/m): dust assumed to be 0.01% (1/d) , sink speed p_rDust (m/d)
         dustsink =  0.0001_wp * 1. / p_rDust
         ! dust (g/day/m2) - > (umolFe/day/m2)
         dustinp = 0.035 * 1.e06 / 55.85
      ENDIF
#endif


   END SUBROUTINE update_bgc_forcings

   !!======================================================================

END MODULE trcsms_my_trc
