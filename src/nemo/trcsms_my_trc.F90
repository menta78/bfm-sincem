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
   USE trcini_my_trc,  ONLY: log_bgc_stats, scale_by_density
   USE trcwri_my_trc,  ONLY: diags_mapping, diags_collect
   USE trcset,         ONLY: trc_set
   USE iom,            ONLY: iom_put
   ! BFM
   USE trcopt_bfm, ONLY: trc_opt_2bds, trc_opt_3bds
   USE global_mem, ONLY: LOGUNIT, bfm_lwp, RLEN, ZERO
   USE constants,  ONLY: SEC_PER_DAY
   USE time,       ONLY: bfmtime

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
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_my_trc')
      !
      IF (kt == nittrc000) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_sms_my_trc:  BFM ecosystem dynamics'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      ! Set diagnostic fields to be saved
      !-------------------------------------------------------
      IF (kt == nittrc000) CALL diags_mapping()
   
      ! Scale initial conditions by seawater density
      !-------------------------------------------------------
      IF (kt == nittrc000 .AND. .NOT. ln_rsttr ) CALL scale_by_density(Kbb, Kmm) 

      ! Save to BFM log global statistics of tracers
      !-------------------------------------------------------
      IF ( (kt-nit000)<100 .OR. MOD(kt,200)==0 .OR. kt==nitend ) &
         call log_bgc_stats(kt, Kmm)

      ! Update bfm internal time
      !-------------------------------------------------------
      bfmtime%stepnow  = kt

      ! Update forcing for bgc
      !-------------------------------------------------------
      CALL update_bgc_forcings(kt, Kbb, Kmm, Krhs)

      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

         ! cycle over land points
         !-------------------------------------------------------
         IF (tmask(ji,jj,1) .EQ. 0) CYCLE 

         ! BFM environmental conditions
         !-------------------------------------------------------
         CALL environmental_conditions(kt, Kmm, ji, jj)

         ! Execute BFM
         !-------------------------------------------------------
         CALL execute_bfm(kt, Kmm, Krhs, ji, jj)

      END_2D 

      ! Vertical Sinking/Settling
      !-------------------------------------------------------
      CALL trc_set(kt, Kmm, Krhs)

      !
      IF( ln_timing )   CALL timing_stop('trc_sms_my_trc')
      !
   END SUBROUTINE trc_sms_my_trc


   SUBROUTINE execute_bfm(kt, Kmm, Krhs, ji, jj)
      !!----------------------------------------------------------------------
      !!                 ***  execute_bfm  ***
      !!
      !! ** Purpose :  Fill BFM arrays of environmental conditions
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      ! NEMO
      ! BFM
      USE api_bfm,    ONLY: SRFindices, BOTindices
      USE mem,        ONLY: D3STATE, D3SOURCE, NO_D3_BOX_STATES, D3STATETYPE, & 
                            D2STATE_BEN, D2SOURCE_BEN, NO_D2_BOX_STATES_BEN
      USE mem_param,  ONLY: CalcTransportFlag, CalcBenthicFlag, CalcPelagicFlag
      USE mem_PelSinkSet, ONLY: sink_var_map, sink_vars, SINKD3STATE
#ifdef INCLUDE_SEAICE
      USE mem,        ONLY: D2STATE_ICE, D2SOURCE_ICE, NO_D2_BOX_STATES_ICE
#endif
#ifdef INCLUDE_PELCO2
      USE mem,        ONLY: D3DIAGNOS, pppH
#endif
      !
      INTEGER, INTENT(in) :: kt         ! ocean time-step index
      INTEGER, INTENT(in) :: Kmm, Krhs  ! time level indices
      INTEGER, INTENT(in) :: ji, jj     ! i and j coordinate from 3D arrays
      INTEGER :: jn, bot
      REAL(wp) :: delt
      !!----------------------------------------------------------------------

      bot = BOTindices(1)
      delt = bfmtime%timestep

      ! Reset flux arrays
      !---------------------------------------------
      CALL ResetFluxes
      
      ! Fill BFM STATE arrays
      !---------------------------------------------
      DO jn = 1, NO_D3_BOX_STATES
         if (var_map(jn) > 0) &
            D3STATE(1:bot,jn) = tr(ji,jj,1:bot,var_map(jn),Kmm)
      END DO
      DO jn = 1, NO_D2_BOX_STATES_BEN
         D2STATE_BEN(:,jn) = tr_b(ji,jj,:,jn)
      END DO
#ifdef INCLUDE_SEAICE
      DO jn = 1, NO_D2_BOX_STATES_ICE
         D2STATE_ICE(:,jn) = tr_i(ji,jj,:,jn)
      END DO
#endif
#ifdef INCLUDE_PELCO2
      D3DIAGNOS(1:bot,pppH) = ph(ji,jj,1:bot)
#endif

      ! Compute Biogeochemical trends
      !---------------------------------------------
      call EcologyDynamics

      ! ODE solver for benthic and not transported pelagic
      ! & assign trend back main arrays
      !---------------------------------------------
      if (CalcPelagicFlag ) then
         do jn = 1, NO_D3_BOX_STATES
            if (var_map(jn) > 0) &
               tr(ji,jj,1:bot,var_map(jn),Krhs) = tr(ji,jj,1:bot,var_map(jn),Krhs) + D3SOURCE(1:bot,jn)
         end do
#ifdef INCLUDE_PELCO2
         ! ph for next iteration
         ph(ji,jj,1:bot) = D3DIAGNOS(1:bot,pppH)
#endif
      end if
      if (CalcBenthicFlag ) then
         DO jn = 1, NO_D2_BOX_STATES_BEN
            tr_b(ji,jj,:,jn) = tr_b(ji,jj,:,jn) + delt*D2SOURCE_BEN(:,jn)
         END DO
      end if
#ifdef INCLUDE_SEAICE
      DO jn = 1, NO_D2_BOX_STATES_ICE
         tr_i(ji,jj,:,jn) = tr_i(ji,jj,:,jn) + delt*D2SOURCE_ICE(:,jn)
      END DO
#endif

      ! Sinking & settling rates
      !---------------------------------------------
      DO jn = 1 , sink_vars 
         sink_rates(ji, jj, 1:bot,jn) = SINKD3STATE(sink_var_map(jn))%sedi(1:bot)
      ENDDO

      ! Diagnostics
      !---------------------------------------------
      CALL diags_collect(kt, ji, jj, bot)

   END SUBROUTINE execute_bfm


   SUBROUTINE environmental_conditions(kt, Kmm, ji, jj)
      !!----------------------------------------------------------------------
      !!                 ***  environmental_conditions  ***
      !!
      !! ** Purpose :  Fill BFM arrays of environmental conditions
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      ! NEMO
      USE par_my_trc, ONLY: bottom_level
      USE oce_trc,    ONLY: jp_tem, jp_sal, ts, rhd, wndm, fr_i, gdept_0, gphit, tmask
      USE par_oce,    ONLY: jpk
      USE phycst,     ONLY: rho0
      USE eosbn2,     ONLY: neos
      USE sbcapr,     ONLY: apr
      USE trc,        ONLY: etot
      USE trcopt,     ONLY: zeps
      ! BFM
      USE mem,        ONLY: xEPS, EIR, ESS, ETW, ESW, EWIND,    &
                            Depth, ERHO, EICE, EPR, NO_BOXES
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
      !!----------------------------------------------------------------------

      ! Mask and indeces
      !-------------------------------------------------------
      SRFindices(1) = 1
      BOTindices(1) = bottom_level(ji, jj)
      bot = BOTindices(1)

      ! Cell Tickness (Depth) & Water Pressure levels (EPR)
      !-------------------------------------------------------
      Depth(:) =  e3t(ji,jj,:,Kmm)
      EPR(:)   = gsw_p_from_z(-gdept(ji,jj,:,Kmm), gphit(ji,jj))

      ! Environmental conditions
      !-------------------------------------------------------
      ESW(:)   =  ts(ji,jj,:,jp_sal,Kmm)
      ETW(:)   =  ts(ji,jj,:,jp_tem,Kmm)
      ! convert to in-situ temperature and practical salinity (TEOS-10 not available)
      if ( neos .eq. 0 ) then
         ETW(1:bot)  = sw_t_from_pt(ESW(1:bot),ETW(1:bot),EPR(1:bot),EPR(1:bot)*0._RLEN)
      endif
      ERHO(:)  =  (rhd(ji,jj,:) + 1._RLEN) * rho0 * tmask(ji,jj,:)
      EWIND(:) =  wndm(ji,jj)
      EICE(:)  =  fr_i(ji,jj)
      !ETAUB missing, require field from nemo (e.g. taubot)

      ! Water column optics
      !-------------------------------------------------------
      EIR(:) = etot(ji,jj,:)
      xEPS(:) = zeps(ji,jj,:)

      ! Gas exchanges
      !-------------------------------------------------------
      AtmSLP%fnow(1) = apr(ji,jj)
      ! broadcast surface atm pressure over the water column
      patm3d = SPREAD(AtmSLP%fnow(1), DIM=1, NCOPIES=SIZE(patm3d))
#ifdef INCLUDE_PELCO2
      AtmCO2%fnow(1) = atm_co2(ji,jj)
      !LEVEL1 'co2', AtmCO2%fnow, 'slp', AtmSLP%fnow
#endif
      !LEVEL1 'ETW ', ETW(1), ', ESW ', ESW(1), ', ERHO ', ERHO(1), &
      !    ', EWIND ', EWIND(1), ', EICE ', EICE(1)

   END SUBROUTINE environmental_conditions


   SUBROUTINE update_bgc_forcings(kt, Kbb, Kmm, Krhs)
      !!----------------------------------------------------------------------
      !!                 ***  update_bgc_forcings  ***
      !!
      !! ** Purpose :  Update forcing fields of BGC processes
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      ! NEMO
      USE trcbc,         ONLY: n_trc_indsbc, sf_trcsbc, rf_trsfac
      USE sbcapr,        ONLY: apr
      USE trc,           ONLY: ln_trc_sbc
      ! BFM
      USE SystemForcing, ONLY: FieldRead
      USE mem_CO2,       ONLY: AtmSLP, patm3d
      USE mem_PAR,       ONLY: ChlAttenFlag
      USE mem,           ONLY: iiPhytoPlankton, iiC, iiL, ppPhytoPlankton, ppN1p, ppN5s
      USE mem_Phyto,     ONLY: p_qlcPPY
      USE constants,     ONLY: E2W
#ifdef INCLUDE_PELCO2
      USE mem,        ONLY: ppO3c, ppO3h
      USE mem_CO2,    ONLY: AtmCO20, AtmCO2
#endif
#ifdef INCLUDE_PELFE
      USE mem_PelChem,   ONLY: p_rDust
      USE mem,           ONLY: ppN7f
#endif
      !
      INTEGER, INTENT(in) :: kt   ! ocean time-step index
      INTEGER, INTENT(in) :: Kbb, Kmm, Krhs ! time level indices
      INTEGER             :: jk, jl, jn
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: dustdep
      REAL(wp), DIMENSION(jpi,jpj)     :: dustinp
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: par_r, par_g, par_b
      REAL(wp)   :: dustsink 
      REAL(wp) :: r_e2w = 1._RLEN / E2W
      !!----------------------------------------------------------------------

      ! Water column optics
      !-------------------------------------------------------
      IF ( ChlAttenFlag == 1 ) THEN
         CALL trc_opt_2bds(Kmm)
      ELSEIF ( ChlAttenFlag == 2 ) THEN
         CALL trc_opt_3bds(Kmm)
      ENDIF
      ! convert from Watt to Einstein
      etot = etot * r_e2w

      ! Atmospheric sea level pressure
      !-------------------------------------------------------
      SELECT CASE ( AtmSLP%init )
        CASE (1) ! Read timeseries in BFM
           CALL FieldRead(AtmSLP)
           apr(:,:) = AtmSLP%fnow(1)
        !CASE (3) should not be needed. maybe for CCSMCOUPLED if still issue with kt < 10 in NEMO coupling
      END SELECT

#ifdef INCLUDE_PELCO2
      ! CO2 atmospheric mixing ratio
      !-------------------------------------------------------
      SELECT CASE ( AtmCO2%init )
        CASE (1) ! Read timeseries in BFM
           CALL FieldRead(AtmCO2)
           atm_co2(:,:) = AtmCO2%fnow(1)
        !   LEVEL1 'bfmtime', bfmtime%stepnow, ' co2', AtmCO2%fnow(1)
        !CASE (3) should not be needed. maybe for CCSMCOUPLED if still issue with kt < 10 in NEMO coupling
      END SELECT
#endif

#ifdef INCLUDE_PELFE
      ! Iron flux from sediments (time-invariant)
      !-------------------------------------------------------
      jl = var_map(ppN7f)
      tr(:,:,:,jl,Krhs) = tr(:,:,:,jl,Krhs) + fesed
#endif

      ! Nutrients from dust deposition below the surface (surface layer applied in trcbc)
      !-------------------------------------------------------
      IF (p_rDust > 0.) THEN
         ! dust sink length scale (1/m): dust assumed to be 0.01% (1/d) , sink speed p_rDust (m/d)
         dustsink =  0.01_wp * 1. / p_rDust
#ifdef INCLUDE_PELFE
         ! Iron
         IF (ln_trc_sbc(jl)) THEN
           dustdep = ZERO
           ! surface input (rf_trsfac = Solub(0.01) * FeinDust(0.035) * g->ug (1e6) / Fe g->mol (55.85))
           jn = n_trc_indsbc(jl)
           dustinp =  rf_trsfac(jn) * sf_trcsbc(jn)%fnow(:,:,1)
           DO jk = 2, jpkm1
              dustdep(:,:,jk) = dustinp(:,:) * dustsink / SEC_PER_DAY * EXP( -gdept(:,:,jk,Kmm) / 540. )
              tr(:,:,jk,jl,Krhs) = tr(:,:,jk,jl,Krhs) + dustdep(:,:,jk)
           END DO
         ENDIF
#endif
         ! Phosphate
         jl = var_map(ppN1p)
         IF (ln_trc_sbc(jl)) THEN
           jn = n_trc_indsbc(jl)
           dustinp =  rf_trsfac(jn) * sf_trcsbc(jn)%fnow(:,:,1)
           DO jk = 2, jpkm1
              dustdep(:,:,jk) = dustinp(:,:) * dustsink / SEC_PER_DAY * EXP( -gdept(:,:,jk,Kmm) / 540. )
              tr(:,:,jk,jl,Krhs) = tr(:,:,jk,jl,Krhs) + dustdep(:,:,jk)
           END DO
         ENDIF
         ! silicate
         jl = var_map(ppN5s)
         IF (ln_trc_sbc(jl)) THEN
           jn = n_trc_indsbc(jl)
           dustinp =  rf_trsfac(jn) * sf_trcsbc(jn)%fnow(:,:,1)
           DO jk = 2, jpkm1
              dustdep(:,:,jk) = dustinp(:,:) * dustsink / SEC_PER_DAY * EXP( -gdept(:,:,jk,Kmm) / 540. )
              tr(:,:,jk,jl,Krhs) = tr(:,:,jk,jl,Krhs) + dustdep(:,:,jk)
           END DO
         ENDIF
      ENDIF

      ! print control on received fluxes
      !-------------------------------------------------------
      !IF ( (kt-nittrc000)<20 .OR. MOD(kt,200)==0 .OR. kt==nitend) THEN
      !   write(LOGUNIT,'(a,i14,a,f10.4,a,f10.3,a,f10.3)') 'envforcing_bfm - Step ', kt,  &
      !      '  Air_xCO2: ',maxval(AtmCO2%fnow), '  SLP:', maxval(AtmSLP%fnow),   &
      !      ' EIR: ',maxval(EIR)
      !ENDIF

   END SUBROUTINE update_bgc_forcings


   !!======================================================================

END MODULE trcsms_my_trc
