MODULE trcini_my_trc
   !!======================================================================
   !!                         ***  MODULE trcini_my_trc  ***
   !! TOP :   initialisation of the MY_TRC tracers
   !!======================================================================
   !! History :        !  2007  (C. Ethe, G. Madec) Original code
   !!                  !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_ini_my_trc   : MY_TRC model initialisation
   !!----------------------------------------------------------------------
#include "do_loop_substitute.h90"
#include "domzgr_substitute.h90"
#include "cppdefs.h"
   ! NEMO
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_my_trc
   USE trcnam_my_trc     ! MY_TRC SMS namelist
   USE lib_fortran,    ONLY: glob_sum
   USE sbcapr,     ONLY: apr
   USE trcopt,     ONLY: trc_opt_alloc, trc_opt_ini, parlux
   !! BFM
   USE constants,  ONLY: SEC_PER_DAY
   USE mem_param,  ONLY: p_atm0
   USE mem,        ONLY: D3STATE, NO_D3_BOX_STATES, D2STATE_BEN, NO_D2_BOX_STATES_BEN, iiPhytoPlankton
   USE mem_PAR,    ONLY: p_PAR, ChlAttenFlag, LightLocationFlag
   USE global_mem, ONLY: RLEN, ONE, ZERO, LOGUNIT, bfm_lwp
   USE api_bfm,    ONLY: bfm_init, in_rst_fname, InitVar, var_names, stBenStateS
   USE mem_globalfun, ONLY: analytical_ic
   USE init_var_bfm_local, ONLY: ini_organic_quotas
   USE mem_PelSinkSet, ONLY: sink_vars
#ifdef INCLUDE_PELFE
   USE mem_PelChem,   ONLY: p_rN7fsed
#endif
#ifdef INCLUDE_SEAICE
   USE api_bfm,    ONLY: stIceStateS
   USE mem,        ONLY: D2STATE_ICE, NO_D2_BOX_STATES_ICE
#endif
#ifdef INCLUDE_PELCO2
   USE mem_CO2,    ONLY: AtmCO20, AtmCO2, AtmSLP
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_my_trc, log_bgc_stats   ! called by trcini.F90 module
   PUBLIC   scale_by_density                ! called by trcsms_my_trc.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_my_trc.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_my_trc( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_my_trc  ***  
      !!
      !! ** Purpose :   initialization for MY_TRC model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER :: ji, jj, jk, jn, jt
      !
      IF( ln_timing )   CALL timing_start('trc_ini_my_trc')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_ini_my_trc : initialize BFM tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      IF( trc_ini_my_trc_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_my_trc: unable to allocate MY_TRC arrays' )
      !
      LEVEL1 ' '
      LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
      LEVEL1 '              BFM INITIALIZATION               '
      LEVEL1 ' '

      ! zeroing support arrays
      !-------------------------------------------------------
      chl_a = ZERO
      sink_rates = ZERO
      ph = ZERO

      ! Initialization from ICs
      !-------------------------------------------------------
      IF (bfm_init == 0) THEN
         ! PELAGIC
         DO jn = 1, NO_D3_BOX_STATES

            Initvar(jn)%varname=var_names(jn)

            SELECT CASE (InitVar(jn) % init)

            CASE (0) ! Homogeneous IC
               InitVar(jn)%unif = minval( D3STATE(:,jn) )

            CASE (1) ! Analytical IC
               ! Check consistency of input z1 and z2
               IF ( InitVar(jn)%anz2 .NE. ZERO .AND. InitVar(jn)%anz2 .LT. InitVar(jn)%anz1 ) THEN
                  STDERR 'ERROR: z2 < z1 in analytical profile creation for tracer ', TRIM(var_names(jn))
                  STDERR '     Check settings in bfm_init_nml.'
                  STOP
               ENDIF
               ! gdept_1d contains the model depth
               D3STATE(:,jn) = analytical_ic( gdept_1d, InitVar(jn)%anz1, &
                       InitVar(jn)%anv1, InitVar(jn)%anz2, InitVar(jn)%anv2 )

            CASE (2) ! IC from file (use nemo trcdta)
               InitVar(jn)%filename='file in namtrc_dta'
            END SELECT
         ENDDO

         ! Initialize internal constitutents quota of functional groups
         !-------------------------------------------------------
         call ini_organic_quotas()

         ! print initial state settings #TODO move this before stepping (here we miss info from trcdta)
         !-------------------------------------------------------
         IF (bfm_lwp) THEN
            WRITE(LOGUNIT,157) 'Init', 'Unif', 'Filename    ', 'Var', 'Anal. Z1', 'Anal. V1', &
                 'Anal. Z2', 'Anal. V2', 'OBC', 'SBC', 'CBC', 'RHO'
            DO jn = 1,NO_D3_BOX_STATES
                 WRITE(LOGUNIT, 158) InitVar(jn)
            !     write(LOGUNIT,*) var_names(m), D3STATE(m,:)
            ENDDO
         ENDIF

!TODO move this to trcsms
!#ifdef INCLUDE_PELCO2
!      ! Scale DIC and ALK from m(g/mol)/kg to m(g/mol)/m3 using NEMO initial in situ density
!      if (bfm_lwp) write(LOGUNIT,*)
!      if (bfm_lwp) write(LOGUNIT,*) 'trc_ini_bfm: Scale DIC and ALK 3D fields using NEMO in-situ density'
!      D3STATE(ppO3c,:) = D3STATE(ppO3c,:) * ERHO(:)
!      D3STATE(ppO3h,:) = D3STATE(ppO3h,:) * ERHO(:)
!#endif

         ! Transfer BFM ICs to passive tracer arrays
         !-------------------------------------------------------
         ! PELAGIC
         DO jn = 1, NO_D3_BOX_STATES
             DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                 tr(ji,jj,:,var_map(jn),Kmm) = D3STATE(:,jn)
             END_2D
             tr(:,:,:,var_map(jn),Kmm) = tr(:,:,:,var_map(jn),Kmm) * tmask(:,:,:)
         END DO
         ! BENTHIC
         DO jn = 1, NO_D2_BOX_STATES_BEN
            DO jk = 1, jpk_b
                tr_b(:,:,jk,jn) = tmask(:,:,1) * D2STATE_BEN(1,jn)
            ENDDO
         END DO
#ifdef INCLUDE_SEAICE
         ! SEAICE
         DO jn = 1, NO_D2_BOX_STATES_ICE
            DO jk = 1, jpk_i
                tr_i(:,:,jk,jn) = tmask(:,:,1) * D2STATE_ICE(1,jn)
            ENDDO
         END DO
#endif
      ENDIF

      ! Restart handled by nemo
      !-------------------------------------------------------
      IF (bfm_init == 1 ) THEN
         call ini_organic_quotas()
      ENDIF

      ! Assign variable names arrays
      !-------------------------------------------------------
      DO jn = 1, NO_D2_BOX_STATES_BEN
         ctrcnm_b(jn) = trim(var_names(jn + stBenStateS - 1))
      END DO
#ifdef INCLUDE_SEAICE
      DO jn = 1, NO_D2_BOX_STATES_ICE
         ctrcnm_i(jn) = trim(var_names(jn + stIceStateS - 1))
      END DO
#endif
         
      ! Initialize sea level pressure for gas exchange
      !-------------------------------------------------------
      IF ( .NOT. ALLOCATED(apr)) ALLOCATE(apr(jpi,jpj))
      IF (AtmSLP%init .NE. 3) apr = p_atm0
#ifdef INCLUDE_PELCO2
      ! Initialize CO2 air concentration
      !-------------------------------------------------------
      IF (AtmCO2%init .NE. 3) atm_co2 = AtmCO20
#endif

      ! Initialise DATA output netcdf file(s)
      !-------------------------------------------------------
      LEVEL1 ''
      IF ( bfm_iomput ) THEN
         LEVEL1 'BFM uses XIOS output system'
         ! Set var_ids values according to the XIOS file_def is done in trc_dia_bfm
      ELSE
         LEVEL1 'BFM output without XIOS not implemented' 
         STOP
      ENDIF 

#ifdef INCLUDE_PELFE
      ! Prescribe iron release from seabed sediments
      !-------------------------------------------------------
      CALL seabed_iron_release
#endif
      !
      IF (bfm_lwp) THEN
            LEVEL1 ' '
            LEVEL1 '         BFM INITIALIZATION ... DONE!          '
            LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
            LEVEL1 ' '
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_ini_my_trc')
      !
157 FORMAT(a4, 1x, a10  , 1x, a25, 1x, a3, 1x, a10  , 1x, a10  , 1x, a10  , 1x, a10  , 3x, a3, 3x, a3, 3x, a3, 3x, a3)
158 FORMAT(i4, 1x, E10.3, 1x, a25, 1x, a3, 1x, E10.3, 1x, E10.3, 1x, E10.3, 1x, E10.3, 3x, L3, 3x, L3, 3x, L3, 3x, L3)
      !
   END SUBROUTINE trc_ini_my_trc


   INTEGER FUNCTION trc_ini_my_trc_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_ini_my_trc_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( tr_b(jpi,jpj,jpk_b,jp_bgc_b) , ctrcnm_b(jp_bgc_b), &
                chl_a(jpi,jpj,jpk),   ph(jpi,jpj,jpk),             &
                sink_rates(jpi,jpj,jpk,sink_vars),                   &
#ifdef INCLUDE_SEAICE
                tr_i(jpi,jpj,jpk_i,jp_bgc_i) , ctrcnm_i(jp_bgc_i), &
#endif
               STAT = trc_ini_my_trc_alloc)
      !
      !
      IF( trc_ini_my_trc_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_ini_my_trc_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_ini_my_trc_alloc


   SUBROUTINE log_bgc_stats(kt, Kmm)
      !!----------------------------------------------------------------------
      !!                 ***  log_bgc_stats  ***
      !!
      !! ** Purpose :  save to log file summary stats of 3D pelagic BGC fields
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      USE api_bfm,    ONLY: var_names, stPelStateS, parallel_rank, NO_D3_BOX_STATES
      !
      INTEGER, INTENT(in) :: kt   ! ocean time-step index
      INTEGER, INTENT(in) :: Kmm  ! time level indices
      INTEGER :: jn, jk, jl
      REAL(wp) :: ztraf, zmin, zmax, zmean
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zvol
      !!----------------------------------------------------------------------

      LEVEL1 'Global statistics on tracer at beginning of step: ' , kt

      DO jk = 1, jpk
         zvol(:,:,jk) = e1e2t(:,:) * e3t(:,:,jk,Krhs) * tmask(:,:,jk)
      END DO

      DO jl = 1, NO_D3_BOX_STATES
         jn = var_map(jl)
         IF (jn > 0) THEN
            ztraf = glob_sum( 'log_bgc_stats', tr(:,:,:,jn,Kmm) * zvol(:,:,:) )
            zmin  = MINVAL( tr(:,:,:,jn,Kmm), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
            zmax  = MAXVAL( tr(:,:,:,jn,Kmm), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
            zmin  = MINVAL( tr(:,:,:,jn,Kmm), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
            zmax  = MAXVAL( tr(:,:,:,jn,Kmm), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
            IF( lk_mpp ) THEN
               CALL mpp_min( 'log_bgc_stats', zmin )      ! min over the global domain
               CALL mpp_max( 'log_bgc_stats', zmax )      ! max over the global domain
            END IF
            zmean  = ztraf / areatot

            IF (bfm_lwp) &
               WRITE(LOGUNIT,9000) jn, trim(var_names(stPelStateS+jl-1)), zmean, zmin, zmax
         END IF
      ENDDO

      LEVEL1 ''
      call FLUSH(LOGUNIT)

9000  FORMAT(' STAT tracer nb :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, '    max :',e18.10)

   END SUBROUTINE log_bgc_stats


   SUBROUTINE seabed_iron_release()
      !!----------------------------------------------------------------------
      !!                 ***  seabed_iron_release  ***
      !!
      !! ** Purpose : Prescribe the iron release from the seabed sediments
      !!
      !! ** Method  : create source array (fesed) with constant rate (p_rN7fsed) 
      !!            - constant relase at each bottom cell
      !!            - use depth varying mask as in Aumont & Bopp (2006)
      !!----------------------------------------------------------------------
      USE iom,        ONLY: iom_open,iom_get,iom_close
      USE par_my_trc, ONLY: bottom_level
      !
      INTEGER    :: ji, jj, jk, numfld
      REAL(RLEN) :: ztraf, zexp, zdexp
      LOGICAL    :: exists
      !!----------------------------------------------------------------------

      ALLOCATE (fesed(jpi,jpj,jpk))
      fesed = ZERO
      IF ( p_rN7fsed > 0. ) THEN
         LEVEL1 ''
         ! check if the mask file exists or create one only for bottom gridcells
         INQUIRE(FILE="bottom_fraction.nc", EXIST=exists)
         if ( exists ) then
            LEVEL1 'trc_ini_my_trc: Read fraction mask from bottom_fraction.nc (varname: btmfrac)'
            CALL iom_open ( TRIM( "bottom_fraction.nc" ), numfld )
            CALL iom_get  ( numfld, 1, TRIM( "btmfrac" ), fesed(:,:,:), 1 )
            CALL iom_close( numfld )
         ELSE
            LEVEL1 'trc_ini_my_trc: apply input to seabottom gridcells'
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls)
                 jk = bottom_level(ji,jj)
                 fesed(ji,jj,jk) = ONE
            END_2D
         ENDIF
         ! iron release dependence on depth (metamodel of Middelburg et al.,1996)
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk )
              zexp  = MIN( 8.,( gdept(ji,jj,jk,Kmm) / 500. )**(-1.5) )
              zdexp = -0.9543 + 0.7662 * LOG( zexp ) - 0.235 * LOG( zexp )**2
              fesed(ji,jj,jk) = fesed(ji,jj,jk) * MIN( 1., EXP( zdexp ) / 0.5 )
         END_3D
         ! total supply (TODO check why this is killing intialization time)
         ztraf = glob_sum('fesed', fesed * spread(e1e2t,3,jpk) * p_rN7fsed * 365. * 1.e-15 * tmask )
         LEVEL1 '  Constant Iron flux from sediments: ', p_rN7fsed
         LEVEL1 '  Total iron load from Sediment [Gmol/y] : ', ztraf
         ! Iron flux into pelagic (umolFe/m2/d to umolFe/m3/s)
         fesed = fesed * p_rN7fsed / ( SEC_PER_DAY * e3t_0(:,:,:) )
      ENDIF

   END SUBROUTINE seabed_iron_release


   SUBROUTINE scale_by_density(Kbb, Kmm)
      !!----------------------------------------------------------------------
      !!                 ***  scale_by_density  ***
      !!
      !! ** Purpose : Scale the value of a tracer by the seawater density
      !!
      !! ** Method  : scale bgc tracers by the value of ocean model density
      !!              based on the settings of Initvar(:)%rho 
      !!----------------------------------------------------------------------
      USE phycst,     ONLY: rho0
      !
      INTEGER, INTENT(in) :: Kbb, Kmm  ! time level indices
      INTEGER    :: ji, jj, jn
      REAL(RLEN) :: ztraf, zexp, zdexp
      LOGICAL    :: exists
      !!----------------------------------------------------------------------
      DO jn = 1, NO_D3_BOX_STATES
          IF ( Initvar(jn)%rho ) THEN 
             LEVEL1 ' Scale by density : ', trim(var_names(jn))
             tr(:,:,:,var_map(jn),Kmm) = tr(:,:,:,var_map(jn),Kmm) * tmask(:,:,:) * (rhd(:,:,:) + 1._RLEN) * rho0
             tr(:,:,:,var_map(jn),Kbb) = tr(:,:,:,var_map(jn),Kmm)
          ENDIF
      END DO
      LEVEL1 ''

   END SUBROUTINE scale_by_density

   !!======================================================================
END MODULE trcini_my_trc
