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
   USE iom,        ONLY: iom_open,iom_get,iom_close
   USE sbcapr,     ONLY: apr
   USE trcopt,     ONLY: trc_opt_alloc, trc_opt_ini, parlux
   !! BFM
   USE constants,  ONLY: SEC_PER_DAY
   USE mem_param,  ONLY: p_atm0
   USE mem,        ONLY: D3STATE, NO_D3_BOX_STATES, D2STATE_BEN, NO_D2_BOX_STATES_BEN, iiPhytoPlankton, Depth
   USE mem_PAR,    ONLY: p_PAR, ChlAttenFlag, LightLocationFlag
   USE global_mem, ONLY: RLEN, ZERO, LOGUNIT, SkipBFMCore, bfm_lwp
   USE api_bfm,    ONLY: bfm_init, in_rst_fname, InitVar, var_names, stBenStateS, SEAmask, ZEROS
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

   PUBLIC   trc_ini_my_trc   ! called by trcini.F90 module

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
      INTEGER :: ji, jj, jk, jn, jt, numfld, ierr
      REAL(RLEN) :: ztraf, zexp, zdexp
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

      Depth(:) =  pack(gdept(:,:,:,Kmm), SEAmask)

      ! Restart handled by nemo
      !-------------------------------------------------------
      if (bfm_init == 1) then 
         ln_rsttr = .true.
         cn_trcrst_in = in_rst_fname

      ! Initialization 
      !-------------------------------------------------------
      else if (bfm_init == 0) then
         ! PELAGIC
         do jn = 1,NO_D3_BOX_STATES

            Initvar(jn)%varname=var_names(jn)

            select case (InitVar(jn) % init)

            case (0) ! Homogeneous IC
               InitVar(jn)%unif = minval( D3STATE(jn,:) )

            case (1) ! Analytical IC
               ! Check consistency of input z1 and z2
               if ( InitVar(jn)%anz2 .ne. ZERO .AND. InitVar(jn)%anz2 .lt. InitVar(jn)%anz1 ) then
                  STDERR 'ERROR: z2 < z1 in analytical profile creation for tracer ', TRIM(var_names(jn))
                  STDERR '     Check settings in bfm_init_nml.'
                  stop
               endif
               ! gdept_1d contains the model depth
               D3STATE(jn,:) = analytical_ic( Depth(:), InitVar(jn)%anz1, &
                       InitVar(jn)%anv1, InitVar(jn)%anz2, InitVar(jn)%anv2 )

            case (2) ! IC from file (use nemo trcdta)
               InitVar(jn)%filename='file in namtrc_dta'
            end select
         end do

         ! Initialize internal constitutents quota of functional groups
         !-------------------------------------------------------
         call ini_organic_quotas()

         ! print initial state settings #TODO move this before stepping (here we miss info from trcdta)
         !-------------------------------------------------------
         if (bfm_lwp) then
            write(LOGUNIT,157) 'Init', 'Unif', 'Filename    ', 'Var', 'Anal. Z1', 'Anal. V1', &
                 'Anal. Z2', 'Anal. V2', 'OBC', 'SBC', 'CBC'
            do jn = 1,NO_D3_BOX_STATES
                 write(LOGUNIT, 158) InitVar(jn)
            !     write(LOGUNIT,*) var_names(m), D3STATE(m,:)
            enddo
         endif

      end if

!TODO move this to trcsms
!#ifdef INCLUDE_PELCO2
!      ! Scale DIC and ALK from m(g/mol)/kg to m(g/mol)/m3 using NEMO initial in situ density
!      if (bfm_lwp) write(LOGUNIT,*)
!      if (bfm_lwp) write(LOGUNIT,*) 'trc_ini_bfm: Scale DIC and ALK 3D fields using NEMO in-situ density'
!      D3STATE(ppO3c,:) = D3STATE(ppO3c,:) * ERHO(:)
!      D3STATE(ppO3h,:) = D3STATE(ppO3h,:) * ERHO(:)
!#endif

      ! Transfer BFM ICs to passive tracers array
      !-------------------------------------------------------
      DO jn = 1, NO_D3_BOX_STATES
         tr(:,:,:,var_map(jn),Kmm) = unpack(D3STATE(jn,:), SEAmask, ZEROS)
      ENDDO


      ! Initialize sea level pressure for gas exchange
      !-------------------------------------------------------
      IF ( .NOT. ALLOCATED(apr)) ALLOCATE(apr(jpi,jpj))
      IF (AtmSLP%init .NE. 3) apr = p_atm0
#ifdef INCLUDE_PELCO2
      ! Initialize CO2 air concentration
      !-------------------------------------------------------
      IF (AtmCO2%init .NE. 3) atm_co2 = AtmCO20
#endif

      ! Zero out fields if cpu is off
      !-------------------------------------------------------
      IF( SkipBFMCore ) tr(:,:,:,:,Kmm) = ZERO


      ! Initialise DATA output netcdf file(s)
      !-------------------------------------------------------
      IF ( bfm_iomput ) THEN
         if (lwp) write(LOGUNIT,*) 'BFM uses XIOS output system via NEMO'
         ! Set var_ids values according to the XIOS file_def is done in trc_dia_bfm
      ELSE
         if (lwp) write(LOGUNIT,*) 'BFM output without XIOS to be implemented' 
         STOP
      ENDIF 

      ! Initialize water column optics
      !-------------------------------------------------------
      IF ( ChlAttenFlag == 2 ) THEN
         ierr = trc_opt_alloc()
         CALL trc_opt_ini
         parlux = p_PAR
         LightLocationFlag = 1 ! light is already integrated over the cell
      ENDIF

#ifdef INCLUDE_PELFE
      ! Iron supply from sediments (see Aumont & Bopp, 2005)
      !-------------------------------------------------------
      ALLOCATE (ironsed(jpi,jpj,jpk))
      ironsed = ZERO
      IF ( p_rN7fsed > 0. ) THEN
         CALL iom_open ( TRIM( "bottom_fraction.nc" ), numfld )
         CALL iom_get  ( numfld, 1, TRIM( "btmfrac" ), ironsed(:,:,:), 1 )
         CALL iom_close( numfld )
         ! iron release dependence on depth (metamodel of Middelburg et al.,1996)
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk )
              zexp  = MIN( 8.,( gdept(ji,jj,jk,Kmm) / 500. )**(-1.5) )
              zdexp = -0.9543 + 0.7662 * LOG( zexp ) - 0.235 * LOG( zexp )**2
              ironsed(ji,jj,jk) = ironsed(ji,jj,jk) * MIN( 1., EXP( zdexp ) / 0.5 )
         END_3D
         ! total supply (TODO this is killing execution .. need to check)
         ztraf = glob_sum('ironsed', ironsed * spread(e1e2t,3,jpk) * p_rN7fsed * 365. * 1.e-15 * tmask )
         LEVEL1 ''
         LEVEL1 'trc_ini_my_trc: Read fraction mask from bottom_fraction.nc (varname: btmfrac)'
         LEVEL1 '  Constant Iron flux from sediments: ', p_rN7fsed
         LEVEL1 '  Total iron load from Sediment [Gmol/y] : ', ztraf
      ENDIF
#endif

      ! Scaling factor for benthic return coefficients !TODO move this into BFM
      !-------------------------------------------------------
      ! Depth dependence from Middelburg et al. (1996) metamodel (see par. 3.4)
      ! if ( p_depscale > ZERO ) then
      !    allocate(rtmp1D(NO_BOXES))
      !    rtmp1D = pack(gdept_0,SEAmask)
      !    DO i = 1, NO_BOXES_XY
      !       j = BOTindices(i)
      !       zexp  = MIN( 8.,( rtmp1D(j) / p_depscale )**(-1.5) )
      !       zdexp = -0.9543 + 0.7662 * LOG( zexp ) - 0.235 * LOG( zexp )**2
      !       RETFAC(i) = MIN( 1., EXP( zdexp ) / 0.5 )
      !    END DO
      !    deallocate(rtmp1D)
      ! endif

      IF (bfm_lwp) THEN
            LEVEL1 ' '
            LEVEL1 '         BFM INITIALIZATION ... DONE!          '
            LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
            LEVEL1 ' '
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_sms_my_trc')
      !
157 FORMAT(a4, 1x, a10  , 1x, a25, 1x, a3, 1x, a10  , 1x, a10  , 1x, a10  , 1x, a10  , 3x, a3, 3x, a3, 3x, a3)
158 FORMAT(i4, 1x, E10.3, 1x, a25, 1x, a3, 1x, E10.3, 1x, E10.3, 1x, E10.3, 1x, E10.3, 3x, L3, 3x, L3, 3x, L3)
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

   !!======================================================================
END MODULE trcini_my_trc
