SUBROUTINE trc_nam_bfm()
   !!======================================================================
   !!                ***  SUBROUTINE trc_nam_bfm  ***
   !! BFM :   Initialize tracers and read namelists  (overwrite namtrc)
   !!======================================================================
#include "cppdefs.h"
   ! NEMO
   USE par_trc, ONLY: jp_bgc, jp_bgc_b, jpk_b
   USE trcnam,  ONLY: sn_tracer
   USE in_out_manager, ONLY: lwp, numout, nitend
   USE dom_oce, ONLY: rn_Dt, narea, nyear, nmonth, nday, tmask
   USE par_oce, ONLY: jpi, jpj, jpk
   USE trc,     ONLY: tr, ln_rsttr, ln_trcdta, ln_trcbc, ln_top_euler, nittrc000, cn_trcrst_in
   USE par_my_trc, ONLY: var_map, jp_bgc_b, bottom_level, bfm_iomput
   ! BFM
   USE constants,  ONLY: SEC_PER_DAY
   USE global_mem, ONLY: RLEN, ZERO, LOGUNIT, bfm_lwp, ALLTRANSPORT
   USE api_bfm,    ONLY: parallel_rank, bio_setup, SEAmask, init_bfm, stPelStateS, stPelStateE, &
                         save_delta, time_delta, out_delta, update_save_delta, &
                         InitVar, var_names, var_long, var_units, SRFindices, BOTindices, &
                         bfm_init, in_rst_fname
   USE time,       ONLY: bfmtime, julian_day, calendar_date
   USE mem,        ONLY: D3STATETYPE, NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z, &
                         NO_BOXES, NO_BOXES_XY, NO_STATES, NO_D3_BOX_STATES, & 
                         NO_BOXES_Z_BEN, NO_BOXES_BEN, NO_STATES_BEN, NO_D2_BOX_STATES_BEN
#ifdef INCLUDE_SEAICE
   USE mem,        ONLY: NO_BOXES_Z_ICE, NO_BOXES_ICE, NO_STATES_ICE
#endif
   
#ifdef CCSMCOUPLED
   USE nemogcm, ONLY: logfile
#endif
   !
   IMPLICIT NONE

   INTEGER     :: yy, mm, dd, hh, nn, jn
   REAL(RLEN)  :: julianday
   LOGICAL     ::  lltrcbc
   !!----------------------------------------------------------------------

   IF(lwp) WRITE(numout,*)
   IF(lwp) WRITE(numout,*) 'trc_nam_bfm : read BFM tracer namelists'
   IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

   ! Assign the rank of the process (meaningful only with key_mpp)
   !-------------------------------------------------------
   parallel_rank = narea - 1

   ! Initialize bfm time
   !-------------------------------------------------------
   WRITE(bfmtime%datestring,'(I4.4,a,I2.2,a,I2.2)') nyear,'-',nmonth,'-',nday
   WRITE(bfmtime%date0,'(I4.4,I2.2,I2.2)') nyear,nmonth,nday
   CALL julian_day(nyear,nmonth,nday,0,0,julianday)
   bfmtime%time0    = julianday
   bfmtime%timeEnd  = julianday + ( ( REAL(nitend - nittrc000, RLEN) ) * rn_Dt ) / SEC_PER_DAY
   bfmtime%step0    = nittrc000 - 1
   bfmtime%timestep = rn_Dt
   bfmtime%stepnow  = nittrc000 - 1
   bfmtime%stepEnd  = nitend
   CALL calendar_date(bfmtime%timeEnd,yy,mm,dd,hh,nn)
   WRITE(bfmtime%dateEnd,'(I4.4,I2.2,I2.2)') yy,mm,dd

   ! Force Euler timestepping for TOP
   !-------------------------------------------------------
   ln_top_euler = .TRUE.
   IF(lwp) WRITE(numout,*) '              force euler integration (ln_top_euler=T)'

   ! Set ocean mask and bottom level
   !-------------------------------------------------------
   ALLOCATE(SEAmask(jpi,jpj,jpk))
   SEAmask = .FALSE.
   !if (lk_c1d) then !TODO check how 1d works
   !   where (tmask(2,2,:) > ZERO)
   !     SEAmask(2,2,:) = .TRUE.
   !   elsewhere
   !     SEAmask(2,2,:) = .FALSE.
   !   end where
   !else
      WHERE (tmask > ZERO)
        SEAmask = .TRUE.
      ELSEWHERE
        SEAmask = .FALSE.
      ENDWHERE
   !end if
   !
   ALLOCATE(bottom_level(jpi,jpj))
   bottom_level = SUM(tmask, DIM=3)
   ALLOCATE(SRFindices(1), BOTindices(1))

   ! Set the dimensions
   !-------------------------------------------------------
   NO_BOXES_X  = 1 ! jpi
   NO_BOXES_Y  = 1 ! jpj
   NO_BOXES_Z  = jpk
   NO_BOXES    = jpk
   NO_BOXES_XY = 1
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES

   NO_BOXES_Z_BEN  = 1
   NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
   NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN
#ifdef INCLUDE_SEAICE
   NO_BOXES_Z_ICE  = 1
   NO_BOXES_ICE = NO_BOXES_XY * NO_BOXES_Z_ICE
   NO_STATES_ICE = NO_BOXES_ICE * NO_D2_BOX_STATES_ICE
#endif

   ! Initialise memory arrays and log
   !-------------------------------------------------------
#ifdef CCSMCOUPLED
   CALL init_bfm( TRIM( logfile(4:len(logfile)) ) )
#else
   CALL init_bfm
#endif

   ! Force BFM restart with NEMO setting
   !-------------------------------------------------------
   IF ( ln_rsttr ) THEN
      bfm_init = 1
      in_rst_fname = TRIM(cn_trcrst_in)
      LEVEL1 'Restart BFM from input file '//TRIM(in_rst_fname)
   ELSE
      bfm_init = 0
      LEVEL1 'BFM start from initial conditions '
   ENDIF
   LEVEL1 ''

   ! Initialise state variable names and diagnostics
   !-------------------------------------------------------
   CALL set_var_info_bfm

   ! Allocate memory and read initial & boundary conditions setting
   !-------------------------------------------------------
   CALL init_var_bfm(bio_setup)

   ! Set output stepping
   !-------------------------------------------------------
   save_delta = bfmtime%step0
   IF ( .NOT. bfm_iomput ) call update_save_delta(out_delta,save_delta,time_delta)

   ! Setup NEMO namtrc settings
   !-------------------------------------------------------
   ALLOCATE(var_map(NO_D3_BOX_STATES)) ! mapping transported vars between nemo and bfm
   var_map = 0
   jp_bgc = 0
   DO jn = stPelStateS, stPelStateE
       IF (D3STATETYPE(jn-stPelStateS+1) == ALLTRANSPORT) THEN
           jp_bgc = jp_bgc + 1
           var_map(jn) =  jp_bgc
           sn_tracer(jp_bgc)%clsname = trim(var_names(jn))
           sn_tracer(jp_bgc)%cllname = trim(var_long(jn))
           sn_tracer(jp_bgc)%clunit = trim(var_units(jn))
           sn_tracer(jp_bgc)%llinit = InitVar(jn)%init == 2
           sn_tracer(jp_bgc)%llsbc = InitVar(jn)%sbc
           sn_tracer(jp_bgc)%llcbc = InitVar(jn)%cbc
           sn_tracer(jp_bgc)%llobc = InitVar(jn)%obc
       ENDIF
   ENDDO
   ln_trcdta = COUNT(sn_tracer(:)%llinit) > 0
   lltrcbc = ( COUNT(sn_tracer(:)%llsbc) + COUNT(sn_tracer(:)%llobc) + COUNT(sn_tracer(:)%llcbc) ) > 0
   IF ( ln_trcbc .AND. .NOT.lltrcbc) ln_trcbc = .FALSE.
   
   !write (LOGUNIT,*) 'map ', var_map
   ! benthic
   jp_bgc_b = NO_D2_BOX_STATES_BEN
   jpk_b = NO_BOXES_Z_BEN
#ifdef INCLUDE_SEAICE
   ! seaice
   jpk_i = NO_BOXES_Z_ICE
#endif

END SUBROUTINE trc_nam_bfm

   



