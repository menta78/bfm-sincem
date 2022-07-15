MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE trcsms_cfc       ! CFCs & SF6 
   USE trcsms_age       ! AGE tracer
   USE trcwri_cfc
   USE trcwri_age
   USE api_bfm,    ONLY: bio_calc
   USE global_mem, ONLY: LOGUNIT
   USE iom
   USE in_out_manager
   USE trcsub
   USE trcrst
   USE trcwri
   USE trcbc,      ONLY: trc_bc_read
   USE trcdiabfm,  ONLY: trc_dia_bfm

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step
   
   !!----------------------------------------------------------------------
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !! Additional software compliant with GPL
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_stp( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      CHARACTER (len=25)    ::  charout
      !!-------------------------------------------------------------------
#ifdef DEBUG
      write(LOGUNIT,*) 'Timestep: ',kt, nit000, nn_dttrc
#endif

      IF( nn_timing == 1 )   CALL timing_start('trc_stp')

      !---------------------------------------------
      ! Check the main BFM flag
      !---------------------------------------------
      IF (bio_calc) THEN
          !
          ! Averaging physical variables for sub-stepping
          IF ( nn_dttrc /= 1 )     CALL trc_sub_stp( kt )
          !
          !---------------------------------------------
          ! Proceed only every nn_dttrc
          !---------------------------------------------   
          IF ( MOD( kt , nn_dttrc ) == 0 ) THEN  
             ! 
                                   CALL trc_bc_read  ( kt )   ! read/update Boundary Conditions
                                   CALL trc_bfm      ( kt )   ! main call to BFM
             ! TOP restart
             IF( lk_age .or. lk_cfc )  &
                         &         CALL trc_rst_opn  ( kt )
             IF( lrst_trc )        CALL trc_rst_cal  ( kt, 'WRITE' )
             ! TOP sms
             IF( lk_cfc        )   CALL trc_sms_cfc  ( kt )   ! surface fluxes of CFC
             IF( lk_age        )   CALL trc_sms_age  ( kt )   ! AGE tracer
             ! 
                                   CALL trc_trp_bfm  ( kt )   ! transport of BFM & TOP tracers
                                   CALL trc_dia_bfm  ( kt )   ! diagnostic output for BFM
             !
             ! TOP output & restart
             IF( lk_iomput     ) THEN 
                IF( lk_age     )   CALL trc_wri_age
                IF( lk_cfc     )   CALL trc_wri_cfc
             ENDIF
             IF( lrst_trc      )   CALL trc_rst_wri  ( kt )
             !
             ! reset physical variables after sub-stepping
             IF( nn_dttrc /= 1 )   CALL trc_sub_reset( kt )

          ENDIF 

      END IF

      IF( nn_timing == 1 )   CALL timing_stop('trc_stp')
      FLUSH(LOGUNIT)

   END SUBROUTINE trc_stp

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif

   !!======================================================================
END MODULE trcstp
