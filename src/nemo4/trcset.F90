MODULE trcset
   !!======================================================================
   !!                         ***  MODULE trcset  ***
   !! TOP :   Vertical sinking of tracers
   !!======================================================================
#include "domzgr_substitute.h90"
#include "do_loop_substitute.h90"
#include"cppdefs.h"
   ! NEMO
   USE oce_trc
   USE trc              ! ocean passive tracers variables
   ! BFM
   USE global_mem, ONLY: LOGUNIT, bfm_lwp, RLEN, ZERO, SkipBFMCore
   use constants,  ONLY: DAY_PER_SEC, SEC_PER_DAY
   USE time,       ONLY: bfmtime

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_set

CONTAINS

   SUBROUTINE trc_set(kt, Kmm, Krhs)
      !!----------------------------------------------------------------------
      !!                 ***  trc_set  ***
      !!
      !! ** Purpose :  Apply vertical sinking & settling to tracers
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      ! NEMO
      USE par_my_trc
      USE iom,        ONLY: iom_put
      ! BFM
      USE api_bfm,    ONLY: var_names
      USE mem_PelSinkSet
      !
      INTEGER, INTENT(in) :: kt         ! ocean time-step index
      INTEGER, INTENT(in) :: Kmm, Krhs  ! time level indices
      INTEGER             :: ji, jj, jk, jn, jl
      REAL(wp)            :: zfact
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: tflux
      !!----------------------------------------------------------------------

      IF( ln_timing )   CALL timing_start('trc_set')
      !
      IF (kt == nittrc000) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_set: Vertical Sinking'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF
      
      ! Loop over tracers
      !-------------------------------------------------------
      DO jn = 1 , sink_vars

         jl = sink_var_map(jn)

         !LEVEL1 'sink of ', var_names(jl), SINKD3STATE(jl)%dosink

         ! Compute vertical sinking with upwind scheme
         !-------------------------------------------------------
         CALL sink_upwind1(-sink_rates(:,:,:,jn), Kmm, Krhs, jl, tflux)
         !tflux = -sink_rates(:,:,:,jn)

         ! Output of fluxes as mmol / m2 / day
         !-------------------------------------------------------
         if (bfm_iomput) &
            call iom_put ("sink"//TRIM(var_names(jl)) , tflux * SEC_PER_DAY)

      ENDDO

      IF( ln_timing )   CALL timing_stop('trc_set')

   END SUBROUTINE trc_set


   SUBROUTINE sink_upwind1(ws_in, Kmm, Krhs, jl, flux)
      !!----------------------------------------------------------------------
      !!                 ***  sink_upwind1  ***
      !!
      !! ** Purpose :  First order upwind scheme for vertical sinking
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      ! NEMO
      USE trc, ONLY: tr
      !  
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: ws_in
      INTEGER, INTENT(in) :: Kmm, Krhs  ! time level indices
      INTEGER, INTENT(in) :: jl         ! tracer index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(out) :: flux
      !
      INTEGER   :: ji, jj, jk
      REAL(wp)  :: Yd, Yc, Yu
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: cu, ws
      !!----------------------------------------------------------------------

      ! Initialize interface fluxes with zero
      !-------------------------------------------------------
      cu   = 0.0_wp
      flux = 0.0_wp
      ws = 0.0_wp
 
      ! Downward shifted velocities array
      !-------------------------------------------------------
      DO jk = 1, jpk-1
         ws(:,:,jk+1) = ws_in(:,:,jk) * DAY_PER_SEC *tmask(:,:,jk+1)
      ENDDO

      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 2, jpkm1 )

         ! Check the velocity direction
         !-------------------------------------------------------
         if (ws(ji,jj,jk) .gt. 0.0_wp) then
            ! positive speed
            Yu = tr(ji,jj,jk+1,jl,Kmm)          ! upstream value
            Yc = tr(ji,jj,jk,jl,Kmm)            ! central value
            if (jk .gt. 1) then
               Yd = tr(ji,jj,jk-1,jl,Kmm)       ! downstream value
            else
               Yd = tr(ji,jj,1,jl,Kmm)
            end if
         else
            ! negative speed
            if (jk .gt. 2) then
               Yu = tr(ji,jj,jk-2,jl,Kmm)       ! upstream value
            else
               Yu = tr(ji,jj,1,jl,Kmm)          ! upstream value
            end if
            Yc = tr(ji,jj,jk-1,jl,Kmm)          ! central value
            Yd = tr(ji,jj,jk,jl,Kmm)            ! downstream value
         end if

         ! compute the mass flow across the central interface
         !-------------------------------------------------------
         cu(ji,jj,jk) = ws(ji,jj,jk) * Yc

      END_3D
   
      ! Force to zero surface and bottom
      !-------------------------------------------------------
      cu(:,:,1) = 0.0_wp
      cu(:,:,jpk) = 0.0_wp

      DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1 )

         ! Add vertical advection trend to tracer RHS
         !-------------------------------------------------------
         tr(ji,jj,jk, jl, Krhs) = tr(ji,jj,jk,jl, Krhs) + (cu(ji,jj,jk+1) - cu(ji,jj,jk)) / e3t(ji,jj,jk,Kmm)

         ! Flux at the cell center (for diagnostics)
         !-------------------------------------------------------
         flux(ji,jj,jk) = 0.5 * ( cu(ji,jj,jk+1) + cu(ji,jj,jk) ) * tmask(ji,jj,jk)

      END_3D


   END SUBROUTINE sink_upwind1

END MODULE trcset
