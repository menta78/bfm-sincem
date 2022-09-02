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
   USE global_mem, ONLY: LOGUNIT, bfm_lwp, RLEN, ZERO
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
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: tflux, zwgt
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

         ! Compute vertical sinking with upwind 1st scheme
         !-------------------------------------------------------
         CALL sink_upwind1(-sink_rates(:,:,:,jn), Kmm, Krhs, jl, tflux)
         !tflux = -sink_rates(:,:,:,jn)

         ! Compute vertical sinking with upwind 3rd scheme
         !-------------------------------------------------------
         !if (jn == 1) call upwind3_weights(zwgt, Kmm)
         !CALL sink_upwind3(-sink_rates(:,:,:,jn), Kmm, Krhs, jl, zwgt, tflux)

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


   SUBROUTINE sink_upwind3(ws_in, Kmm, Krhs, jl, zwgt, flux)
      !!----------------------------------------------------------------------
      !!                 ***  sink_upwind3  ***
      !!
      !! ** Purpose :  Third order upwind scheme for vertical sinking
      !!
      !! ** Method  : Sets additional vertical velocity field and computes
      !!            resulting advection using a conservative 3rd upwind
      !!            scheme with QUICKEST TVD limiter, based on GOTM module
      !             adv_center.F90 (www.gotm.net).
      !!            Currently assuming zero flux at sea surface and sea floor.
      !!----------------------------------------------------------------------
      ! NEMO
      USE trc, ONLY: tr
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: ws_in, zwgt
      INTEGER, INTENT(in) :: Kmm, Krhs  ! time level indices
      INTEGER, INTENT(in) :: jl         ! tracer index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(out) :: flux
      !
      INTEGER   :: ji, jj, jk
      REAL(wp)  :: cmax_no, z2dt
      LOGICAL   :: ll_euler
      REAL(wp), DIMENSION(jpk) :: tr_it, tr_u, tr_d, tr_c, tr_slope, c_no, &
                                 flux_lim, phi_lim, x_fac, w_if, flux_if
      !!----------------------------------------------------------------------

      ! Set time step
      !-------------------------------------------------------
      ll_euler = (l_1st_euler .OR. ln_top_euler)
      IF( ll_euler ) THEN
         z2dt = rdt
      ELSE
         z2dt = 2._wp * rdt
      ENDIF

      ! Initialize interface fluxes with zero
      !-------------------------------------------------------
      flux = 0.0_wp

      ! Interior vertical velocities and include them in source array.
      !-------------------------------------------------------
      DO jj = 1,jpj
         DO ji = 1,jpi
            IF ( mbkt(ji,jj)>1 ) THEN
               jk = mbkt(ji,jj)
               ! surface
               w_if(1) = 0._wp
               ! inner
               w_if(2:jk) = (zwgt(ji,jj,2:jk)*ws_in(ji,jj,1:jk-1)&
                  +(1._wp-zwgt(ji,jj,1:jk-1))*ws_in(ji,jj,2:jk))*DAY_PER_SEC
               ! sea floor and below
               w_if(jk+1:) = 0._wp
               !
               c_no(2:jk) = abs(w_if(2:jk))*z2dt/ &
                   ( 0.5_wp*(e3t(ji,jj,2:jk,Kmm) + &
                   e3t(ji,ji,1:jk-1,Kmm)) )
               cmax_no = MAXVAL(c_no(2:jk))

               ! number of iterations:
               IF ( sn_cfctl%l_prtctl .AND. (cmax_no .GT. 1) ) THEN
                  WRITE(numout,*) 'sink_upwind3:'
                  WRITE(numout,*) '   Maximum Courant number is ', cmax_no
               ENDIF

               tr_it(1:jk) = tr(ji,jj,1:jk,jl,Kmm)

               ! Compute slope ratio
               !-------------------------------------------------------
               IF (jk.gt.2) THEN !More than 2 vertical wet points
                  IF (jk.gt.3) THEN
                     ! upward movement
                     WHERE (w_if(3:jk-1) .GE. 0._wp)
                        tr_u(3:jk-1)=tr_it(4:jk)
                        tr_c(3:jk-1)=tr_it(3:jk-1)
                        tr_d(3:jk-1)=tr_it(2:jk-2)
                     ! downward movement
                     ELSEWHERE
                        tr_u(3:jk-1)=tr_it(1:jk-3)
                        tr_c(3:jk-1)=tr_it(2:jk-2)
                        tr_d(3:jk-1)=tr_it(3:jk-1)
                     ENDWHERE
                  ENDIF
                  IF (w_if(2) .GE. 0._wp) THEN
                     tr_u(2)=tr_it(3)
                     tr_c(2)=tr_it(2)
                     tr_d(2)=tr_it(1)
                  ELSE
                     tr_u(2)=tr_it(1)
                     tr_c(2)=tr_it(1)
                     tr_d(2)=tr_it(2)
                  ENDIF
                  IF (w_if(jk).ge.0._wp) THEN
                     tr_u(jk)=tr_it(jk)
                     tr_c(jk)=tr_it(jk)
                     tr_d(jk)=tr_it(jk-1)
                  ELSE
                     tr_u(jk)=tr_it(jk-2)
                     tr_c(jk)=tr_it(jk-1)
                     tr_d(jk)=tr_it(jk)
                  ENDIF
               ELSE !only 2 vertical wet points, i.e. only 1 interface
                  IF (w_if(jk).ge.0._wp) THEN
                     tr_u(2)=tr_it(2)
                     tr_c(2)=tr_it(2)
                     tr_d(2)=tr_it(1)
                  ELSE
                     tr_u(2)=tr_it(1)
                     tr_c(2)=tr_it(1)
                     tr_d(2)=tr_it(2)
                  ENDIF
               ENDIF
               WHERE (abs(tr_d(2:jk)-tr_c(2:jk)) .GT. 1.e-10_wp)
                  tr_slope(2:jk)= &
                     (tr_c(2:jk)-tr_u(2:jk))/ &
                     (tr_d(2:jk)-tr_c(2:jk))
               ELSEWHERE
                     tr_slope(2:jk)=SIGN(1._wp,w_if(2:jk))* &
                        (tr_c(2:jk)-tr_u(2:jk))*1.e10_wp
               ENDWHERE

               ! QUICKEST flux limiter
               !-------------------------------------------------------
               x_fac(2:jk)=(1._wp-2._wp*c_no(2:jk))/6._wp
               phi_lim(2:jk)=(0.5_wp+x_fac(2:jk)) + &
                  (0.5_wp-x_fac(2:jk))*tr_slope(2:jk)
               flux_lim(2:jk)=max( 0._wp, &
                  min( phi_lim(2:jk),2._wp/(1._wp-c_no(2:jk)), &
                  2._wp*tr_slope(2:jk)/(c_no(2:jk)+1.e-10_wp)) )

               ! Compute limited flux
               !-------------------------------------------------------
               flux_if(1) = 0._wp
               flux_if(2:jk) = w_if(2:jk)* &
                  ( tr_c(2:jk) + &
                     0.5_wp*flux_lim(2:jk)*(1._wp-c_no(2:jk))* &
                     (tr_d(2:jk)-tr_c(2:jk)) )

               tr(ji,jj,1:jk-1,jl,Krhs) = tr(ji,jj,1:jk-1,jl,Krhs) &
                  + flux_if(2:jk)/e3t(ji,jj,1:jk-1,Kmm)
               tr(ji,jj,2:jk,jl,Krhs) = tr(ji,jj,2:jk,jl,Krhs) &
                  - flux_if(2:jk)/e3t(ji,jj,2:jk,Kmm)
               ! interpolate flux to cell centre:
               flux(ji,jj,1:jk-1) = 0.5 * (flux_if(1:jk-1) + flux_if(2:jk))

            END IF ! Level check
         END DO
      END DO

     
      END SUBROUTINE sink_upwind3


   SUBROUTINE upwind3_weights(zwgt, Kmm)
      !!----------------------------------------------------------------------
      !!                 ***    ***
      !!
      !! ** Purpose : Cell centered weights for third order upwind scheme
      !!
      !! ** Method  : Linearly interpolate to velocities at the interfaces 
      !!            between vertical layers. Note:
      !!            - interface k sits between cell centre k and k-1,
      !!            - k [1,jpk] increases downwards
      !!            - upward velocity is positive, downward velocity is negative
      !!            - surface and bottom velocities =0 (no friction)
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) :: zwgt
      INTEGER, INTENT(in) :: Kmm  ! time level indices
      !
      INTEGER :: ji, jj, jk
      !!----------------------------------------------------------------------
      zwgt(:,:,:) = 0._wp
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Compute only on (ji,jj) gridcells with layers > 1
            !-------------------------------------------------------
            IF ( mbkt(ji,jj)>1 ) THEN
               jk = mbkt(ji,jj)
               ! surface
               zwgt(ji,jj,1) = 0._wp
               ! inner
               zwgt(ji,jj,2:jk) = &
                  e3t(ji,jj,2:jk,Kmm) / ( e3t(ji,jj,1:jk-1,Kmm) + e3t(ji,jj,2:jk,Kmm) )
               ! sea floor and below
               zwgt(ji,jj,jk+1:) = 0._wp
            ELSE
               zwgt(ji,jj,:) = 0._wp
            ENDIF
         END DO
      END DO

   END SUBROUTINE upwind3_weights


END MODULE trcset
