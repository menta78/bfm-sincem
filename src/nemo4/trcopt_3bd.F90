SUBROUTINE trc_opt_3bd(Kmm, ezchl)
   !!---------------------------------------------------------------------
   !!                    ***  ROUTINE trcopt_3bd  ***
   !!
   !! ** Purpose : Compute the light availability in the water column
   !!              depending on depth and chlorophyll concentration
   !!
   !! ** Method  : 3-band (Morel)
   !!---------------------------------------------------------------------
#include "do_loop_substitute.h90"
#include "domzgr_substitute.h90"
#include "cppdefs.h"
   !
   USE oce_trc
   USE trc,        ONLY: etot, tr
   USE trcopt,     ONLY: zeps
   !
   USE global_mem, ONLY: LOGUNIT, bfm_lwp, RLEN, ZERO, ONE
   USE mem_PAR,    ONLY: p_PAR, p_eps0, p_epsESS, p_epsR6, xepsRGB, p_epsIR
   USE mem,        ONLY: ppR6c, iiC, iiL, iiPhytoPlankton, ppPhytoPlankton
   USE mem_Param,  ONLY: p_small
   USE mem_Phyto,  ONLY: p_qlcPPY, p_epsChla
   USE constants,  ONLY: E2W
   !
   IMPLICIT NONE

   INTEGER, INTENT(in) :: Kmm
   REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: ezchl
   INTEGER :: ji, jj, jk, jl, irgb
   REAL(wp), DIMENSION(jpi,jpj,jpk) :: zchl, ze1, ze2, ze3
   REAL(wp), DIMENSION(jpi,jpj,jpk,3) :: eps_bgr ! 1: Blue, 2: Green, 3: Red
   REAL(wp) :: ztmp, par_rgb
   ! 
   par_rgb = p_PAR / 3._RLEN
   !
   ! Abiotic extinction factors (not in cpl v3.6)
   !-------------------------------------------------------
   !zeps = p_eps0 + p_epsR6 * tr(:,:,:,ppR6c,Kmm)
   
   ! RGB attenuation coefficients
   !-------------------------------------------------------
   zchl = ZERO
   DO jl = 1 , iiPhytoPlankton
      IF ( ppPhytoPlankton(jl, iiL) == 0) THEN
         zchl = zchl + tr(:,:,:,ppPhytoPlankton(jl, iiC), Kmm) * p_qlcPPY(jl)
      ELSE
         zchl = zchl + tr(:,:,:,ppPhytoPlankton(jl, iiL), Kmm)
      ENDIF
   ENDDO

   DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1 )
      ztmp = MIN(  10. , MAX( 0.05, zchl(ji,jj,jk) )  )
      irgb = NINT( 41 + 20.* LOG10( ztmp ) + p_small )
      !LEVEL1 'opt', ji, jj, jk, irgb, zchl(ji,jj,jk)
      !
      eps_bgr(ji,jj,jk,1) = xepsRGB(1,irgb)
      eps_bgr(ji,jj,jk,2) = xepsRGB(2,irgb)
      eps_bgr(ji,jj,jk,3) = xepsRGB(3,irgb)
   END_3D
   
   ! Photosynthetically active radiation (Watt/m2)
   !-------------------------------------------------------
   ze1(:,:,1) = par_rgb * (qsr(:,:) + p_small)
   ze2(:,:,1) = ze1(:,:,1)
   ze3(:,:,1) = ze1(:,:,1)
   
   DO jk = 1, jpk-1
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         ze1(ji,jj,jk+1) = ze1(ji,jj,jk) * EXP( -eps_bgr(ji,jj,jk,1) * e3w(ji,jj,jk,Kmm) )
         ze2(ji,jj,jk+1) = ze2(ji,jj,jk) * EXP( -eps_bgr(ji,jj,jk,2) * e3w(ji,jj,jk,Kmm) )
         ze3(ji,jj,jk+1) = ze3(ji,jj,jk) * EXP( -eps_bgr(ji,jj,jk,3) * e3w(ji,jj,jk,Kmm) )
      END_2D
   ENDDO
   etot = ze1 + ze2 + ze3

   ! weighted broadband diffuse attenuation coefficient
   !-------------------------------------------------------
   zeps = (ze1 * eps_bgr(:,:,:,1) + ze2 * eps_bgr(:,:,:,2) + ze3*eps_bgr(:,:,:,3) ) / (etot + p_small )

   IF (ln_qsr_bio) THEN
      etot3(:,:,1) = (ONE-p_PAR) * qsr(:,:)
      ze1(:,:,1) = par_rgb * qsr(:,:)
      ze2(:,:,2) = ze1(:,:,1)
      ze3(:,:,3) = ze1(:,:,1)
      DO jk = 1, jpk-1
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            ! Infrared
            etot3(ji,jj,jk+1) = etot3(ji,jj,jk) * EXP( -p_epsIR * e3t(ji,jj,jk,Kmm) )
            ! BGR
            ze1(ji,jj,jk+1) = ze1(ji,jj,jk) * EXP( -eps_bgr(ji,jj,jk,1) * e3t(ji,jj,jk,Kmm) )
            ze2(ji,jj,jk+1) = ze2(ji,jj,jk) * EXP( -eps_bgr(ji,jj,jk,2) * e3t(ji,jj,jk,Kmm) )
            ze3(ji,jj,jk+1) = ze3(ji,jj,jk) * EXP( -eps_bgr(ji,jj,jk,3) * e3t(ji,jj,jk,Kmm) )
         END_2D
      ENDDO
      etot3 = etot3 + ze1 + ze2 + ze3
   ENDIF

END SUBROUTINE trc_opt_3bd
