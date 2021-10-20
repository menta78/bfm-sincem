SUBROUTINE trc_opt_2bd(kt, Kmm)
   !!---------------------------------------------------------------------
   !!                    ***  ROUTINE trcopt_2bd  ***
   !!
   !! ** Purpose : Compute the light availability in the water column
   !!              depending on depth and chlorophyll concentration
   !!
   !! ** Method  : 2-band
   !!---------------------------------------------------------------------
#include "do_loop_substitute.h90"
#include "domzgr_substitute.h90"
   !
   USE oce_trc
   USE trc,        ONLY: etot, tr
   USE trcopt,     ONLY: zeps
   !
   USE global_mem, ONLY: LOGUNIT, bfm_lwp, RLEN, ZERO, ONE
   USE mem_PAR,    ONLY: p_PAR, p_eps0, p_epsESS, p_epsR6
   USE mem,        ONLY: ppR6c
   USE mem_Phyto,  ONLY: p_qlcPPY, p_epsChla
   USE constants,  ONLY: E2W
   !
   INTEGER, INTENT(in) :: kt
   INTEGER :: ji, jj, jk, jl
   REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpar, ztotv
   REAL(wp) :: r_e2w = 1._RLEN / E2W
   
   !
   ! Abiotic extinction factors
   !-------------------------------------------------------
   zeps = p_eps0 + p_epsR6 * tr(:,:,:,ppR6c,Kmm)
   
   ! Broadband linear attenuation  
   !-------------------------------------------------------
   DO jl = 1 , iiPhytoPlankton
      IF ( ppPhytoPlankton(jl, iiL) == 0) THEN
         zeps = zeps + p_epsChla(jl) * tr(:,:,:,ppPhytoPlankton(jl, iiC), Kmm) * p_qlcPPY(jl)
      ELSE
         zeps = zeps + p_epsChla(jl) * tr(:,:,:,ppPhytoPlankton(jl, iiL), Kmm)
      ENDIF
   ENDDO 

   ! Photosynthetically active radiation (Watt/m2)
   !-------------------------------------------------------
   etot(:,:,1) = p_PAR * (qsr(:,:)+p_small)
   
   DO jk = 1, jpk-1
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         etot(ji,jj,jk+1) = etot(ji,jj,jk) * EXP( -zeps(ji,jj,jk) * e3w(ji,jj,jk,Kmm)  )
      END_2D
   ENDDO

   IF (ln_qsr_bio) THEN
      etot3(:,:,1) = (ONE-p_PAR) * qsr(:,:)
      ztotv(:,:,1) = p_PAR * qsr(:,:)
      DO jk = 1, jpk-1
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            ! infrared
            etot3(ji,jj,jk+1) = etot3(ji,jj,jk) * EXP( -p_epsIRi * e3t(ji,jj,jk,Kmm)  )
            ! visible
            ztotv(ji,jj,jk+1) = ztotv(ji,jj,jk) * EXP( -ztotv(ji,jj,jk) * e3t(ji,jj,jk,Kmm)  )
         END_2D
      ENDDO
      etot3 = etot3 + zepsv
   ENDIF

END SUBROUTINE trc_opt_2bd
