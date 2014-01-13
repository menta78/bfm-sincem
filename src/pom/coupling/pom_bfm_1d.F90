!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pom_bfm.f
!
! !INTERFACE:
       subroutine pom_bfm_1d
!
! !DESCRIPTION:
!  BFM coupling with POM
!  Time-marching routines: trend computations (bio+transport)
!  and integration
!
!
! !USES:

       use api_bfm
       use Service
       use constants,  only:SEC_PER_DAY
       use POM
       use Mem
       IMPLICIT NONE

! !INCLUDES:
!       include time.h
!
! !INPUT PARAMETERS:
   
!
! !REVISION HISTORY:
!  Original author(s): M. Butenschoen,  M. Vichi, 
!                      M. Zavatarelli, L. Polimene 
!
! !LOCAL VARIABLES:
!----------------
       real(RLEN) :: TT
       save TT
       logical ::  first
       save first
       data first /.true./
       integer :: ntime,k
!----------------


!EOP
!-----------------------------------------------------------------------
!BOC
!-------------------------------------------------------
! Pass physical variables into bfm
! compute extinction coefficients
! compute vertical light distribution
!-------------------------------------------------------
       call env_forcing_pom_bfm_1d
!       call flush(6)
!-------------------------------------------------------
! calculate biological processes 
!-------------------------------------------------------
       call EcologyDynamics
!       call flush(6)
!-------------------------------------------------------
!VERTICAL DIFFUSION   and Integration of BFM state Vars    
! SOURCE SPLITTING 
!-------------------------------------------------------
!      do k = 1, kb-1
!     write(6,*) 'PO4 BEF VERTDIFF', D3state (2,k), N1p(k),D3SOURCE(2,k)
!     enddo
        call vertdiff_bio_1d 
!       call flush(6)
!-------------------------------------------------------
!  leap frog integration of 2d state var's
!-------------------------------------------------------
#ifdef INCLUDE_BEN
       call lf2d
!        call flush(6)
#endif
!-------------------------------------------------------
! update time
!-------------------------------------------------------
        ntime=iint
!-------------------------------------------------------
! Write output
!-------------------------------------------------------
       if(first) then
           TT=time-time0-(dti/SEC_PER_DAY)
           out_delta=savef
!           TT = 0.0
         first=.false.
       endif
          TT = TT + dti/SEC_PER_DAY
!          write(6,*) 'TIME IN DIA', time, time0, TT
!          call flush(6)
       call pom_dia_bfm(time,ntime,TT)
!         call flush(6)
!-------------------------------------------------------
! Reset trend arrays
!-------------------------------------------------------
       call ResetFluxes
!         call flush(6)
       return

       end subroutine pom_bfm_1d

!EOC
!-----------------------------------------------------------------------


