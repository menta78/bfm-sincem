#include "INCLUDE.h"
 SUBROUTINE LF2D
!
#ifdef POM_BFM
!
#ifdef INCLUDE_BEN
 use global_mem, ONLY:RLEN,ZERO
! use Mem, ONLY:D2STATE,D2SOURCE,NO_BOXES_XY,NO_D2_BOX_STATES
 use Mem, ONLY:D2STATE_BEN,NO_BOXES_XY,NO_D2_BOX_STATES_BEN,D2SOURCE_BEN

!#ifndef ONESOURCE
#ifdef EXPLICIT_SINK
 use Mem, Only: D2SINK_BEN
#endif

 use api_bfm, ONLY:D2STATEB_BEN
 use Service, ONLY: smothass,deltat
 implicit none
 integer m, n, i
 real(RLEN) :: dtie2
 real(RLEN) :: ccc_tmp2D(NO_D2_BOX_STATES_BEN,NO_BOXES_XY)
! real(RLEN) :: ccc_tmp2D(NO_D2_BOX_STATES_BEN,NO_BOXES_BEN)
 real(RLEN) :: p_small=1.0D-80
!
  ccc_tmp2D=ZERO
!
 dtie2 = deltat*2.0
   PRINT*,'ccc_tmp2D',size(ccc_tmp2D)
!
#ifndef D1SOURCE
!        
   do n = 1, NO_D2_BOX_STATES_BEN
!#ifdef ONESOURCE
#ifndef EXPLICIT_SINK
!          ccc_tmp2D(:,:)=ccc_tmp2D(:,:)+D2SOURCE_BEN(n,:)
#else
          ccc_tmp2D(:,:)=ccc_tmp2D(:,:)+(D2SOURCE_BEN(:,n,:)-D2SINK_BEN(:,n,:))
#endif
   end do 
#endif
!
#ifdef D1SOURCE
          ccc_tmp2D=D2STATEB_BEN+(D2SOURCE_BEN*dtie2)
#else
          ccc_tmp2D=D2STATEB_BEN+(ccc_tmp2D*dtie2)
#endif
!
!
!       -----MIX THE TIME STEP-----
!
!call ftrace_region_begin("LF2D3")
!
          do  n = 1,NO_D2_BOX_STATES_BEN
          D2STATE_BEN(n,:)=D2STATE_BEN(n,:)+.5*smothass*(ccc_tmp2D(n,:)+ &
                       D2STATEB_BEN(n,:)-2.*D2STATE_BEN(n,:))
          end do
!
!call ftrace_region_end("LF2D3")
!
!                =====================================================
!                RESTORE TIME SEQUENCE
!                =====================================================
!         D2STATEB=D2STATE
!         D2STATE=ccc_tmp2D
!
!call ftrace_region_begin("LF2D4")
!
          do  n = 1,NO_D2_BOX_STATES_BEN
          ccc_tmp2D(n,:)=max(p_small,ccc_tmp2D(n,:))
          end do
!
          do  n = 1,NO_D2_BOX_STATES_BEN
          D2STATEB_BEN(n,:)=D2STATE_BEN(n,:)
          D2STATE_BEN(n,:)=ccc_tmp2D(n,:)
          end do
!
!
!call ftrace_region_end("LF2D4")
!
!
#endif
!
#endif
      return
      end
