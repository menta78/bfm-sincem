#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: LF2D
!
!DESCRIPTION    
!
!  This routine calculates and solves for the leap-frog time step 
!  applied to 2-dimensional fields
!
!                                          Marco.Zavatarelli@unibo.it
! !INTERFACE
  SUBROUTINE LF2D
!
! !USES:
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
 use mem_Param, ONLY: p_small
!
 use global_mem, ONLY:RLEN,ZERO
!
#ifdef INCLUDE_BEN
 use Pom, ONLY: smoth,dti
!
 use mem , ONLY:D2STATE_BEN,NO_BOXES_XY,NO_D2_BOX_STATES_BEN,D2SOURCE_BEN, &
               NO_BOXES_Z_BEN, NO_STATES_BEN, NO_BOXES_BEN
!
 use api_bfm, ONLY:D2STATEB_BEN
!
#ifdef EXPLICIT_SINK
 use Mem, Only: D2SINK_BEN
#endif
!
!-------------------------------------------------------------------------!
!BOC
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
  IMPLICIT NONE
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Local Variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
! -----COUNTERS-----
!
  integer m, n, i
!
! -----TWICE THE TIME STEP-----
!
  real(RLEN) :: dti2
!
! -----TEMPORARY STORAGE-----
!
  real(RLEN) :: ccc_tmp2D(NO_D2_BOX_STATES_BEN,NO_BOXES_BEN)
!
! -----ZEROING-----
!
  ccc_tmp2D(:,:)=ZERO
!
! -----TWICE THE TIME STEP-----
!
  dti2 = dti*2.0_RLEN

!
! -----START LOOP OVER BENTHIC STATE VAR'S-----
! -----TO COMPOSE SOURCE TERM WHEN D1SOURCE IS NOT DEFINED -----
!
#ifndef D1SOURCE
  
   do n = 1, NO_D2_BOX_STATES_BEN

#ifndef EXPLICIT_SINK
       ccc_tmp2D(:,:)=ccc_tmp2D+D2SOURCE_BEN
#else
        ccc_tmp2D(:,:)=ccc_tmp2D(:,:)+(D2SOURCE_BEN(:,n,:)-D2SINK_BEN(:,n,:))
#endif
   end do 
#endif
!
!-----LEAP FROG INTEGRATION-----
!
#ifdef D1SOURCE
!
          ccc_tmp2D=D2STATEB_BEN+(D2SOURCE_BEN*dti2)
!
#else
!
          ccc_tmp2D=D2STATEB_BEN+(ccc_tmp2D*dti2)
!
#endif
!
!         
! -----CLIPPING (IF NEEDED....)-----
!
      do  n = 1,NO_D2_BOX_STATES_BEN
          ccc_tmp2D(n,:)=max(p_small,ccc_tmp2D(n,:))
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   Mix the time step
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!
     D2STATE_BEN(n,:)=D2STATE_BEN(n,:)+0.5_RLEN*smoth* &
                      (ccc_tmp2D(n,:)+D2STATEB_BEN(n,:)- &
                      2.0_RLEN*D2STATE_BEN(n,:))
!
! -----RESTORE TIME SEQUENCE-----
!
     D2STATEB_BEN(n,:)=D2STATE_BEN(n,:)
     D2STATE_BEN(n,:)=ccc_tmp2D(n,:)
     end do
!
#endif
      return
      end subroutine LF2D

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

