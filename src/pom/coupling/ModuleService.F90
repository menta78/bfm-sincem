!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: Service
!
!DESCRIPTION    
! List of Fortran parameters
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!ONE-DIMENSIONAL BFM-POM SYSTEM
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !INTERFACE
  MODULE  Service
!
!
!................................................Marco.Zavatarelli@unibo.it
!................................................G.Mussap@sincem.unibo.it
!                                                    2014
! -----DEFINITION AND ALLOCATION OF PARAMETERS, SCALAR AND ARRAYS USED---
! -----FOR THE BFM-POM COUPLING-----
! 
! !USES:
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
  use POM,ONLY: KB,ilong
  use global_mem, ONLY:RLEN
!
!-------------------------------------------------------------------------!
!BOC
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Implicit typing is never allowed
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
  IMPLICIT NONE
! 
!     -----SURFACE NUTRIENTS-----
!
      real(RLEN)                        :: PO4SURF,NO3SURF,NH4SURF,SIO4SURF
!
!     -----SUSPENDED INORGANIC MATTER PROFILE-----
!
      real(RLEN),public,dimension(KB-1) :: ISM
!
!     ------FREQUENCY OF AVERAGING FOR OUTPUTi (IN DAYS)-----
!
      integer(ilong)                    :: savef,nitend
!
      real(RLEN),public                 :: deltat
!
!
!     -----THESE ARE THE PATHWAYS FOR THE IC AND FORCING FILES (READ TROUGH NML-----
!
     character(200)                     :: wind_input,      &
                                           surfaceS_input,  &
                                           radiance_input,  &
                                           ism_input,       &
                                           Sal_input,       &
                                           Temp_input,      &
                                           Sprofile_input,  &
                                           Tprofile_input,  &
                                           heat_input,      &
                                           surfNut_input,   &
                                           read_restart
!
 end module Service
!
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
