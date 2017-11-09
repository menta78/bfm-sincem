!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: vdiff_SOS
!
!DESCRIPTION	
! 
!  This routine calculates the vertical diffusivity of BFM
!  biochemical components and integrates BFM state var's 
!  with Source Splitting (SoS) method.
!
! !INTERFACE
!
  SUBROUTINE vdiff_SOS
!
! !USES:
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
 use global_mem, ONLY:RLEN,ZERO
!
 use Service
!
 use api_bfm, ONLY:D3STATEB
!
 use constants, ONLY: SEC_PER_DAY
!
 use mem_Param, ONLY: AssignAirPelFluxesInBFMFlag,                    &
                      AssignPelBenFluxesInBFMFlag,p_small 
!
 use mem_Settling
!
 use Mem, ONLY:D3STATE,D3SOURCE,NO_D3_BOX_STATES,                     &
               ppO2o,                                                 &
               ppO3c,                                                 &
               ppN1p,ppN3n,ppN4n,ppN5s,                               &
               ppR6c,ppR6n,ppR6p,ppR6s,                               &
               ppP1c,ppP1n,ppP1p,ppP1s,ppP1l,                         &
               ppP2c,ppP2n,ppP2p,ppP2l,                               &
               ppP3c,ppP3n,ppP3p,ppP3l,                               &
               ppP4c,ppP4n,ppP4p,ppP4l,                               &
               sediR6,sediPPY,                                        &
               iiP1,iiP2,iiP3,iiP4,                                   &
               n1p,n3n,n4n,n5s,                                       &
               jsurO2o,jsurO3c,                                       &
               Depth,                                                 &
               iiPhytoPlankton,                                       &
               jbotN1p,jbotN4n,jbotN3n,                               &
               jbotO2o,jbotO3c
!
 use POM, ONLY:SMOTH,KB,H,DTI,DZR,NRT,NBCBFM,UMOLBFM,NTP,GRAV,        &
               RHO
!
!-------------------------------------------------------------------------!
!
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
!-----COUNTER & FLAGS-----
!
 integer                  :: k,m,n,nbc
!
!-----BFM STATE VAR. @ time t-DTI, t, T+DTI RESPECTIVELY-----
!
 real(RLEN)               :: fbio(KB), ffbio(KB), fbbio(KB) 
!
!-----SURFACE FLUX STORAGE FOR NUT'S O2 & CO2-----
!
 real(RLEN)               :: surflux
!
!-----BOTTOM FLUX STORAGE FOR NUT'S O2 & CO2-----
!
 real(RLEN)               :: botflux       
!
!-----RELAXATION VELOCITY FOR NUT'S-----
!
 real(RLEN)               :: trelax
!
!-----SEDIMENTATION VELOCITY-----
!
 real(RLEN)               :: sink(KB)
!
!-----TWICE THE TIME STEP-----
!
 real(RLEN)               :: dti2  
!
!
    dti2 = DTI*2.

! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Relaxation of nutrients at surface
! -=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
    trelax=NRT/SEC_PER_DAY
!
!-----LOOP OVER BFM STATE VAR'S-----
!
  do m = 1 , NO_D3_BOX_STATES
!
!    -----ZEROING-----
!
          surflux    = ZERO 
          botflux    = ZERO
          fbio(:)    = ZERO
          fbbio(:)   = ZERO
          ffbio(:)   = ZERO
          sink(:)    = ZERO
!
!         -----LOAD BFM STATE VAR.-----
!       
          do k = 1 , KB - 1
!
              fbio(k)  = D3STATE(m,k)
              fbbio(k) = D3STATEB(m,k)
!
          end do

!
          fbio(kb) =fbio(kb-1)
          fbbio(kb)=fbbio(kb-1)
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Nutrients surface fluxes: 
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
          select case ( m )
!
                 case (ppO2o)
!
                       surflux = -(jsurO2o(1)/SEC_PER_DAY)
!
                       botflux = jbotO2o(1)/SEC_PER_DAY
!
                 case (ppO3c)
!
                       surflux= -(jsurO3c(1)/SEC_PER_DAY)
!
                       botflux = jbotO3c(1)/SEC_PER_DAY
!
                 case (ppN1p)
!
                      surflux = -(PO4SURF-n1p(1))*trelax
!
                      botflux = jbotN1p(1)/SEC_PER_DAY
!
                 case (ppN3n)
!
                      surflux = -(NO3SURF-n3n(1))*trelax

                      botflux = jbotN3n(1)/SEC_PER_DAY
!
                 case (ppN4n)
!
                      surflux = -(NH4SURF-n4n(1))*trelax
!
                      botflux = jbotN4n(1)/SEC_PER_DAY
!
                 case (ppN5s)
!
                      surflux = -(SIO4SURF-n5s(1))*trelax
!
                case default
!
                      surflux = ZERO
!
          end select
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Sedimentation R6
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
          select case ( m )
!
                 case (ppR6c:ppR6s)
!
                       do k = 1 , KB - 1
!
!                         --------- Calculate sink via Stoke's Law ---------------------
!                          sink(k) =-((GRAV*(dia_R6**2.0))/(18.0*visc))*(rhos/(RHO(k)+1.0)-1.0)
!                         --------------------------------------------------------------
                           sink(k) = -sediR6(k)/SEC_PER_DAY
!
                       end do
!
                       ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                       ! Set the  detritus sinking vel. at   
                       ! the water-sediment interface 
                       ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

                       sink(kb)=-p_burvel_R6/SEC_PER_DAY
                       if(ABS(sink(kb-1)).lt.ABS(sink(kb))) &
                       sink(kb) = sink(kb-1)

!
                  end select
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Sedimentation Phytoplankton 
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          if (m.GE.ppP1c .AND. m.LE.ppP4l) then
!
!
                     select case (m)
!
                               case (ppP1c:ppP1s)
!
                                     n = iiP1
!
                               case (ppP2c:ppP2l)
!
                                     n = iiP2
!
                               case (ppP3c:ppP3l)
!
                                     n = iiP3
!
                               case (ppP4c:ppP4l)
!
                                     n = iiP4
!
                     end select
!
                         do k = 1 , KB - 1
!
                            sink(k) = -sediPPY(n,k)/SEC_PER_DAY
!
                         enddo
!
                         ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                         ! Set the plankton sinking vel. at
                         ! the water-sediment interface
                         ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

                       sink(kb)=-p_burvel_PI/SEC_PER_DAY
                       if(ABS(sink(kb-1)).lt.ABS(sink(kb))) &
                       sink(kb) = sink(kb-1)

!-----------------------------------------------------------------

!
        end if
!
!                 -----SINKING: UPSTREAM VERTICAL ADVECTION-----
!
                  call adverte(fbbio,fbio,ffbio,sink)
!
!                 -----SOURCE SPLITTING LEAPFROG INTEGRATION-----
!
      do K=1,KB-1
!
        ffbio(k)=fbbio(k)+DTI2*((ffbio(k)/H)+D3SOURCE(m,k))
!
      end do
!
!     -----COMPUTE VERTICAL DIFFUSION AND TERMINATE INTEGRATION-----
!     -----IMPLICIT LEAPFROGGING-----
!
       CALL PROFTS(ffbio,surflux,botflux,ZERO,ZERO,NBCBFM,DTI2,NTP,UMOLBFM)
!
!     -----CLIPPING......IF NEEDED-----
!
      do k=1, KB-1
!
            ffbio(k)=max(p_small,ffbio(k))
!
      end do
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Mix the time step and restore time sequence
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
     do n = 1, KB-1
!
      D3STATEB(m,n)=fbio(n)+0.5_RLEN*smoth*(ffbio(n)+ fbbio(n)-2.0_RLEN*fbio(n))
      D3STATE(m,n)=ffbio(n)
!
     enddo
!

   enddo ! loop over NO_D3_BOX_STATES
!
      end subroutine vdiff_sos
!
!
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: adverte
!
!DESCRIPTION    
!
!    SUBROUTINE TO HANDLE THE SINKING OF BFM STATE VAR'S
!    SINKING IS TREATED AS DOWNWARD VERTICAL ADVECTION
!    COMPUTED WITH UPSTREAM FINITE DIFFERENCES.
!
!                                          Marco.Zavatarelli@unibo.it
!
! !INTERFACE
!
  SUBROUTINE adverte(FB,F,FF,W)
!
! !USES:
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use POM, ONLY:KB, DZR, H
!
      use global_mem, ONLY:RLEN
!-------------------------------------------------------------------------!
!
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
 real(RLEN) :: FB(KB),F(KB),FF(KB)
 real(RLEN) :: W(KB)
 real(RLEN) :: DTI2
 integer :: k
!
       F(KB)=F(KB-1)
       FB(KB)=FB(KB-1)
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
! Calculate vertical advection. Mind downward velocities are negative!
! Upwind scheme:
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
!
      FF(1)=DZR(1)*F(1)*W(2)

      do K=2,KB-1
!
         FF(K)=DZR(K)*(F(K)*W(K+1)-F(K-1)*W(K)) 
!
      end do
!
      return
!
    end subroutine adverte 
!
!EOC
