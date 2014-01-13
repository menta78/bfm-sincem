  SUBROUTINE vertdiff_bio_1d
!
 use global_mem, ONLY:RLEN,ZERO
 use Service
 use api_bfm, ONLY:D3STATEB
 use constants
 use mem_Param
!
#ifdef INCLUDE_BEN   
 use Mem, ONLY:D3STATE,D3SOURCE,NO_BOXES,NO_BOXES_Z,NO_D3_BOX_STATES,ppO2o, &
   ppN1p,ppN3n,ppN4n,ppN5s,ppR6c,ppR6n,ppR6p,ppR6s,ppP1c,ppP1n,ppP1p, &
   ppP1s,ppP1l,sediR6,sediPPY,iiP1,n1p,n3n,n4n,n5s,D2SOURCE_BEN
#else
 use Mem, ONLY:D3STATE,D3SOURCE,NO_BOXES,NO_BOXES_Z,NO_D3_BOX_STATES,ppO2o, &
   ppN1p,ppN3n,ppN4n,ppN5s,ppR6c,ppR6n,ppR6p,ppR6s,ppP1c,ppP1n,ppP1p, &
   ppP1s,ppP1l,sediR6,sediPPY,iiP1,n1p,n3n,n4n,n5s
#endif
 use POM, ONLY:SMOTH,KB,DZ,H,DTI,DZR,time 
!
 implicit none
 integer :: ll,k,m,llneg,n,nbc
 real(RLEN) :: fbio(KB), ffbio(KB), fbbio(KB),fmean(KB) 
 real(RLEN) :: surflux  
 real(RLEN) :: sink(KB),dti2,trelax  
!
        dti2 = DTI*2.
!        trelax=dz(1)/SEC_PER_DAY
        trelax=0.22/SEC_PER_DAY
        do m = 1 , NO_D3_BOX_STATES
          surflux = ZERO 
          fbio    = ZERO
          fbbio   = ZERO
          ffbio   = ZERO
          fmean   = ZERO
          sink    = ZERO
!          if(m.eq.2) then
!            do k = 1, kb-1
!            write(6,*), 'PO4 start loop ', D3STATE(m,k),n1p(k), D3SOURCE(m,k)
!            enddo
!          endif
          do k = 1 , KB - 1
              fbio(k) = D3STATE(m,k)
              fbbio(k) = D3STATEB(m,k)
          end do
              fbio(kb)=fbio(kb-1)
              fbbio(kb)=fbbio(kb-1)
!
! Surface fluxes Nutrients:
          select case ( m )
          case (ppO2o)
            surflux = .00
            NBC = 1
          case (ppN1p)
!            surflux = -(pon1p-n1p(1))*trelax
 !           NBC = 2
             surflux=pon1p
             NBC=2
          case (ppN3n)
!            surflux = -(pon3n-n3n(1))*trelax
!            NBC= 2
             surflux= pon3n
             NBC=2
          case (ppN4n)
!            surflux = -(pon4n-n4n(1))*trelax
!            NBC = 2
             surflux=pon4n
             NBC=2
          case (ppN5s)
!            surflux = -(pon5s-n5s(1))*trelax
!            NBC = 2
             surflux=pon5s
             NBC=2
          case default
             surflux = 0.
             NBC = 1
          end select
            
!    Sedimentation:
          select case ( m )
          case (ppR6c,ppR6n,ppR6p,ppR6s)
            do k = 1 , KB - 1
                      sink(k) = -sediR6(k)/SEC_PER_DAY
            end do

!            case (ppP1c,ppP1n,ppP1p,ppP1l,ppP1s) 
          case (ppP1c,ppP1n,ppP1p,ppP1s)
            do k = 1 , KB - 1
                      sink(k) = -sediPPY(iiP1,k)/SEC_PER_DAY
!                      PRINT *,'k,iiP1,sediPPY(k,iiP1)',k,iiP1,sediPPY(k,iiP1)
!                      sink(k)=0.            (G)
            end do
          end select
!                      PRINT*,'sink',sink

! SINKING:
!
        CALL ADVERTE(fbbio,fbio,DTI2,ffbio,DZR,sink,H,m)
!
! calculate transport due to vertical diffusion.
! 
        call profe(ffbio,surflux,ZERO,surflux,nbc,dti2,m)

!       -----COMPUTE RATE OF CHANGE DUE TO PHYSICAL PROCESSES-----
!
          do k=1,NO_BOXES_Z + 1
             ffbio(k)=max(p_small,ffbio(k))
          end do
!
!         -----MIX THE TIME STEP AND RESTORE TIME SEQUENCE-----
!
          do n = 1, NO_BOXES
 !            D3STATEB(m,n)=fbio(n)+0.5*smoth*0.5*smoth* &
!                            (ffbio(n)+ fbbio(n)-2.*fbio(n))
             D3STATEB(m,n)=fbio(n)+0.5*smoth*(ffbio(n)+ fbbio(n)-2.*fbio(n))
             D3STATE(m,n)=ffbio(n)
          enddo
                     if(m.eq.2) then
            do k = 1, kb-1
            enddo
          endif

      enddo
      end
!
      SUBROUTINE ADVERTE(FB,F,DTI2,FF,DZR,W,H,N)
!
! Momme Butenschoen, September 2006
! Universita' di Bologna
!
!    advection subroutine of POM adapted to include ecological sources
!    by source splitting mechanism
!
      use constants
      use global_mem, ONLY:RLEN
      use mem, ONLY:D3SOURCE
      use POM, ONLY:KB
      IMPLICIT NONE
      REAL(RLEN) :: FB(KB),F(KB),FF(KB)
      REAL(RLEN) :: DZR(KB),W(KB)
      REAL(RLEN) :: H,DTI2
      INTEGER :: k,n
!
       F(KB)=F(KB-1)
       FB(KB)=FB(KB-1)
!
!****** DO VERTICAL ADVECTION **********
! mind downward velocities are negative!
! central scheme:
!      FF(1)=.5E0*DZR(1)*(F(1)+F(2))*W(2)
!
!      DO K=2,KB-1
!         FF(K)=-.5E0*DZR(K)*((F(K-1)+F(K))*W(K) &
!                          -(F(K)+F(K+1))*W(K+1))
!      END DO
! upwind scheme:
      FF(1)=DZR(1)*F(1)*W(2)

      do K=2,KB-1
         FF(K)=DZR(K)*(F(K)*W(K+1) &
                          -F(K-1)*W(K))
      end do

!******STEP FORWARD IN TIME **********
!
      do K=1,KB-1
!       if(n.eq.11.and.k.eq.1) then
!        write(6,*)' D3S INT. P1C ',  FB(1),D3SOURCE(n,k), D3SOURCE(n,k)*DTI2/SEC_PER_DAY
!       endif
!       FF(K)=FB(K)+(DTI2*FF(K)/H)+D3SOURCE(n,k)*DTI2/SEC_PER_DAY
        FF(K)=FB(K)+(DTI2*FF(K)/H)+D3SOURCE(n,k)*DTI2

      end do
!
      return

!
      END

