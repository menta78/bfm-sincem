
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                            DECK DENS                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ...THIS SUBROUTINE COMPUTES DENSITY-1....                            c
!    T = POTENTIAL TEMPERATURE                                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE DENS(T,S,ZZ,DT,RHO,KB)

      use global_mem,ONLY: RLEN
      implicit none
!     .. Scalar Arguments ..
      REAL(RLEN)  :: DT
      INTEGER :: KB
!     ..
!     .. Array Arguments ..
      REAL(RLEN) :: RHO(KB),S(KB),T(KB),ZZ(KB)
!     ..
!     .. Local Scalars ..
      REAL(RLEN) :: CR,GRAV,P,RHOR,SR,TR
      INTEGER :: K
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Data statements ..
      DATA GRAV/9.806/
!     ..

      DO K = 1,KB - 1
          TR = T(K)
          SR = S(K)
!      Here, the (approximate) pressure is in units of bars.
          P = -GRAV*1.025*ZZ(K)*DT*0.01
          RHOR = 999.842594 + 6.793952E-2*TR - 9.095290E-3*TR**2 + &
                1.001685E-4*TR**3 - 1.120083E-6*TR**4 + &
                6.536332E-9*TR**5
          RHOR = RHOR + (0.824493-4.0899E-3*TR+7.6438E-5*TR**2- &
                8.2467E-7*TR**3+5.3875E-9*TR**4)*SR + &
                (-5.72466E-3+1.0227E-4*TR-1.6546E-6*TR**2)* &
                (ABS(SR))**1.5 + 4.8314E-4*SR**2 
!      For shallow water the pressure dependency can be neglected
!      in which case it should also be omitted in PROFQ
          CR = 1449.1 + .0821*P + 4.55*TR - .045*TR**2 + 1.34* (SR-35.)
!          RHOR = RHOR + 1.E5*P/CR**2* (1.-2.0*P/CR**2)
         RHO(K) = (RHOR-1000.)*1.E-3
      END DO
      RHO(KB) = RHO(KB-1)

      RETURN
      END
