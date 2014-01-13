
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                          DECK MLDPTH                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MLDPTH(ZZ,T,KB,ZZMLD)
!
      use global_mem,ONLY: RLEN
!
      IMPLICIT NONE
!     .. Parameters ..
      REAL(RLEN) :: ZERO
      PARAMETER (ZERO=1.e-6)
!     ..
!     .. Scalar Arguments ..
      REAL(RLEN) :: ZZMLD
      INTEGER :: KB
!     ..
!     .. Array Arguments ..
      REAL(RLEN) :: T(KB),ZZ(KB)
!     ..
!     .. Local Scalars ..
      INTEGER :: K
!     ..
      DO K = 1,KB - 1
          IF (T(1).GT. (T(K)+.2)) GO TO 100
          ZZMLD = ZZ(K) - (T(K) + 0.2 - T(1))* (ZZ(K) - ZZ(K+1))/ 
     &         (T(K)-T(K+1)+ZERO)
  100 END DO
      RETURN
      END
