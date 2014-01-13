
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                            DECK PROFV                                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     c
!        DT2*(KM*V')'-V= -VB                                           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE PROFV(DT2)

      use global_mem,ONLY: RLEN
      use POM,ONLY: H,KB,A,C,KM,DZ,DZZ,VH,VHP,WVSURF,VF,UB,VB,UMOL, &
                    WVBOT,Z,ZZ
      implicit none
!     .. Scalar Arguments ..
      REAL(RLEN) :: DT2
!     ..
!     .. Local Scalars ..
      REAL(RLEN) :: CBC,DH,UMOL1
      INTEGER :: K,KI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LOG,MAX,SQRT
!     ..
      UMOL1 = 0.0007
      DH = H
      DO K = 2,KB - 1
          A(K-1) = -DT2* (KM(K)+UMOL)/ (DZ(K-1)*DZZ(K-1)*DH*DH)
          C(K) = -DT2* (KM(K)+UMOL)/ (DZ(K)*DZZ(K-1)*DH*DH)
      END DO
      VH(1) = A(1)/ (A(1)-1.)
      VHP(1) = (-DT2*WVSURF/ (-DZ(1)*DH)-VF(1))/ (A(1)-1.)
   98 CONTINUE
      DO K = 2,KB - 2
          VHP(K) = 1./ (A(K)+C(K)* (1.-VH(K-1))-1.)
          VH(K) = A(K)*VHP(K)
          VHP(K) = (C(K)*VHP(K-1)-VF(K))*VHP(K)
      END DO
  104 CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
      CBC = CBC*SQRT((.25* (UB(KB-1)+UB(KB-1)+UB(KB-1)+UB(KB-1)))**2+ &
           VB(KB-1)**2)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        TO RESTORE BOTTOM B.L. DELETE NEXT LINE                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      VF(KB-1) = (C(KB-1)*VHP(KB-2)-VF(KB-1))/ &
                (CBC*DT2/ (-DZ(KB-1)*DH)-1.- (VH(KB-2)-1.)*C(KB-1))
      DO K = 2,KB - 1
          KI = KB - K
          VF(KI) = VH(KI)*VF(KI+1) + VHP(KI)
      END DO
   92 WVBOT = -CBC*VF(KB-1)
      DO K = 1,KB
          VH(K) = 0.
          VHP(K) = 0.
          A(K) = 0.
          C(K) = 0.
      END DO

      RETURN
      END
