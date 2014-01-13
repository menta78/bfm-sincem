
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                        DECK PROFU                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE PROFU(DT2)
      use global_mem,ONLY: RLEN
      use POM,ONLY: H,KB,A,C,KM,DZ,DZZ,VH,VHP,WUSURF,UF,UB,VB,WUBOT, &
                    UMOL,Z,ZZ
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!         DT2*(KM*U')' - U= -UB                                        *
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   85 DH = H
      DO K = 2,KB - 1
          A(K-1) = -DT2* (KM(K)+UMOL)/ (DZ(K-1)*DZZ(K-1)*DH*DH)
          C(K) = -DT2* (KM(K)+UMOL)/ (DZ(K)*DZZ(K-1)*DH*DH)
      END DO
      VH(1) = A(1)/ (A(1)-1.)
      VHP(1) = (-DT2*WUSURF/ (-DZ(1)*DH)-UF(1))/ (A(1)-1.)
      DO K = 2,KB - 2
          VHP(K) = 1./ (A(K)+C(K)* (1.-VH(K-1))-1.)
          VH(K) = A(K)*VHP(K)
          VHP(K) = (C(K)*VHP(K-1)-UF(K))*VHP(K)
      END DO
      CBC = MAX(.0025,.16/LOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
      CBC = CBC*SQRT(UB(KB-1)**2+ (.25* (VB(KB-1)+VB(KB-1)+VB(KB-1)+ &
           VB(KB-1)))**2)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        TO RESTORE BOTTOM B.L. DELETE NEXT LINE                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      UF(KB-1) = (C(KB-1)*VHP(KB-2)-UF(KB-1))/ &
                (CBC*DT2/ (-DZ(KB-1)*DH)-1.- (VH(KB-2)-1.)*C(KB-1))
      DO K = 2,KB - 1
          KI = KB - K
          UF(KI) = VH(KI)*UF(KI+1) + VHP(KI)
      END DO
   92 WUBOT = -CBC*UF(KB-1)
      DO K = 1,KB
          VH(K) = 0.
          VHP(K) = 0.
          A(K) = 0.
          C(K) = 0.
      END DO

      RETURN
      END
