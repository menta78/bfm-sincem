
      SUBROUTINE PROFS(FF,FB,WFSURF,WADV,SWRAD,FSURF,NBC,DT2)
!
      use global_mem,ONLY: RLEN
      use POM,ONLY: UMOL,H,KB,A,C,KH,DZ,DZZ,VH,VHP,Z
      implicit none
!     .. Scalar Arguments ..
      REAL(RLEN) :: DT2,FSURF,SWRAD,WFSURF
      INTEGER :: NBC
!     ..
!     .. Array Arguments ..
      REAL(RLEN) :: FB(KB),FF(KB),WADV(KB)
!     ..
!     .. Local Scalars ..
      REAL(RLEN) :: DH,KAPPA,PR,UMOLPR
      INTEGER :: K,KI,NTP
!     ..
!     .. Local Arrays ..
      REAL(RLEN) :: EXTC(5),RAD(KB),TR(5),NADV(KB)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC EXP
!     ..
!     .. Data statements ..
      DATA KAPPA/0.4/,PR/1./
      DATA TR/.32,.31,.29,.26,.24/
      DATA EXTC/.037,.042,.056,.073,.127/
!     ..

      UMOLPR = 9.e-6
!      UMOLPR = UMOL

      NTP = 5
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     NTP =              1      2       3       4       5              c
!    JERLOV TYPE   =     I      IA      IB      II      III            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   NBC=1: SURF. B.C. IS WFSURF. NO RADIATIVE PENETRATION.             c
!   NBC=2; SURF. B.C. IS WFSURF+SWRAD*(1.-TR). WITH SWRAD*TR PENETRATION.
!   NBC=3; SURF. B.C. IS TSURF. NO SW RADIATIVE PENETRATION.           c
!   NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE NEG. VALUES WHEN WATER    c
!   COLUMN IS WARMING.                                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DH = H
      DO K = 2,KB - 1
          A(K-1) = -DT2* (KH(K)+UMOLPR)/ (DZ(K-1)*DZZ(K-1)*DH*DH)
          C(K) = -DT2* (KH(K)+UMOLPR)/ (DZ(K)*DZZ(K-1)*DH*DH)
      END DO

      GO TO (50,51,52) NBC

   50 CONTINUE
      VH(1) = A(1)/ (A(1)-1.)
      VHP(1) = -DT2*WFSURF/ (-DZ(1)*DH) - FB(1)
      VHP(1) = VHP(1)/ (A(1)-1.)
      DO K = 1,KB
          RAD(K) = 0.
      END DO
      GO TO 53

   51 CONTINUE
      NADV = -dt2*WADV
      VH(1) = A(1)/ (A(1)-1.)
      VHP(1) = -DT2* WFSURF/ (-DZ(1)*DH) - NADV(1) - FB(1)
      VHP(1) = VHP(1)/ (A(1)-1.)
!     DO K = 1,KB
!         RAD(K) = SWRAD*TR(NTP)*EXP(EXTC(NTP)*Z(K)*DH)
!     END DO
      GO TO 53

   52 CONTINUE
      VH(1) = 0.
      VHP(1) = FSURF
      DO K = 1,KB
          RAD(K) = 0.
      END DO

   53 CONTINUE

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     c
!        DT2*(KH*FF')' -FF = FB                                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NADV = -dt2*WADV
      DO K = 1,KB - 1
          FB(K) = FB(K) + NADV(K)
      END DO
      DO K = 2,KB - 2
          VHP(K) = 1./ (A(K)+C(K)* (1.-VH(K-1))-1.)
          VH(K) = A(K)*VHP(K)
          VHP(K) = (C(K)*VHP(K-1) - FB(K))*VHP(K)
!          VHP(K) = (C(K)*VHP(K-1) - NADV(K) - FB(K))*VHP(K)
      END DO
      FF(KB-1) = FB(KB-1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   INSTEAD OF MATCHING SOLUTION TO A LOWER LAYER VALUE OF F(KB-1),    c
!   ONE MAY IMPOSE AN ADIABATIC BOTTOM B.C. BY REMOVING C'S FROM       c
!   COL. 1 OF THE NEXT TWO LINES.                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FF(KB-1) = (C(KB-1)*VHP(KB-2)-FB(KB-1))/ &
                (C(KB-1)* (1.-VH(KB-2))-1.)

   99 CONTINUE
      DO K = 2,KB - 1
          KI = KB - K
          FF(KI) = VH(KI)*FF(KI+1) + VHP(KI)
      END DO

      DO K = 1,KB
          VH(K) = 0.0
          VHP(K) = 0.0
          A(K) = 0.0
          C(K) = 0.0
      END DO

      RETURN
      END
