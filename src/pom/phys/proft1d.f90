      SUBROUTINE PROFT1D(F,FB,WFSURF,SWRAD,FSURF,NBC,DT2)
!
!   Momme Butenschon, February 2005
!   Dipartimento di Fisica
!   Universita di Bologna
!
!
!	INCLUDE 'comblk98.h'
      use global_mem,ONLY: RLEN
      use POM,ONLY:KB,H,UMOL,DZ,DZZ,KH,Z
      implicit none
      integer,parameter :: kbm1=kb-1
      integer,parameter :: kbm2=kb-2
      REAL(RLEN) :: f,rad,r,ad1,ad2,WFSURF,DH,SWRAD,FSURF,DT2,ee,gg,fb 
      REAL(RLEN) :: a,c,umolpr
      DIMENSION F(KB),A(KB),C(KB),ee(KB),gg(KB),fb(KB)
      DIMENSION RAD(KB),R(5),AD1(5),AD2(5)  
      integer :: ntp,k,nbc,ki,i,j
!
! Irradiance parameters after Paulson and Simpson, JPO, 1977, 952-956.
!
      NTP=2
!       NTP         =     1      2       3       4       5
!   JERLOV TYPE     =     I      IA      IB      II      III
      DATA R   /        .58 ,   .8  ,  .67  ,  .77  ,  .78   /
      DATA AD1 /        .35 ,   .2  ,  1.0  ,  1.5  ,  1.4   /
      DATA AD2 /        23. ,   5.88235  ,  17.  ,  14.  ,  7.9   /
!
      UMOLPR=1.e-5
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!         DTI2*(KH*F')'-F=-FB                                          *
!                                                                      *
!***********************************************************************
      DH=H
      DO 20 K=2,KBM1
      A(K-1)=-DT2*(KH(K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*DH &
          *DH)
      C(K)=-DT2*(KH(K)+UMOLPR)/(DZ(K)*DZZ(K-1)*DH &
          *DH)
   20 CONTINUE
!-----------------------------------------------------------------------
!   NBC=1: SURF. B.C. IS WFSURF. NO SW RADIATIVE PENETRATION.
!   NBC=2; SURF. B.C. IS WFSURF. SWRAD PENETRATES WATER COLUMN
!   NBC=3; SURF. B.C. IS TSURF. NO SW RADIATIVE PENETRATION.
!   NBC=4; SURF. B.C. IS TSURF. SWRAD PENETRATES WATER COLUMN
!
! NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE NEG. VALUES WHEN WATER COLUMN IS
!            WARMING.
!-----------------------------------------------------------------------
!------------------------------------------------------------------
!     Penetrative Radiation Calculation. At the bottom any 
!     unattenuated radiation is deposited in the bottom layer.
!------------------------------------------------------------------
        DO 512 K=1,KB
  512   RAD(K)=0.
      IF(NBC.EQ.2.OR.NBC.EQ.4.) THEN   
        DO 511 K=1,KBM1
        RAD(K)=SWRAD &
          *(   R(NTP)*EXP(Z(K)*DH/AD1(NTP)) &
                  +(1.-R(NTP))*EXP(Z(K)*DH/AD2(NTP))   )
  511   CONTINUE
      ENDIF
      GO TO (50,51,52,52), NBC
   50 CONTINUE
      EE(1)=A(1)/(A(1)-1.0)
      GG(1)=-DT2*(WFSURF+SWRAD)/(-DZ(1)*DH)-FB(1)
  500 GG(1)=GG(1)/(A(1)-1.0)
      GO TO 53
!
   51 CONTINUE
      EE(1)=A(1)/(A(1)-1.0)
      GG(1)=DT2*(WFSURF  &
           +RAD(1)-RAD(2)) &
          /(DZ(1)*DH)-FB(1)
  510 GG(1)=GG(1)/(A(1)-1.0)
      GO TO 53
!
   52 CONTINUE
      EE(1)=0.
  520 GG(1)=FSURF
!----------------------------------------------------------------------
   53 CONTINUE
!
!----------------------------------------------------------------------
      DO 101 K=2,KBM2
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1))-1.)
      EE(K)=A(K)*GG(K)
      GG(K)=(C(K)*GG(K-1)-FB(K) &
          +DT2*(RAD(K)-RAD(K+1))/(DH*DZ(K)))*GG(K)
  101 CONTINUE
!-----  BOTTOM ADIABATIC B.C. ------------------------------------------
  102 F(KBM1)=((C(KBM1)*GG(KBM2)-FB(KBM1) &
             +DT2*(RAD(KBM1)-RAD(KB))/(DH*DZ(KBM1))) &
               /(C(KBM1)*(1.-EE(KBM2))-1.))
!----------------------------------------------------------------------
      DO 105 K=2,KBM1
      KI=KB-K
      F(KI)=(EE(KI)*F(KI+1)+GG(KI))
  105 CONTINUE
!
      RETURN
      END
