      MODULE POM
!
      use global_mem,ONLY: RLEN
!
      implicit none
      integer,parameter        :: KB=31
      integer,parameter        :: ilong=selected_int_kind(12)
      real(RLEN),dimension(KB) :: TF,T,TB,SF,S,SB, &
                                  UF,U,UB,VF,VB,V,Q2F,Q2,Q2B,Q2LF,Q2L,&
                                  Q2LB,RHO,WT 
      real(RLEN),dimension(KB) :: Z,ZZ,DZ,DZZ,DZR, &
                                  A,C,VH,VHP
      real(RLEN),public,dimension(KB) :: GM,GH,SM,SH,KN,SPROD,BPROD,&
                                         PROD,DTEF,KM,KH,KQ,L,TINIZ,SINIZ
      real(RLEN)               :: WUSURF,WVSURF,WUBOT,WVBOT, &
                                  WTSURF,SWRAD,WSSURF,UMOL,DTI,SMOTH 
      real(RLEN)               :: H,D,DT,rlat,rlon, &
                                  TSURF,SSURF,alat,alon,cor
      integer                  :: ISHIFT,INSHIFT,IDSHIFT,IINT,&
                                  iend
      real(RLEN)               :: time, time0
      integer(ilong)           :: idays
!
      END
