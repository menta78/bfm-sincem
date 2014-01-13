      subroutine set_initial_conditions
!   Momme Butenschon, May 2004
!   Dipartimento di Fisica
!   Universita di Bologna
!
      use global_mem, ONLY:RLEN
      use Mem 
      use POM
      use Service 
      implicit none
      integer :: k,ib,iup
      real(RLEN) :: silt
      integer,parameter :: KBE=KB
      real(RLEN) :: dd, d1,d2,d1cc,d1cn,d1cp, d1cs, d1ci
      real(RLEN),parameter :: p_nRc=0.0126,p_pRc=0.7862e-3,p_sRc=0.0118 &
                              ,p_iRc=1./25.
!
! HERE INITIAL CONDITIONS ARE HARDWIRED.
! IF I,C. ARE TO BE PROVIDED VIA FILE READING
! PUT HERE THE "READ" STATEMENT AND DELETE THE UNNECESSARY CODE.
!
!---------------------------------------------
! Initialise river nutrients
!---------------------------------------------
      pon1p = 4.6
!      pon3n = 150.0
!      pon4n = 21.0
!      pon5s = 120.0
!
!-------------------------------
! Initialization of sediment concentration (kg/m3)
!
!       IN SEDIMENT TRANSPORT MODEL
!
      silt = 0.001
!
!       IN BFM (SEDIMENT CONC. IN MG/M3)
!
      do k = 1 , NO_BOXES
         ess(k) = silt*1.E6
      end do
!
!-----------------------------------------------------
! Definition of biogeochemical  global variables
!-----------------------------------------------------
 
! IrrOPT  in the equation of Steele and light
      ib=0
      do k = 1, NO_BOXES 
        eir(k) = 4.
 
!-----------------------------------------------------
! Definition of general pelagic state variables:
!-----------------------------------------------------
 
! pelagic gases
 
        O2o(k) = 280.0
!         G2o=O2o(NO_BOXES)*dz(NO_BOXES)*h
 
! pelagic nutrients (mMol /m3)
        n1p(k)=0.036
        n3n(k)=10.7
        n4n(k)=1.08
        n5s(k)=6.78
        dd=abs(z(k)*h)
        d1=40.
        if(dd.le.d1) then         
!         N1p(k) = 0.025
!         N3n(k) = 1.
!         N4n(k) = 0.25
!         N5s(k) = 1.0
         endif
        d2=200.
        if(dd.gt.d1.and.dd.le.d2) then
!         n1p(k)=0.025+(((0.1-0.025)/(d2-d1))*(dd-d1))
!         n3n(k)=1.+(((4.0-1.)/(d2-d1))*(dd-d1))
!         n5s(k)=1.+(((4.-1.0)/(d2-d1))*(dd-d1))
        endif 
        d1=200.
        d2=400.
        if(dd.gt.d1.and.dd.le.d2) then
!         n1p(k)=0.1+(((0.075-0.1)/(d2-d1))*(dd-d1))
!         n3n(k)=4.00+(((3.0-5.0)/(d2-d1))*(dd-d1))
        endif
!
!
        O4n(k)=5.
        N6r(k)=0.03
! 
! pelagic detritus  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
 
        R6c(k) = 1.0
        R6n(k) = r6c(k)*p_nRc
        R6p(k) = r6c(k)*p_pRc
        R6s(k) = r6c(k)*p_sRc
!
        R2c(k) = 200.
!!!         R7c(k) = 0.

! dissolved organic matter
!         R1c(k) = 1080.
!         R1n(k) = 14.
!         R1p(k) = 0.81
!         dd=abs(ze(k)*he)
!         d1 = 50.
!         d2 = 100.
!         if(dd.gt.d1.and.dd.le.d2) then
!         R1c(k) = 1080.+(((840.-960.)/(d2-d1))*(dd-d1))
!         R1n(k) = 14.+(((10.5-12.)/(d2-d1))*(dd-d1))
!         R1p(k) = 0.81+(((0.63-0.72)/(d2-d1))*(dd-d1))
!         end if 
!         if(dd.gt.d2) then
!         R1c(k) = 840.
!         R1n(k) = 10.5
!         R1p(k) = 0.63    
!         end if
!         r1c(k)=r1c(k)*0.5
!         r1n(k)=r1n(k)*0.5
!         r1p(k)=r1p(k)*0.5
!
        R1c(k) = 10.
        R1n(k) = R1c(k)*p_nRc*0.5
        R1p(k) = R1c(k)*p_pRc*0.5
!-----------------------------------------------------
! State variables for phytoplankton model
!-----------------------------------------------------
 
! pelagic diatoms  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
        P1c(k)=1.e-24
        p1n(k)=1.e-24
        p1p(k)=1.e-24
        p1s(k)=1.e-24
        p1l(k)=1.e-24
        d1=100.
        d1cc=10.
        d1cn=d1cc*p_nRc
        d1cp=d1cc*p_pRc
        d1cs=d1cc*p_sRc
        d1ci=d1cc*p_iRc
        dd=abs(z(k)*h)
!        if(dd.le.d1) then
        P1c(k) = d1cc
        P1n(k) = d1cn
        P1p(k) = d1cp
        P1s(k) = d1cs
        P1l(k) = d1ci
!        endif
        d2=800.
        if(dd.gt.d1.and.dd.le.d2) then
!         p1c(k)=d1cc+(((1.e-24-d1cc)/(d2-d1))*(dd-d1))
!         p1n(k)=d1cn+(((1.e-24-d1cn)/(d2-d1))*(dd-d1))
!         p1p(k)=d1cp+(((1.e-24-d1cp)/(d2-d1))*(dd-d1))
!         p1s(k)=d1cs+(((1.e-24-d1cs)/(d2-d1))*(dd-d1))
!         p1l(k)=d1ci+(((1.e-24-d1ci)/(d2-d1))*(dd-d1)) 
         endif
! pelagic flagellates  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
        P2c(k) = 1.e-24
        P2n(k) = 1.e-24
        P2p(k) = 1.e-24
        P2l(k) = 1.e-24
        d1cc = 5.
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_nRc*.5
        dd=abs(z(k)*h)
!        if(dd.le.d1) then
        P2c(k) = d1cc
        P2n(k) = d1cn
        P2p(k) = d1cp
        P2l(k) = d1ci
!        endif
        if(dd.gt.d1.and.dd.le.d2) then
!        p2c(k)=d1cc+(((1.e-24-d1cc)/(d2-d1))*(dd-d1))
!        p2n(k)=d1cn+(((1.e-24-d1cn)/(d2-d1))*(dd-d1))
!        p2p(k)=d1cp+(((1.e-24-d1cp)/(d2-d1))*(dd-d1))
!        p2l(k)=d1ci+(((1.e-24-d1ci)/(d2-d1))*(dd-d1))
        endif
 
! picophytoplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
 
        P3c(k) = 1.e-24
        P3n(k) = 1.e-24
        P3p(k) = 1.e-24
        P3l(k) = 1.e-24
        d1cc = 5.
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
        dd=abs(z(k)*h)
!        if(dd.le.d1) then
        P3c(k) = d1cc
        P3n(k) = d1cn
        P3p(k) = d1cp
        P3l(k) = d1ci
!#ifdef NO_BACT
!          P3c(k) = 1.e-24
!          P3n(k) = 1.e-24
!          P3p(k) = 1.e-24
!          P3l(k) = 1.e-24
!
!#endif
!         endif
        if(dd.gt.d1.and.dd.le.d2) then
!        p3c(k)=d1cc+(((1.e-24-d1cc)/(d2-d1))*(dd-d1))
!        p3n(k)=d1cn+(((1.e-24-d1cn)/(d2-d1))*(dd-d1))
!        p3p(k)=d1cp+(((1.e-24-d1cp)/(d2-d1))*(dd-d1))
!        p3l(k)=d1ci+(((1.e-24-d1ci)/(d2-d1))*(dd-d1))
        endif
!#ifdef NO_BACT
!         P3c(k) = 1.e-24
!         P3n(k) = 1.e-24
!         P3p(k) = 1.e-24
!         P3l(k) = 1.e-24
!
!#endif           
! pelagic inedibles  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
        P4c(k) = 1.
        P4n(k) = 0.
        P4p(k) = 0.
        P4l(k) = 0.
 
!-----------------------------------------------------
! State variables for mesozooplankton model
!-----------------------------------------------------
 
! carnivorous mesozooplankton ( mg C/m3 )
        Z3c(k) = 5.
        Z3n(k) = Z3c(k)*p_nRc
        Z3p(k) = Z3c(k)*p_pRc
! omnivorous mesozooplankton ( mg C/m3 )
        Z4c(k) = 5.
        Z4n(k) = Z4c(k)*p_nRc
        Z4p(k) = Z4c(k)*p_pRc
 
!-----------------------------------------------------
! State variables for microzooplankton model
!-----------------------------------------------------
 
! pelagic microzooplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
 
        Z5c(k) = 5.
!        Z5n(k) = Z5c(k)*p_nRc
!        Z5p(k) = Z5c(k)*p_pRc
         Z5n(k) = Z5c(k)*1.67d-2
         Z5p(k) = Z5c(k)*1.85d-3
         
!#ifdef NO_BACT
!         Z5c(k)=1.e-24
!         Z5n(k)=1.e-24
!         Z5p(k)=1.e-24
!#endif
! heterotrophic flagellates (respectively mg C/m3 mMol N/m3 mMOL P/m3)
 
        Z6c(k) = 5.
!        Z6n(k) = Z6c(k)*p_nRc
!        Z6p(k) = Z6c(k)*p_pRc
        Z6n(k) = Z6c(k)*1.67d-2
        Z6p(k) = Z6c(k)*1.85d-3
 
!#ifdef NO_BACT
!        Z6c(k)=1.e-24
!        Z6n(k)=1.e-24
!        Z6p(k)=1.e-24
!#endif
!-----------------------------------------------------
! State variables for pelagic bacteria model B1
!-----------------------------------------------------
! pelagic bacteria  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
 
        B1c(k) = 5.
        B1n(k) = B1c(k)*p_nRc
        B1p(k) = B1c(k)*p_pRc
#ifdef NO_BACT
        B1c(k) = 1.e-24
        B1n(k) = 1.e-24
        B1p(k) = 1.e-24
#endif
!
      enddo
!-----------------------------------------------------
! State variables for the benthic modules
!-----------------------------------------------------
      iup=KBE-1
!!!             H1c(1) = 120.
!!!             H1n(1)=2.
!!!             H1p(1)=.15
!!!             H2c(1) = 120.
!!!            H2n(1)=2.
!!!             H2p(1)=.15
! benthic detritus  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
!            Q1c(1) = 19.4193
!!!            Q1c(1) = R1c(iup)*Depth(iup)
!!!            Q1n(1) = R1n(iup)*Depth(iup)
!!!            Q1p(1) = R1p(iup)*Depth(iup)
!            Q1n(1) = 0.253463
!            Q1p(1) = .018929
!            Q1s(1)=0.
!!!            Q11c(1)=1.e-24
!!!            Q11n(1)=1.e-24
!!!            Q11p(1)=1.e-24
!!!            Q11s(1)=0.
!!!            Q6c(1) = 1000.
!            Q6n(1) = 13.6147
!            Q6p(1) = 1.01714
!            Q6s(1) = 88.1956
!!!             q6c(1)=r6c(iup)*Depth(iup)
!!!             q6n(1)=r6n(iup)*Depth(iup)
!!!             q6p(1)=r6p(iup)*Depth(iup)
!!!             q6s(1)=r6s(iup)*Depth(iup)
!             Q9c(1) = 22930.0
!             Q9n(1) = 158.0
!             Q9p(1) = 10.21
!             Q9s(1)=0.
!             Q17c(1) = 49.3475
!             Q17n(1) = 0.386752
!             Q17p(1) = .043545
! gases
!             G2o(1) = 0.67
!            G3c(1) = 221.50
!            G4n(1) = 37.0
! benthic nutrients
!            K5s(1) = 20.75
!             k5s(1)=n5s(iup)*Depth(iup)
!             K15s(1) = 1.E-15
!            K25s(1) = 1.E-15
!!!             K6r(1) = .85
!             K26e(1) = 1.39
!             K4n(1) = 1.15
!!!             k4n(1)=n4n(iup)*Depth(iup)
!!!             K14n(1) = 1.E-15
!             k24n(1) = 1.E-15
!             K1p(1) = 9.039
!!!              k1p(1)=N1p(iup)*Depth(iup)
!!!             K11p(1) = 1.E-15
!             K21p(1) = 1.E-15
!            K3n(1) = 1.35264
!             k3n(1)=n3n(iup)*Depth(iup)
!             k13n(1) = 1.E-15
!            K23n(1) = 1.E-15
! layers
!!!             D1m(1) = 0.01
!!!             D2m(1) = 0.1
!            xd2m(1) = 0.3
!            D3m(1) = 0.164747
!            D4m(1) = 0.185556
!             D5m(1) = 0.460259
!!!             D6m(1) = .020832
!!!             D7m(1) = .017886
!!!             D8m(1) = .017438
!!!             D9m(1) = .021504
! zoobenthos
!            Y1c(1) = 30.0
!            Y2c(1) = 610.0
!            Y3c(1) = 140.0
!            Y4c(1) = 100.0
!            Y5c(1) = 30.0
!            Y1n(1) = Y1c(1)/(7.6*12.)
!            Y2n(1) = Y2c(1)/(7.6*12.)
!            Y3n(1) = Y31(1)/(7.6*12.)
!            Y4n(1) = Y4c(1)/(7.6*12.)
!            Y5n(1) = Y5c(1)/(7.6*12.)
!            Y1p(1) = Y1c(1)/(106.*12.)
!            Y2p(1) = Y2c(1)/(106.*12.)
!            Y3p(1) = Y3c(1)/(106.*12.)
!            Y4p(1) = Y4c(1)/(106.*12.)
!            Y5p(1) = Y5c(1)/(106.*12.)
!!!             Y1c(1) = 30.
!!!             Y2c(1) = 2100.
!!!             Y3c(1) = 1400.
!!!             Y4c(1) = 200.
!!!             Y5c(1) = 70.
!!!             Y1n(1) = .357
!!!             Y2n(1) = 24.99
!!!             Y3n(1) = 14.
!!!             Y4n(1) = 2.38
!!!            Y5n(1) = .833
!!!             Y1p(1) = .02373
!!!             Y2p(1) = 1.661
!!!           Y3p(1) = 1.107
!!!            Y4p(1) = .158
!!!             Y5p(1) = .055
      return
      end
