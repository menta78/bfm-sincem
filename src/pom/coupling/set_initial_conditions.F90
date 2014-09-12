      subroutine set_initial_conditions

!   Momme Butenschon, May 2004
!   Dipartimento di Fisica
!   Universita di Bologna
!
!      use global_mem, ONLY:RLEN
      use global_mem
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

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Definition of Initial Pelagic (D3) state variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     real(RLEN) :: O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, B1n0, B1p0, &
      P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, P3n0, P3p0, &
      P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, Z4p0, Z5c0, &
      Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, R6c0, R6n0, &
      R6p0, R6s0, O3c0, O3h0


     namelist /bfm_init_nml/ O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, &
      B1n0, B1p0, P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, &
      P3n0, P3p0, P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, &
      Z4p0, Z5c0, Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, &
      R6c0, R6n0, R6p0, R6s0, O3c0, O3h0

#ifdef INCLUDE_BEN
     real(RLEN) :: G3c0, G3h0, G13c0, G13h0, G23c0, G23h0, Y1c0, Y1n0, Y1p0, &
      Y2c0, Y2n0, Y2p0, Y3c0, Y3n0, Y3p0, Y4c0, Y4n0, Y4p0, Y5c0, Y5n0, Y5p0, &
      Q6c0, Q6n0, Q6p0, Q6s0, Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, H1c0, H1n0, &
      H1p0, H2c0, H2n0, H2p0, K1p0, K11p0, K21p0, K4n0, K14n0, K24n0, K6r0, &
      K16r0, K26r0, K3n0, K5s0, G2o0, G4n0, D1m0, D2m0, D6m0, D7m0, D8m0, D9m0


     namelist /bfm_init_nml_ben/ G3c0, G3h0, G13c0, G13h0, G23c0, G23h0, Y1c0, &
      Y1n0, Y1p0, Y2c0, Y2n0, Y2p0, Y3c0, Y3n0, Y3p0, Y4c0, Y4n0, Y4p0, Y5c0, &
      Y5n0, Y5p0, Q6c0, Q6n0, Q6p0, Q6s0, Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, &
      H1c0, H1n0, H1p0, H2c0, H2n0, H2p0, K1p0, K11p0, K21p0, K4n0, K14n0, K24n0, &
      K6r0, K16r0, K26r0, K3n0, K5s0, G2o0, G4n0, D1m0, D2m0, D6m0, D7m0, D8m0, &
      D9m0
#endif

     open(NMLUNIT,file='BFM_General.nml',status='old',action='read',err=100)
     read(NMLUNIT,nml=bfm_init_nml,err=101)
#ifdef INCLUDE_BEN
     read(NMLUNIT,nml=bfm_init_nml_ben,err=102)
#endif
     close(NMLUNIT)
!
! HERE INITIAL CONDITIONS ARE HARDWIRED.
! IF I,C. ARE TO BE PROVIDED VIA FILE READING
! PUT HERE THE "READ" STATEMENT AND DELETE THE UNNECESSARY CODE.
!
!---------------------------------------------
! Initialise river nutrients
!---------------------------------------------
      pon1p = 4.6
      pon3n = 150.0
      pon4n = 21.0
      pon5s = 120.0
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
        write(6,*) 'CalcPhyto(1)',CalcPhytoPlankton(1)
        write(6,*) 'CalcPhyto(2)',CalcPhytoPlankton(2)
        write(6,*) 'CalcPhyto(3)',CalcPhytoPlankton(3)
        write(6,*) 'CalcPhyto(4)',CalcPhytoPlankton(4)
        write(6,*) 'CalcBact',CalcPelBacteria(1)
        write(6,*) 'CalcMeso(1)',CalcMesoZooPlankton(1)
        write(6,*) 'CalcMeso(2)',CalcMesoZooPlankton(2)
        write(6,*) 'CalcMicro(1)',CalcMicroZooPlankton(1)
        write(6,*) 'CalcMicro(2)',CalcMicroZooPlankton(2)
 
! IrrOPT  in the equation of Steele and light
      ib=0
      do k = 1, NO_BOXES 
        eir(k) = 4.
 
!-----------------------------------------------------
! Definition of general pelagic state variables:
!-----------------------------------------------------
 
! pelagic gases
 
        O2o(k) =O2o0
 
! pelagic nutrients (mMol /m3)
        n1p(k)=N1p0
        n3n(k)=N3n0
        n4n(k)=N4n0
        n5s(k)=N5s0

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
        O4n(k)=O4n0
        N6r(k)=N6r0
! 
! pelagic detritus  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
 
        R6c(k) = R6c0 
        R6n(k) = r6c(k)*p_nRc
        R6p(k) = r6c(k)*p_pRc
        R6s(k) = r6c(k)*p_sRc
!
        R2c(k) = R2c0

!!!         R7c(k) = 0.

! dissolved organic matter
         R1c(k) = R1c0
         R1n(k) = R1c(k)*p_nRc*0.5
         R1p(k) = R1c(k)*p_pRc*0.5

!         dd=abs(ze(k)*he)
!         d1 = 50.
!         d2 = 100.
         if(dd.gt.d1.and.dd.le.d2) then
!         R1c(3) = 1080.+(((840.-960.)/(d2-d1))*(dd-d1))
!         R1n(k) = 14.+(((10.5-12.)/(d2-d1))*(dd-d1))
!         R1p(k) = 0.81+(((0.63-0.72)/(d2-d1))*(dd-d1))
         end if 
         if(dd.gt.d2) then
!         R1c(k) = 840.
!         R1n(k) = 10.5
!         R1p(k) = 0.63    
         end if
!         r1c(k)=r1c(k)*0.5
!         r1n(k)=r1n(k)*0.5
!         r1p(k)=r1p(k)*0.5
!
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
!        d1cc=10.
        d1cc=P1c0
        d1cn=d1cc*p_nRc
        d1cp=d1cc*p_pRc
        d1cs=d1cc*p_sRc
        d1ci=d1cc*p_iRc
        dd=abs(z(k)*h)
!        if(dd.le.d1) then
      if (CalcPhytoPlankton(1)) then
        P1c(k) = d1cc
        P1n(k) = d1cn
        P1p(k) = d1cp
        P1s(k) = d1cs
        P1l(k) = d1ci
      else
        P1c(k) = 1.e-24
        P1n(k) = 1.e-24
        P1p(k) = 1.e-24
        P1s(k) = 1.e-24
        P1l(k) = 1.e-24
      end if
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
!        d1cc = 5.
        d1cc = P2c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
        dd=abs(z(k)*h)
!        if(dd.le.d1) then
      if (CalcPhytoPlankton(2)) then
        P2c(k) = d1cc
        P2n(k) = d1cn
        P2p(k) = d1cp
        P2l(k) = d1ci
      else
         P2c(k) = 1.e-24
         P2n(k) = 1.e-24
         P2p(k) = 1.e-24
         P2l(k) = 1.e-24
      end if
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
!        d1cc = 5.
        d1cc = P3c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
        dd=abs(z(k)*h)
!        if(dd.le.d1) then
      if (CalcPhytoPlankton(3)) then
        P3c(k) = d1cc
        P3n(k) = d1cn
        P3p(k) = d1cp
        P3l(k) = d1ci
      else
         P3c(k) = 1.e-24
         P3n(k) = 1.e-24
         P3p(k) = 1.e-24
         P3l(k) = 1.e-24
      end if
!
!#endif
!         endif
        if(dd.gt.d1.and.dd.le.d2) then
!        p3c(k)=d1cc+(((1.e-24-d1cc)/(d2-d1))*(dd-d1))
!        p3n(k)=d1cn+(((1.e-24-d1cn)/(d2-d1))*(dd-d1))
!        p3p(k)=d1cp+(((1.e-24-d1cp)/(d2-d1))*(dd-d1))
!        p3l(k)=d1ci+(((1.e-24-d1ci)/(d2-d1))*(dd-d1))
        endif
!#endif           

! Large phytoplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

        P4c(k) = 1.e-24
        P4n(k) = 1.e-24
        P4p(k) = 1.e-24
        P4l(k) = 1.e-24
!        d1cc = 5.
        d1cc = P4c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
        dd=abs(z(k)*h)
!        if(dd.le.d1) then
      if (CalcPhytoPlankton(4)) then
        P4c(k) = d1cc
        P4n(k) = d1cn
        P4p(k) = d1cp
        P4l(k) = d1ci
      else
         P4c(k) = 1.e-24
         P4n(k) = 1.e-24
         P4p(k) = 1.e-24
         P4l(k) = 1.e-24
      end if
!-----------------------------------------------------
! State variables for mesozooplankton model
!-----------------------------------------------------
! carnivorous mesozooplankton ( mg C/m3 )
      if (CalcMesoZooPlankton(1)) then
        Z3c(k) = Z3c0
        Z3n(k) = Z3c(k)*p_nRc
        Z3p(k) = Z3c(k)*p_pRc
      else
         Z3c(k) = 1.e-24
         Z3n(k) = 1.e-24
         Z3p(k) = 1.e-24
      end if
! omnivorous mesozooplankton ( mg C/m3 )
      if (CalcMesoZooPlankton(2)) then
        Z4c(k) = Z4c0
        Z4n(k) = Z4c(k)*p_nRc
        Z4p(k) = Z4c(k)*p_pRc
      else
         Z4c(k) = 1.e-24
         Z4n(k) = 1.e-24
         Z4p(k) = 1.e-24
       end if

!-----------------------------------------------------
! State variables for microzooplankton model
!-----------------------------------------------------

! pelagic microzooplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
      if (CalcMicroZooPlankton(1)) then 
        Z5c(k) = Z5c0
!        Z5n(k) = Z5c(k)*p_nRc
!        Z5p(k) = Z5c(k)*p_pRc
         Z5n(k) = Z5c(k)*1.67d-2
         Z5p(k) = Z5c(k)*1.85d-3 
      else
         Z5c(k) = 1.e-24
         Z5n(k) = 1.e-24
         Z5p(k) = 1.e-24
      end if

! heterotrophic flagellates (respectively mg C/m3 mMol N/m3 mMOL P/m3)
      if (CalcMicroZooPlankton(2)) then
         Z6c(k) = Z6c0
!        Z6n(k) = Z6c(k)*p_nRc
!        Z6p(k) = Z6c(k)*p_pRc
         Z6n(k) = Z6c(k)*1.67d-2
         Z6p(k) = Z6c(k)*1.85d-3
      else
         Z6c(k) = 1.e-24
         Z6n(k) = 1.e-24
         Z6p(k) = 1.e-24
      end if
!-----------------------------------------------------
! State variables for pelagic bacteria model B1
!-----------------------------------------------------
! pelagic bacteria  (respectively mg C/m3 mMol N/m3 mMOL P/m3)

      if (CalcPelBacteria(1)) then 
        B1c(k) = B1c0
        B1n(k) = B1c(k)*p_nRc
        B1p(k) = B1c(k)*p_pRc
      else
         B1c(k) = 1.e-24
         B1n(k) = 1.e-24
         B1p(k) = 1.e-24
      end if

      enddo
!-----------------------------------------------------
! State variables for the benthic modules
!-----------------------------------------------------
      iup=KBE-1
#ifdef INCLUDE_BEN
! zoobenthos
            Y1c(1) = Y1c0
            Y2c(1) = Y2c0
            Y3c(1) = Y3c0
            Y4c(1) = Y4c0
            Y5c(1) = Y5c0

             H1c(1) = H1c0
             H2c(1) = H2c0
!
! benthic nutrients
!
             K5s(1) = 20.75
!             k5s(1)=n5s(iup)*Depth(iup)
!             k5s(1) = K5s0
              K6r(1) = K6r0
!             k4n(1)=n4n(iup)*Depth(iup)
             k4n(1)=K4n0
             K14n(1) = K14n0
             k24n(1) = K24n0
!              k1p(1)=N1p(iup)*Depth(iup)
!             K11p(1) = 1.E-15
!             K21p(1) = 1.E-15
              k1p(1)=K1p0
             K11p(1) = K11p0
             K21p(1) = K21p0
!            k3n(1)=n3n(iup)*Depth(iup)
             K3n(1) = K3n0
!             k13n(1) = 1.E-15
!            K23n(1) = 1.E-15

! benthic detritus  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
!            Q1c(1) = R1c(iup)*Depth(iup)
            Q1c(1) = Q1c0
!!!            Q1n(1) = R1n(iup)*Depth(iup)
!!!            Q1p(1) = R1p(iup)*Depth(iup)
            Q11c(1)=Q11c0

#ifdef IMPFLUX
            Q6c(1) = 1.E9 
            Q6n(1) = 1.E9
            Q6p(1) = 1.E9
            Q6s(1) = 1.E9
#else
!             Q6c(1)=r6c(iup)*Depth(iup)
!             Q6n(1)=r6n(iup)*Depth(iup)
!             Q6p(1)=r6p(iup)*Depth(iup)
!             Q6s(1)=r6s(iup)*Depth(iup)

            Q6c(1) = Q6c0
            Q6n(1) = Q6n0
            Q6p(1) = Q6p0
            Q6s(1) = Q6s0

#endif
! gases
            G2o(1) = G2o0
            G3c(1) = G3c0
!           G4n(1) = 37.0
!
! layers

             D1m(1) = D1m0
             D2m(1) = D2m0
             D6m(1) = D6m0
             D7m(1) = D7m0
             D8m(1) = D8m0
             D9m(1) = D9m0

#endif
      return
 100 call error_msg_prn(NML_OPEN,"InitParam.f90","BFM_General.nml")
 101 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters")
 102 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters_ben")

      end
