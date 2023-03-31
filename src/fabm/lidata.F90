      subroutine lidata(nlt,nchl)

      USE ogs_bfm_shared, ONLY: lam,lam1,lam2,aw,bw,bbw,ac,ac_ps,bc,bbc,acdom,apoc,bpoc,bbpoc
      USE fabm_types,     ONLY: rk
      IMPLICIT NONE 
!  Reads in radiative transfer data: specifically 
!  water data (seawater absorption, total and back scattering coefficients,
!  and chl-specific absorption and total scattering data for 
!  several phytoplankton groups). CDOM and POC can be read from files
!  or computed from parameters
!  PAR (350-700) begins at index 3,  and ends at index 17.
      integer, intent(IN) :: nlt,nchl 
      character*80 title
      character*80 cfle
      character cacbc*11,cabw*20,cacbpoc*10,cacdom*11
      double precision lambda1,lambda2,saw,sbw,sbbw,sac,sac_ps,sbc,sbb,sapoc,sbpoc,sbbpoc,sacdom
      character*4 cdir
      data cdir /'bcs/'/
      data cacbc,cabw,cacbpoc,cacdom /'acbc25b.dat','abw25_boundaries.dat','poc25b.dat','cdom25b.dat'/
      integer    :: i, n, nl,lambda

!  Water data files, a, b and bb in m-1
      cfle = cdir//cabw
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,5
       read(4,'(a80)')title
       write(6,'(a80)')title
      enddo
      do nl = 1,nlt
       read(4,20)lambda, lambda1, lambda2, saw,sbw,sbbw
       write(6,20)lambda, lambda1, lambda2, saw,sbw,sbbw
       lam(nl) = REAL(lambda,rk)
       lam1(nl) = lambda1
       lam2(nl) = lambda2
       aw(nl) = saw
       bw(nl) = sbw
       bbw(nl) = sbbw
       write(*,*)lam(nl)
      enddo
      close(4)
20    format(i5,f8.1,f8.1,f14.4,f8.4,f9.5)
 
!  Phytoplankton group chl-specific absorption and c-specific scattering 
!  data.  Chl-specific absorption data is normalized to 440 nm (?);
!  convert here to actual ac*(440)

      cfle = cdir//cacbc
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,6
       read(4,'(a80)')title
      enddo
      do n = 1,nchl
       read(4,'(a80)')title
       do nl = 1,19
        read(4,30)lambda,sac,sac_ps,sbc,sbb
        ac(n,nl) = sac
        ac_ps(n,nl) = sac_ps
        bc(n,nl) = sbc
        bbc(n,nl) = sbb
        write(*,*)lambda,ac(n,nl),ac_ps(n,nl),bc(n,nl),bbc(n,nl)
       enddo
       do nl = 20,nlt
        ac(n,nl) = 0.0
        ac_ps(n,nl) = 0.0
        bc(n,nl) = 0.0
        bbc(n,nl) = 0.0
       enddo
      enddo
      close(4)
30    format(i4,2f10.4,1f10.5,1f10.4)
!30    format(i4,4f10.4)      

!  POC absorption, scattering and back scattering normalized to mgC/m3
      cfle = cdir//cacbpoc
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,6
       read(4,'(a80)')title
      enddo
      do nl = 1,nlt
       read(4,40)lambda,sapoc,sbpoc,sbbpoc
       write(*,*) lambda,sapoc,sbpoc,sbbpoc
       apoc(nl) = sapoc
       bpoc(nl) = sbpoc
       bbpoc(nl) = sbbpoc
      enddo
      close(4)
40    format(i5,3f10.2)

!  CDOM absorption, m2 mgC-1
      cfle = cdir//cacdom
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,6
       read(4,'(a80)')title
      enddo
      do nl = 1,nlt
       read(4,50)lambda,sacdom
       write(*,*) lambda,sacdom
       acdom(1,nl) = sacdom
       acdom(2,nl) = sacdom
       acdom(3,nl) = sacdom       
      enddo
      close(4)
50    format(i5,1E14.6)

      
      return
      end
