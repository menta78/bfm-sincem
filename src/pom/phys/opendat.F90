   subroutine opendat(WSU1,WSV1,SSS1,SLUX1,ISM1,SINIZ1,TINIZ1,SINIZ,&
                 TINIZ,SWRAD1,WTSURF1,QCORR1,NO3_1,NH4_1,PO4_1,&
                 SIO4_1)
! 
! !USES:
      use global_mem,ONLY:RLEN,LOGUNIT,error_msg_prn,NML_OPEN,NML_READ
      use Service, only: wind_input,surfaceS_input,radiance_input,&
                         ism_input,Sal_input,Temp_input,Sprofile_input,&
                         Tprofile_input,heat_input,surfNut_input,ISM
!    
      use pom, ONLY: KB
      IMPLICIT NONE 
      INTEGER :: RLENGTH
      real(RLEN) :: WSU1,WSV1,SSS1,SLUX1,SWRAD1,WTSURF1,QCORR1, &
                    NO3_1,NH4_1,PO4_1,SIO4_1
      integer,parameter  :: namlst=10
      REAL(RLEN),dimension(KB) :: SINIZ1,TINIZ1,SINIZ,TINIZ
      INTEGER,PARAMETER :: N_COMP=KB-1 
      REAL(RLEN),dimension(N_COMP) :: ISM1

       namelist /pom_input/ wind_input,surfaceS_input,radiance_input, &
                          ism_input,Sal_input,Temp_input,Sprofile_input, &
                          Tprofile_input,heat_input,surfNut_input
!
       open(namlst,file='pom_input.nml',status='old',action='read',err=100)
       read(namlst,nml=pom_input, err=102)
       close(namlst)
       inquire(IOLENGTH=rlength) WSU1,WSV1
       open(11,file=wind_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 11 done'
       inquire(IOLENGTH=rlength) SSS1
       open(13,file=surfaceS_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 13 done'
       inquire(IOLENGTH=rlength) SLUX1
       open(17,file=radiance_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 17 done'
       inquire(IOLENGTH=rlength) ISM1(1)
       open(19, file=ism_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 19 done'
       inquire(IOLENGTH=rlength) SINIZ1(1)
       open(20, file=Sal_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 20 done'
       inquire(IOLENGTH=rlength) TINIZ1(1)
       open(15, file=Temp_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 15 done'
       inquire(IOLENGTH=rlength) SINIZ(1)
       open(29,file=Sprofile_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 29 done'
       inquire(IOLENGTH=rlength) TINIZ(1)
       open(10,file=Tprofile_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 10 done'
       inquire(IOLENGTH=rlength) SWRAD1,WTSURF1,QCORR1
       open(21, file=heat_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 21 done'
       inquire(IOLENGTH=rlength) NO3_1,NH4_1,PO4_1,SIO4_1
       open(18, file=surfNut_input, form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 18 done'
       write(6,*) 'open units done'
!
     return
!
100   call error_msg_prn(NML_OPEN,"opendat.F90","pom_input.nml")
102   call error_msg_prn(NML_READ,"opendat.F90","pom_input")
      end subroutine opendat
!EOC
