#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Light and other environmental forcing used in the BFM
!
! !INTERFACE
   subroutine envforcing_bfm(kt)
!
! !DESCRIPTION
!
! !USES
! BFM modules
   use constants,  only: E2W
   use global_mem, only: RLEN,ZERO,LOGUNIT,ONE
   use mem_param,  only: p_small, p_atm0
   use mem_PAR,    only: ChlAttenFlag, P_PARRGB, P_PAR, &
                         R_EPS, B_EPS, G_EPS, P_EPSIR,  &
                         EIRR, EIRB, EIRG
   use mem,        only: xEPS, ESS, ETW, ESW, EWIND,    &
                         Depth, EIR, ERHO, EICE, EPR,   &
                         NO_BOXES, NO_BOXES_XY
   use api_bfm
   use SystemForcing, only : FieldRead
   use sw_tool, only: sw_t_from_pt
   use time, only: bfmtime
#ifdef INCLUDE_PELCO2
   use mem,        only: ppO3c, ppO3h, ppN6r
   use mem_CO2,    only: AtmCO20, AtmCO2, AtmSLP, patm3d
   USE trcbc,      only: sf_trcsbc, n_trc_indsbc
#endif
! OPA modules
   use oce_trc
   use trc_oce, only: etot3
   use eosbn2,  only: nn_eos
   use sbc_oce, only: atm_co2
   use sbcapr,  only: apr

IMPLICIT NONE
! OPA domain substitutions
#include "domzgr_substitute.h90"
!
! !INPUT PARAMETERS:
   integer, intent( IN ) ::  kt  ! ocean time-step index
!
! !OUTPUT PARAMETERS:

! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer             :: i,j,k,n
   real(RLEN),allocatable,dimension(:,:,:) :: zepsv,zpar
   real(RLEN),allocatable,dimension(:,:,:) :: zparb,zparg,zparr
   real(RLEN),allocatable,dimension(:,:,:) :: zepsb,zepsg,zepsr
   real(RLEN),allocatable,dimension(:,:,:) :: zetotb,zetotg,zetotr,zetotv
!EOP
!-----------------------------------------------------------------------
!BOC
 
   !---------------------------------------------
   ! Assign temperature, salinity and density
   !---------------------------------------------
   ETW(:) = pack(tsn(:,:,:,jp_tem),SEAmask)
   ESW(:) = pack(tsn(:,:,:,jp_sal),SEAmask)

   ! convert to in-situ temperature and practical salinity
   ! TEOS-10 conversions not yet available
   if ( nn_eos .eq. 0 ) then
      ETW(:)  = sw_t_from_pt(ESW,ETW,EPR,EPR*0.)
   endif

   !---------------------------------------------
   ! Assign in-situ density
   !---------------------------------------------
   ERHO(:) = pack( (rhd(:,:,:) + 1._RLEN) * rau0, SEAmask)

   !---------------------------------------------
   ! Assign wind speed
   !---------------------------------------------
   EWIND(:) = pack(wndm(:,:),SRFmask(:,:,1) )
   !---------------------------------------------
   ! Assign Sea-ice cover
   !---------------------------------------------
   EICE(:) = pack(fr_i(:,:),SRFmask(:,:,1) )

#ifdef INCLUDE_PELCO2
   !---------------------------------------------
   ! Assign atmospheric CO2 fields
   !---------------------------------------------
   !
   ! CO2 atmospheric mixing ratio
   !
   SELECT CASE ( AtmCO2%init )
     CASE (1) ! Read timeseries in BFM
       call FieldRead(AtmCO2)
     CASE (2) ! Read Boundary Conditions using NEMO fldread (namtrc_bc)
       n = n_trc_indsbc(ppO3c)
       AtmCO2%fnow = pack( sf_trcsbc(n)%fnow(:,:,1),SRFmask(:,:,1) )
     CASE (3) ! Receive Boundary Conditions from NEMO
       AtmCO2%fnow = pack( atm_co2(:,:),SRFmask(:,:,1) )
#ifdef CCSMCOUPLED
       if ( kt < 10 .AND. bfm_init /= 1 ) then                 ! ugly patch to get values in first step of CESM
          WHERE ( AtmCO2%fnow == 0. ) ; AtmCO2%fnow = AtmCO20 ; END WHERE ;
       endif
#endif
   END SELECT
   !
   ! Atmospheric sea level pressure (MFS index jp_msl  = 4)
   SELECT CASE ( AtmSLP%init )
     CASE (1) ! Read timeseries in BFM
       call FieldRead(AtmSLP)
     CASE (2) ! Read Boundary Conditions using NEMO fldread (namtrc_bc)
       n = n_trc_indsbc(ppO3h)
       AtmSLP%fnow = pack( sf_trcsbc(n)%fnow(:,:,1),SRFmask(:,:,1) )
     CASE (3) ! Receive Boundary Conditions from NEMO
       AtmSLP%fnow = pack( apr(:,:),SRFmask(:,:,1) )
#ifdef CCSMCOUPLED
       if ( kt < 10 .AND. bfm_init /= 1 ) then                 ! ugly patch to get values in first step of CESM
         WHERE ( AtmSLP%fnow == 0. ) ; AtmSLP%fnow = p_atm0 ; END WHERE ;
       endif
#endif
   END SELECT
   ! broadcast surface atm pressure over the water column (NO_BOXES)
   patm3d = pack(SPREAD(unpack(AtmSLP%fnow,SRFmask(:,:,1),ZERO),DIM=3,Ncopies=SIZE(SEAmask,3)),SEAmask)

   ! print control on received fluxes
   IF ( (kt-nit000)<20 .OR. MOD(kt,200)==0 .OR. kt==nitend) THEN
      write(LOGUNIT,'(a,i14,a,f10.4,a,f10.3)') 'envforcing_bfm - Step ', kt,  & 
         '  Air_xCO2: ',maxval(AtmCO2%fnow), '  SLP:', maxval(AtmSLP%fnow) 
   ENDIF
#endif

   !---------------------------------------------
   ! Compute the light climate
   ! Note that in BFM light is defined at the
   ! top of each cell (W grid)
   ! and extinction coefficients are in the
   ! middle of the cell (T grid)
   !---------------------------------------------
   ! Update the extinction coefficient
   ! both Broadband or 3-band
   !---------------------------------------------
   call CalcVerticalExtinction( )

   !---------------------------------------------
   ! Assign surface PAR to the top layer
   ! for the BFM (convert W/m2 to PAR in uE,
   ! add parametric zero for nighttime)
   ! Initialise the bioshading array if ln_qsr_bio
   ! using irradiance (not PAR!) 
   ! also including the IR extinction
   !---------------------------------------------
   select case (ChlAttenFlag) 
   case (2) ! 3-band
      allocate(zparr(jpi,jpj,jpk))
      zparr(:,:,1) = p_PARRGB*(qsr(:,:)+p_small)/E2W 
      allocate(zparg(jpi,jpj,jpk))
      zparg(:,:,1) = zparr(:,:,1)
      allocate(zparb(jpi,jpj,jpk))
      zparb(:,:,1) = zparr(:,:,1)
      allocate(zepsb(jpi,jpj,jpk))
      zepsb(:,:,:) = unpack(B_eps(:),SEAmask,ZEROS)
      allocate(zepsr(jpi,jpj,jpk))
      zepsr(:,:,:) = unpack(R_eps(:),SEAmask,ZEROS)
      allocate(zepsg(jpi,jpj,jpk))
      zepsg(:,:,:) = unpack(G_eps(:),SEAmask,ZEROS)
      if (ln_qsr_bio) then
         allocate(zetotb(jpi,jpj,jpk))
         allocate(zetotg(jpi,jpj,jpk))
         allocate(zetotr(jpi,jpj,jpk))
         etot3(:,:,1) = (ONE-p_PAR)*qsr(:,:)   ! infrared
         zetotb(:,:,1)= p_PARRGB*qsr(:,:) ! blue
         zetotg(:,:,1)= p_PARRGB*qsr(:,:) ! green
         zetotr(:,:,1)= p_PARRGB*qsr(:,:) ! red
      end if
   case default ! broadband
      allocate(zpar(jpi,jpj,jpk))
      zpar(:,:,1) = p_PAR*(qsr(:,:)+p_small)/E2W 
      allocate(zepsv(jpi,jpj,jpk)) ! temporary for visible extinction
      zepsv(:,:,:) = unpack(xEPS(:),SEAmask,ZEROS)
      if (ln_qsr_bio) then
         allocate(zetotv(jpi,jpj,jpk))
         etot3(:,:,1) = (ONE-p_PAR)*qsr(:,:)  ! infrared
         zetotv(:,:,1)= p_PAR*qsr(:,:)        ! visible
      end if
   end select
   !---------------------------------------------
   ! Light field in the interior
   ! Distinguish the broadband and 3-band
   ! and store the array etot3 to be used by 
   ! NEMO when ln_qsr_bio
   ! Note that NEMO prescribes absorption
   ! of shortwave in the first 400 m.
   ! This is not done here and results may 
   ! slightly differ at depth when coupled
   !---------------------------------------------
   select case (ChlAttenFlag) 
   case (2) ! 3-band
      do k = 1,jpkm1
         do j = 1,jpj
            do i = 1,jpi
               zparb(i,j,k+1) = zparb(i,j,k)*exp(-zepsb(i,j,k)*fse3w(i,j,k))
               zparg(i,j,k+1) = zparg(i,j,k)*exp(-zepsg(i,j,k)*fse3w(i,j,k))
               zparr(i,j,k+1) = zparr(i,j,k)*exp(-zepsr(i,j,k)*fse3w(i,j,k))
            end do 
         end do 
      end do 
      EIRB(:) = pack(zparb(:,:,:),SEAmask)
      EIRG(:) = pack(zparg(:,:,:),SEAmask)
      EIRR(:) = pack(zparr(:,:,:),SEAmask)
      EIR(:) = EIRB(:) + EIRG(:) + EIRR(:)
      ! weighted broadband diffuse attenuation coefficient for diagnostics
      xEPS(:) = (EIRB(:)*B_eps(:) + EIRG(:)*G_eps(:) + EIRR(:)*R_eps(:))/EIR(:)
      if (ln_qsr_bio) then
         do k = 2,jpk
            do j = 1,jpj
               do i = 1,jpi
                  etot3(i,j,k) = etot3(i,j,k-1)*exp(-p_epsIR*fse3t(i,j,k-1))       ! infrared
                  zetotb(i,j,k)= zetotb(i,j,k-1)*exp(-zepsb(i,j,k)*fse3t(i,j,k-1)) ! blue
                  zetotg(i,j,k)= zetotg(i,j,k-1)*exp(-zepsg(i,j,k)*fse3t(i,j,k-1)) ! green
                  zetotr(i,j,k)= zetotr(i,j,k-1)*exp(-zepsr(i,j,k)*fse3t(i,j,k-1)) ! red
               end do 
            end do 
         end do 
         etot3(:,:,:) = etot3(:,:,:)+zetotb(:,:,:)+zetotg(:,:,:)+zetotr(:,:,:)
      end if
   case default ! broadband
      do k = 1,jpkm1
         do j = 1,jpj
            do i = 1,jpi
               zpar(i,j,k+1) = zpar(i,j,k)*exp(-zepsv(i,j,k)*fse3w(i,j,k))
            end do 
         end do 
      end do 
      EIR(:) = pack(zpar(:,:,:),SEAmask)
      if (ln_qsr_bio) then
         do k = 2,jpk
            do j = 1,jpj
               do i = 1,jpi
                  etot3(i,j,k) = etot3(i,j,k-1)*exp(-p_epsIR*fse3t(i,j,k-1))         ! infrared
                  zetotv(i,j,k)= zetotv(i,j,k-1)*exp(-zepsv(i,j,k-1)*fse3t(i,j,k-1)) ! visible
               end do 
            end do 
         end do 
         etot3(:,:,:) = etot3(:,:,:)+zetotv(:,:,:)
      end if
   end select

   !---------------------------------------------
   ! Deallocate temporary arrays
   !---------------------------------------------
   select case (ChlAttenFlag) 
   case (2) ! 3-band
      deallocate(zparb)
      deallocate(zparg)
      deallocate(zparr)
      deallocate(zepsb)
      deallocate(zepsg)
      deallocate(zepsr)
   case default
      deallocate(zpar)
      deallocate(zepsv)
   end select
   if (ln_qsr_bio) then
     select case (ChlAttenFlag) 
     case (2) ! 3-band
        deallocate(zetotb)
        deallocate(zetotr)
        deallocate(zetotg)
     case default
        deallocate(zetotv)
     end select
   end if

   end subroutine envforcing_bfm
!EOC
!-----------------------------------------------------------------------
