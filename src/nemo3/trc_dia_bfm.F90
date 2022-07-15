#include "cppdefs.h"
MODULE trcdiabfm
!   !!======================================================================
!   !!                       ***  MODULE  trcdiabfm  ***
!   !! Ocean passive tracers:  save output for passives tracers
!   !!======================================================================

   IMPLICIT NONE
   
   PRIVATE  

   PUBLIC  trc_dia_bfm

#if defined key_iomput
   LOGICAL,PUBLIC    :: bfm_iomput=.TRUE. ! use xios in nemo
#else
   LOGICAL,PUBLIC    :: bfm_iomput=.FALSE.
#endif

CONTAINS
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_dia_bfm.F90
!
! !INTERFACE:
   subroutine trc_dia_bfm(kt)
!
! !DESCRIPTION:
!
! !USES:
   use oce_trc
   use global_mem, only: RLEN, LOGUNIT, bfm_lwp, SkipBFMCore
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use api_bfm,    only: out_delta, save_delta, time_delta, &
                         update_save_delta, unpad_out, &
                         stStart, stEnd, var_names, var_ids
   use time,       only: bfmtime
#ifdef INCLUDE_PELCO2
   use constants, ONLY:MW_C
   use mem_CO2, ONLY: CloseCO2
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: O3h, O3c, DIC, ALK, ERHO
#endif
#endif
  use iom, ONLY: iom_use

   implicit none
!
! !INPUT PARAMETERS:
   integer,intent(IN)     :: kt
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): Marcello Vichi (CMCC-INGV)
!  2014     : Tomas Lovato (CMCC)
!
! !LOCAL VARIABLES:
   real(RLEN)             :: localtime  !time in seconds
   integer                :: m
!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Re-initialize BFM var_ids values according to the XIOS file_def
   IF ( bfm_iomput .and. kt == nit000 ) THEN
      do m = stStart, stEnd
         var_ids(m) = -1
         !if ( iom_use( TRIM(var_names(m)) ) ) write(LOGUNIT,*) 'iomuse -> ' , TRIM(var_names(m))
         if ( iom_use( TRIM(var_names(m)) ) ) var_ids(m) = 100
      enddo
    ENDIF
   !---------------------------------------------
#ifdef INCLUDE_PELCO2
   ! Update DIC and alkalinity from model units to diagnostic output
   ! after transport of model state variables O3c and O3h
   ! mg C/m3 --> umol/kg
   ! mmol eq/m3 --> umol/kg
   DIC(:) = O3c(:)/MW_C/ERHO(:)*1000.0_RLEN
   ALK(:) = O3h(:)/ERHO(:)*1000.0_RLEN
#endif
   IF ( bfm_iomput ) THEN
      call bfm_iom(kt) 
   ELSE
      ! Update means
      call calcmean_bfm(ACCUMULATE)
   ENDIF 
   !---------------------------------------------
   ! Write diagnostic output
   !---------------------------------------------
   if ( bfmtime%stepnow .eq. save_delta ) then
      localtime = (time_delta - real(bfmtime%step0,RLEN)) * bfmtime%timestep
      IF ( .NOT. bfm_iomput ) THEN
         if ( lwp ) then
            write(numout,*) 'trc_dia_bfm : write NetCDF passive tracer concentrations at ', kt, 'time-step'
            write(numout,*) '~~~~~~ '
         end if
         call calcmean_bfm(MEAN)
         call save_bfm(localtime)
      ENDIF
      if ( unpad_out ) then
        localtime = real((bfmtime%stepnow - bfmtime%step0),RLEN) * bfmtime%timestep
        call write_rst_bfm(localtime)
      end if
      call update_save_delta(out_delta,save_delta,time_delta)
   end if
   !---------------------------------------------
   ! Save Restart and Close the files
   !---------------------------------------------
   if ( bfmtime%stepnow == bfmtime%stepEnd ) then
      if ( .NOT. unpad_out ) then
        localtime = real((bfmtime%stepEnd - bfmtime%step0),RLEN) * bfmtime%timestep
        call write_rst_bfm(localtime)
        IF ( .NOT. bfm_iomput ) call close_ncdf(ncid_bfm)
      endif
     !close systemforcings
#ifdef INCLUDE_PELCO2
      call CloseCO2()
#endif
      ! clear main memory
      IF ( .NOT. SkipBFMCore ) call ClearMem

      LEVEL1 ' '
      LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
      LEVEL1 '             EXPERIMENT FINISHED               '
      LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
      LEVEL1 ' '

   else
   !---------------------------------------------
   ! Reset the arrays for next step
   !---------------------------------------------
      call ResetFluxes
   endif
   
   return
   end subroutine trc_dia_bfm

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Create the restart NetCDF file
!
! !INTERFACE:
  subroutine write_rst_bfm(savetime)
!
! !DESCRIPTION:
! Wrapper subroutine to store restart file of BFM variables
!
! !REVISION HISTORY:
!  Original author(s): Tomas Lovato (2014)
!
! !USES:
   use global_mem, only: RLEN, LOGUNIT, bfm_lwp
   use time,       only: bfmtime
   use api_bfm,    only: out_rst_fname, parallel_rank,  &
                         ocepoint, surfpoint, botpoint
   use netcdf_bfm, only: init_netcdf_rst_bfm, save_rst_bfm, &
                         ncid_rst, close_ncdf
   use dom_oce
!
   IMPLICIT NONE
  ! * Substitutions
#include "domzgr_substitute.h90"
!
! !INPUT PARAMETERS:
  real(RLEN),intent(in)     :: savetime
! !LOCAL VARIABLES:
  character(len=PATH_MAX)   :: thistime, thisfname
  character(LEN=4)          :: str
!EOP
!-----------------------------------------------------------------------
!BOC
  LEVEL1 'write_rst_bfm: Create restart file at timestep: ',bfmtime%stepNow
  !
  ! variable parallel_rank must have been assigned previously
  ! in the coupling with the ocean model
  write(str,'(I4.4)') parallel_rank
  !
  ! Compose restart filename
  write(thistime,'(I8.8)') bfmtime%stepNow
  thisfname=TRIM(out_rst_fname)//'_'//TRIM(thistime)//'_restart_bfm_'//str

  ! Create the restart file
   call init_netcdf_rst_bfm(thisfname,TRIM(bfmtime%datestring),0,  &
             lat2d=gphit,lon2d=glamt,z=fsdept(1,1,:),  &
             oceanpoint=ocepoint,                      &
             surfacepoint=surfpoint,                   &
             bottompoint=botpoint,                     &
             mask3d=tmask)

  ! write restart data
  call save_rst_bfm(savetime)

  ! close the file
  call close_ncdf(ncid_rst)

  LEVEL1 'write_rst_bfm: restart file creation ... DONE!'
  
  return
  end subroutine write_rst_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Save data using NEMO Xios infrastructure
!           This require remapping data from BFM to NEMO memory structure
!
! !INTERFACE:
  subroutine bfm_iom(kstp)
! 
! USES only 
   use iom, only: iom_put

   use mem, only: NO_BOXES,D3STATE,D3DIAGNOS,D3FLUX_FUNC,D2DIAGNOS,  &
            D2STATE_BEN,D2DIAGNOS_BEN,D2DIAGNOS_BEN,D2FLUX_FUNC_BEN
#if defined INCLUDE_SEAICE
   use mem, only: D2STATE_ICE,D2DIAGNOS_ICE,D2DIAGNOS_ICE,D2FLUX_FUNC_ICE
#endif
   use global_mem, only: RLEN, LOGUNIT, bfm_lwp
   use api_bfm, only: var_names, var_ids
   use api_bfm, only: SEAmask, BOTmask, SRFmask, c1dim, ZEROS
   use api_bfm, only: stStart, stEnd,                                &
#if defined INCLUDE_SEAICE
           & stIceStart, stIceEnd, stIceStateS, stIceStateE,         &
           & stIceDiag2dS, stIceDiag2dE, stIceFlux2dS, stIceFlux2dE, &
#endif
           & stBenStart, stBenEnd, stBenStateS, stBenStateE,         &
           & stBenDiag2dS, stBenDiag2dE, stBenFlux2dS, stBenFlux2dE, &
           &  stPelStart, stPelFluxE, stPelStateS, stPelStateE,      &
           & stPelDiagS, stPelDiagE, stPelFluxS, stPelFluxE,         &
           & stPelDiag2dS, stPelDiag2dE, stPelSurS, stPelSurE,       &
           & stPelBotS, stPelRivE

   IMPLICIT NONE

   INTEGER,intent(in)    :: kstp

   INTEGER               :: n, idx_tmp

   !---------------------------------------------
   ! Pelagic 3D variables
   !---------------------------------------------
   do n = stPelStart , stPelFluxE
      if ( var_ids(n) > 0 ) then  
         !-- Store pelagic state variables
         if ( n >= stPelStateS .AND. n <= stPelStateE ) then
            idx_tmp=n-stPelStateS+1
            call iom_put ( TRIM(var_names(n)) , unpack(D3STATE(idx_tmp,:),SEAmask,ZEROS) )
         end if
         !-- Store pelagic diagnostics
         if ( n >= stPelDiagS .AND. n <= stPelDiagE ) then
            idx_tmp=n-stPelDiagS+1
            call iom_put ( TRIM(var_names(n)) , unpack(D3DIAGNOS(idx_tmp,:),SEAmask,ZEROS) )
         end if
         !-- Store pelagic fluxes
         if ( n >= stPelFluxS .AND. n <= stPelFluxE ) then
            idx_tmp=n-stPelFluxS+1
            call correct_flux_output(1,idx_tmp,1,NO_BOXES,c1dim)
            call iom_put ( TRIM(var_names(n)) , unpack(c1dim,SEAmask,ZEROS) )
         endif
      endif
   enddo

   !---------------------------------------------
   ! Pelagic 2D variables
   !---------------------------------------------
   do n = stPelDiag2dS , stPelRivE
      if ( var_ids(n) > 0 ) then
         ! Store pelagic 2D diagnostics
         if ( n >= stPelDiag2dS .AND. n <= stPelDiag2dE ) then
            idx_tmp=n-stPelDiag2dS+1
            call iom_put ( TRIM(var_names(n)) , unpack(D2DIAGNOS(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
         ! Store pelagic 2D diagnostics at surface
         if ( n >= stPelSurS .AND. n <= stPelSurE) then
            idx_tmp=n-stPelDiag2dS+1
            call iom_put ( TRIM(var_names(n)) , unpack(D2DIAGNOS(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
         ! Store pelagic 2D diagnostics at bottom
         if ( n >= stPelBotS .AND. n <= stPelRivE) then
            idx_tmp=n-stPelDiag2dS+1
            call iom_put ( TRIM(var_names(n)) , unpack(D2DIAGNOS(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
            !iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2DIAGNOS(idx_tmp,:))
         end if
      end if
   end do

#if defined INCLUDE_SEAICE
   !---------------------------------------------
   ! 2D Seaice variables
   !---------------------------------------------
   do n = stIceStart , stIceEnd
      if ( var_ids(n) > 0 ) then
         ! Store seaice 2D state
         if ( n >= stIceStateS .AND. n <= stIceStateE) then
            idx_tmp=n-stIceStateS+1
            call iom_put( TRIM(var_names(n)) , unpack(D2STATE_ICE(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
         ! Store seaice 2D diagnostics
         if ( n >= stIceDiag2dS .AND. n <= stIceDiag2dE ) then
            idx_tmp=n-stIceDiag2dS+1
            call iom_put( TRIM(var_names(n)) , unpack(D2DIAGNOS_ICE(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
         ! Store seaice 2D flux
         if ( n >= stIceFlux2dS .AND. n <= stIceFlux2dE ) then
            idx_tmp=n-stIceFlux2dS+1
            call iom_put( TRIM(var_names(n)) , unpack(D2FLUX_FUNC_ICE(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
      end if
   end do
#endif

   !---------------------------------------------
   ! 2D Benthic variables
   !---------------------------------------------
   do n = stBenStart , stBenEnd
      if ( var_ids(n) > 0 ) then
         ! Store benthic 2D state
         if ( n >= stBenStateS .AND. n <= stBenStateE) then
            idx_tmp=n-stBenStateS+1
            call iom_put( TRIM(var_names(n)) , unpack(D2STATE_BEN(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
         ! Store benthic 2D diagnostics
         if ( n >= stBenDiag2dS .AND. n <= stBenDiag2dE ) then
            idx_tmp=n-stBenDiag2dS+1
            call iom_put( TRIM(var_names(n)) , unpack(D2DIAGNOS_BEN(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
         ! Store benthic 2D flux
         if ( n >= stBenFlux2dS .AND. n <= stBenFlux2dE ) then
            idx_tmp=n-stBenFlux2dS+1
            call iom_put( TRIM(var_names(n)) , unpack(D2FLUX_FUNC_BEN(idx_tmp,:),SRFmask(:,:,1),ZEROS(:,:,1)) )
         end if
      end if
   end do

  return
  end subroutine bfm_iom
!EOC

END MODULE trcdiabfm
