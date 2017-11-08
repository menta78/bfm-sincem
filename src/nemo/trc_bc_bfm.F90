SUBROUTINE trc_bc_bfm ( kt, m )
   !!----------------------------------------------------------------------
   !!                  ***  ROUTINE trc_bc_bfm  ***
   !!                   
   !! ** Purpose :   Compute the tracer surface boundary condition trend of
   !!      (concentration/dilution effect) and add it to the general 
   !!       trend of tracer equations.
   !!
   !! ** Method :
   !!      * concentration/dilution effect:
   !!            The surface freshwater flux modify the ocean volume
   !!         and thus the concentration of a tracer as :
   !!            tra = tra + emp * trn / e3t   for k=1
   !!         where emp, the surface freshwater budget (evaporation minus
   !!         precipitation ) given in kg/m2/s is divided
   !!         by 1000 kg/m3 (density of plain water) to obtain m/s.
   !!
   !!         Runoff is treated separately according to the choice 
   !!         of online/offline coupling and river load
   !!
   !! ** Action  : - Update the 1st level of tra with the trend associated
   !!                with the tracer surface boundary condition 
   !!
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !! Adapted to usage with the BFM by M. Vichi (CMCC-INGV)
   !! Added runoff input of tracer concentrations
   !! Added surface input of tracer concentrations
   !!----------------------------------------------------------------------

   ! NEMO
   USE oce_trc              ! ocean dynamics and active tracers variables
   USE sbcrnf               ! contains river runoff (kg/m2/s) and river mask
   USE trc                  ! ocean  passive tracers variables
   USE fldread
   USE trcbc
   USE iom                  !  print control for debugging
   ! BFM
   use api_bfm
   use mem
   use global_mem,          only: LOGUNIT
   use sbc_oce,             only: ln_rnf
   use constants,           only: SEC_PER_DAY
   use print_functions
#ifdef INCLUDE_PELFE
   use mem_PelChem,         only: p_rDust
#endif
 
   ! substitutions
#  include "top_substitute.h90"

   !! * Arguments
   INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
   integer, intent(IN)     ::  m      ! BFM variable index

   !! * Local declarations
   INTEGER  ::   ji, jj, jk, jn             ! dummy loop indices
   REAL(wp) ::   ztra, zsrau, zse3t     ! temporary scalars
   TYPE(FLD), DIMENSION(1) ::   sf_dta  ! temporary array of information on the field to read
   REAL(wp), POINTER, DIMENSION(:,:  ) :: zsfx, field
   REAL(wp), POINTER, DIMENSION(:,:,:) :: IronDep
   REAL(wp) :: dustsink,  dustinp
   CHARACTER (len=25) :: charout
   !!---------------------------------------------------------------------
   !
   IF( nn_timing == 1 )  CALL timing_start('trc_bc_bfm')
   !
   ! Allocate temporary workspace
   CALL wrk_alloc( jpi, jpj,      zsfx  )
   CALL wrk_alloc( jpi, jpj, jpk, IronDep  )
#ifdef DEBUG
   CALL wrk_alloc( jpi, jpj,      field  )
   field = 0.0_wp
#endif
   !!----------------------------------------------------------------------
   ! initialization of density and scale factor
   zsrau = 1._wp / rau0
#ifdef INCLUDE_PELFE
   ! dust temporary array
   if ( m .eq. ppN7f ) IronDep = 0._wp
#endif

   ! Coupling online : 
   ! 1) constant volume: river runoff is added to the horizontal divergence (hdivn) in the subroutine sbc_rnf_div 
   !    one only consider the concentration/dilution effect due to evaporation minus precipitation + 
   !    freezing/melting of sea-ice
   ! 2) variable volume : all freshwater fluxes are added to the volume change (no additional dilution)
   ! Coupling offline : runoff are in emp which contains E-P-R
   ! If there is no input of biogeochemical variables associated to the river, it is assumed that
   ! there is no dilution associated with the river runoff (the freshwater has the same concentration of the sea)
   zsfx(:,:) = emp(:,:)                                    ! - standard case: on/offline coupling with const vol
   IF( .NOT. lk_offline .AND. lk_vvl )  zsfx(:,:) = 0._wp  ! - online coupling with variable volume
   IF( ln_rnf .AND. .NOT. ln_trc_cbc(m) )  &               ! - remove river dilution effect in the absence
      zsfx(:,:) = zsfx(:,:) + rnf(:,:)                     ! of a river load

   ! Concentration and dilution effect on tra 
   DO jj = 2, jpj
      DO ji = fs_2, fs_jpim1   ! vector opt.
           zse3t = 1. / fse3t(ji,jj,1)
           ztra = zsrau * zse3t * tmask(ji,jj,1) * zsfx(ji,jj) * trn(ji,jj,1,m) 
           ! add the trend to the general tracer trend
           tra(ji,jj,1,m) = tra(ji,jj,1,m) + ztra
      END DO
   END DO

    ! read and add surface input flux if needed
    IF (ln_trc_sbc(m)) THEN
       jn = n_trc_indsbc(m)
       DO jj = 2, jpj
          DO ji = fs_2, fs_jpim1   ! vector opt.
             zse3t = 1. / fse3t(ji,jj,1)
             ! The units in BFM input files are 1/day
             tra(ji,jj,1,m) = tra(ji,jj,1,m) + rf_trsfac(jn) * sf_trcsbc(jn)%fnow(ji,jj,1) &
                              * zse3t / SEC_PER_DAY
#ifdef INCLUDE_PELFE
             ! Store surface flux to output in IRONDEPO
             if ( m .eq. ppN7f ) &
             IronDep(ji,jj,1) = rf_trsfac(jn) * sf_trcsbc(jn)%fnow(ji,jj,1) * zse3t / SEC_PER_DAY
#endif
          END DO
       END DO

#ifdef INCLUDE_PELFE
       ! dust dissolution below the surface (Moore et al., 2008; Aumont et al.,2015)
       if ( m .eq. ppN7f .AND. p_rDust > 0. ) then
          ! dust sink length scale (1/m): dust assumed to be 0.01% (1/d) , sink speed p_rDust (m/d)
          dustsink =  0.0001_wp * 1. / p_rDust
          if (lwp) write(numout,*) " dust sink length scale : ", dustsink, p_rDust
          DO jj = 2, jpt
             DO ji = fs_2, fs_jpim1   ! vector opt.
                !  dust (g/day/m2) - > (umolFe/day/m2) 
                dustinp = 0.035 * 1.e06 / 55.85 * sf_trcsbc(jn)%fnow(ji,jj,1) 
                DO jk = 2, jpkm1
                   IronDep(ji,jj,jk) = dustinp * dustsink / SEC_PER_DAY * EXP( -fsdept(ji,jj,jk) / 540. )
                   tra(ji,jj,jk,m) = tra(ji,jj,jk,m) + IronDep(ji,jj,jk)
                END DO  
             END DO
          END DO
          IRONDEPO(:) = pack(IronDep,SEAmask)
       endif
#endif
    END IF

    ! Add mass from prescribed river concentration if river values are given
    ! An istantaneous mixing in the cell volume is assumed, the time unit of BFM input files must be 1/day
    IF (ln_rnf .AND. ln_trc_cbc(m)) THEN 
       jn = n_trc_indcbc(m)
        DO jj = 2, jpj
           DO ji = fs_2, fs_jpim1   ! vector opt.
              zse3t = 1. / (e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
              ztra = rn_rfact * rf_trcfac(jn) * sf_trccbc(jn)%fnow(ji,jj,1) * zse3t / SEC_PER_DAY
              tra(ji,jj,1,m) = tra(ji,jj,1,m) + ztra
#ifdef DEBUG
              field(ji,jj) = ztra
#endif
           END DO
        END DO
    END IF

#ifdef DEBUG
    charout = TRIM( var_names(stPelStateS+m-1) )
    WRITE(LOGUNIT,*) ''//charout//' trends in trc_sbc'
    WRITE(LOGUNIT,*)
    WRITE(LOGUNIT,*)'  level = 1'
    CALL prxy( LOGUNIT, 'trend at level = 1',field(:,:), jpi, 1, jpj, 1, ZERO)
    CALL prxy( LOGUNIT, 'trn at level = 1',trn(:,:,1,m), jpi, 1, jpj, 1, ZERO)
    CALL wrk_dealloc( jpi, jpj,      field  )       
#endif

   CALL wrk_dealloc( jpi, jpj,      zsfx  )       
   CALL wrk_dealloc( jpi, jpj, jpk, IronDep  )
    IF( nn_timing == 1 )  CALL timing_stop('trc_bc_bfm')

   END SUBROUTINE trc_bc_bfm

