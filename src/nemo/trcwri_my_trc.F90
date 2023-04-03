MODULE trcwri_my_trc
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!     trc_wri_my_trc   :  outputs of concentration fields
   !!======================================================================
#include "cppdefs.h"
#if defined key_top && defined key_xios
   !!----------------------------------------------------------------------
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   ! NEMO
   USE par_trc         ! passive tracers common variables
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   ! BFM
   USE global_mem, ONLY: LOGUNIT, bfm_lwp
   USE time,       ONLY: bfmtime
   USE api_bfm,    ONLY: stStart, stEnd, var_names,                 &
#if defined INCLUDE_SEAICE
          & stIceStart, stIceEnd, stIceStateS, stIceStateE,         &
          & stIceDiag2dS, stIceDiag2dE, stIceFlux2dS, stIceFlux2dE, &
#endif
          & stBenStart, stBenEnd, stBenStateS, stBenStateE,         &
          & stBenDiag2dS, stBenDiag2dE, stBenFlux2dS, stBenFlux2dE, &
          &  stPelStart, stPelFluxE, stPelStateS, stPelStateE,      &
          & stPelDiagS, stPelDiagE, stPelFluxS, stPelFluxE,         &
          & stPelDiag2dS, stPelDiag2dE, stPelSurS, stPelSurE,       &
          & stPelBotS, stPelBotE
#ifdef INCLUDE_PELCO2
   use mem_CO2, ONLY: CloseCO2
#endif

   IMPLICIT NONE
   PRIVATE

   INTEGER, ALLOCATABLE, DIMENSION(:) :: id_dia2d ! BFM indexes of 2D diags
   INTEGER, ALLOCATABLE, DIMENSION(:) :: id_dia3d ! BFM indexes of 3D diags

   PUBLIC trc_wri_my_trc, diags_mapping, diags_collect

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_my_trc.F90 14239 2020-12-23 08:57:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_wri_my_trc( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: Kmm   ! time level indices
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      !-------------------------------------------------------
      ! Pelagic
      DO jn = 1, jp_bgc
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, tr(:,:,:,jn,Kmm) )
      END DO
      ! Benthic
      DO jn = 1, jp_bgc_b
         cltra = TRIM( ctrcnm_b(jn) )                  ! short title for tracer
         CALL iom_put( cltra, tr_b(:,:,:,jn) )
      END DO
#ifdef INCLUDE_SEAICE
      ! Seaice
      DO jn = 1, jp_bgc_i
         cltra = TRIM( ctrcnm_i(jn) )                  ! short title for tracer
         CALL iom_put( cltra, tr_i(:,:,:,jn) )
      END DO
#endif
      ! Write diagnostics
      !---------------------------------------------
      ! 3-D
      DO jn = 1 , jp_dia3d
         cltra = TRIM( var_names(id_dia3d(jn)) )
         CALL iom_put( cltra, trc3d(:,:,:,jn) )
      ENDDO
      ! 2-D
      DO jn = 1 , jp_dia2d
         cltra = TRIM( var_names(id_dia2d(jn)) )
         CALL iom_put( cltra, trc2d(:,:,jn) )
      ENDDO

      ! Clear BFM memory
      !---------------------------------------------
      IF ( bfmtime%stepnow == bfmtime%stepEnd ) THEN
         !close systemforcings
#ifdef INCLUDE_PELCO2
         call CloseCO2()
#endif
         ! clear main memory
         call ClearMem

         LEVEL1 ' '
         LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
         LEVEL1 '             EXPERIMENT FINISHED               '
         LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
         LEVEL1 ' '
      ENDIF
      !
   END SUBROUTINE trc_wri_my_trc


   SUBROUTINE diags_collect(kt, ji, jj, bot)
      !!----------------------------------------------------------------------
      !!                     ***  diags_mapping  ***
      !!
      !! ** Purpose : Fill 2D and 3D diagnostics with data from BFM 
      !!
      !! ** Method  : 
      !!----------------------------------------------------------------------
      USE api_bfm, ONLY: c1dim
      USE mem, ONLY: NO_BOXES,D3STATE,D3DIAGNOS,D3FLUX_FUNC,D2DIAGNOS,  &
            D2STATE_BEN,D2DIAGNOS_BEN,D2DIAGNOS_BEN,D2FLUX_FUNC_BEN
#if defined INCLUDE_SEAICE
      USE mem, ONLY: D2STATE_ICE,D2DIAGNOS_ICE,D2DIAGNOS_ICE,D2FLUX_FUNC_ICE
#endif
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   ji, jj   ! dummy loop index
      INTEGER, INTENT(in) ::   bot  ! bottom level id
      INTEGER :: jn, jl, idx
      !!----------------------------------------------------------------------

      ! 3D diagnostics
      !---------------------------------------------
      DO jn = 1 , jp_dia3d
         jl = id_dia3d(jn)
         IF ( jl >= stPelDiagS .AND. jl <= stPelDiagE ) THEN
            idx = jl - stPelDiagS + 1 
            trc3d(ji, jj, 1:bot, jn) = D3DIAGNOS(1:bot,idx)
         ENDIF
         IF ( jl >= stPelFluxS .AND. jl <= stPelFluxE ) then
            idx = jl - stPelFluxS + 1
            call correct_flux_output(1,idx,1,NO_BOXES,c1dim)
            trc3d(ji, jj, 1:bot, jn) = c1dim(1:bot) 
         endif
      ENDDO

      ! 2D diagnostics
      !---------------------------------------------
      DO jn = 1 , jp_dia2d
         jl = id_dia2d(jn)
         IF ( jl >= stPelDiag2dS .AND. jl <= stPelBotE ) THEN
            idx = jl - stPelDiag2dS + 1
            trc2d(ji, jj, jn) = D2DIAGNOS(1,idx)
         ENDIF
         IF ( jl >= stBenDiag2dS .AND. jl <= stBenDiag2dE ) THEN
            idx = jl - stBenDiag2dS + 1
            trc2d(ji, jj, jn) = D2DIAGNOS_BEN(1,idx)
         ENDIF
         IF ( jl >= stBenFlux2dS .AND. jl <= stBenFlux2dE ) THEN
            idx = jl - stBenFlux2dS + 1
            trc2d(ji, jj, jn) = D2FLUX_FUNC_BEN(1,idx)
         ENDIF
#if defined INCLUDE_SEAICE
         IF ( jl >= stIceDiag2dS .AND. jl <= stIceDiag2dE ) THEN
            idx = jl - stIceDiag2dS + 1
            trc2d(ji, jj, jn) = D2DIAGNOS_ICE(1,idx)
         ENDIF
         If ( jl >= stIceFlux2dS .AND. jl <= stIceFlux2dE ) THEN
            idx = jl - stIceFlux2dS + 1
            trc2d(ji, jj, jn) = D2FLUX_FUNC_ICE(1,idx)
         ENDIF
#endif
      ENDDO

   END SUBROUTINE diags_collect


   SUBROUTINE diags_mapping()
      !!----------------------------------------------------------------------
      !!                     ***  diags_mapping  ***
      !!
      !! ** Purpose : Map 2D and 3D diagnostics based on XIOS I/O list
      !!
      !! ** Method  : Loop over BFM diagnostics and check if are needed for xios
      !!----------------------------------------------------------------------
      USE IOM,     ONLY: iom_use
      !
      USE global_mem, ONLY: ZERO
      !
      INTEGER, ALLOCATABLE, DIMENSION(:) :: id_2d, id_3d
      LOGICAL, DIMENSION(stEnd, 2) :: pp_index
      INTEGER :: jn, j_2d, j_3d
      !!----------------------------------------------------------------------

      LEVEL1 'Collect 2D & 3D BFM diagnostics'

      ! Inquire XIOS output list
      !-------------------------------------------------------
      pp_index(:,:) = .FALSE.
      DO jn = 1, stEnd
         IF (iom_use(TRIM(var_names(jn)))) THEN
            IF ( (jn >= stPelDiag2dS .AND. jn <= stPelBotE) .OR. &
#if defined INCLUDE_SEAICE
               (jn >= stIceDiag2dS .AND. jn <= stIceFlux2dE) .OR. &
#endif
               (jn >= stBenDiag2dS .AND. jn <= stBenFlux2dE) ) THEN
               pp_index(jn,2) = .TRUE.
            ELSEIF (jn >= stPelDiagS .AND. jn <= stPelFluxE) THEN
               pp_index(jn,1) = .TRUE.
            END IF
         ENDIF
      END DO

      ! set diagnostics indexes
      !-------------------------------------------------------
      jp_dia3d = COUNT(pp_index(:,1))
      jp_dia2d = COUNT(pp_index(:,2))
      ALLOCATE(id_dia3d(jp_dia3d), id_dia2d(jp_dia2d))
      ALLOCATE(trc3d(jpi,jpj,jpk,jp_dia3d), trc2d(jpi,jpj,jp_dia2d))
      trc3d = ZERO
      trc2d = ZERO
      !
      j_2d = 1
      j_3d = 1
      !
      DO jn = 1, stEnd
         IF (pp_index(jn,1) .eqv. .TRUE.) THEN
            id_dia3d(j_3d) = jn
            j_3d = j_3d + 1
         ENDIF
         IF (pp_index(jn,2) .eqv. .TRUE.) THEN
            id_dia2d(j_2d) = jn
            j_2d = j_2d + 1
         ENDIF
      ENDDO

      LEVEL1 ''

   END SUBROUTINE diags_mapping

#else

CONTAINS

   SUBROUTINE trc_wri_my_trc
      !
   END SUBROUTINE trc_wri_my_trc

#endif

END MODULE trcwri_my_trc
