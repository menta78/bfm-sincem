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
   USE api_bfm, ONLY: stStart, stEnd, var_names,                    &
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
   PRIVATE

   INTEGER, ALLOCATABLE, DIMENSION(:) :: id_dia2d ! BFM indexes of 2D diags
   INTEGER, ALLOCATABLE, DIMENSION(:) :: id_dia3d ! BFM indexes of 3D diags

   PUBLIC trc_wri_my_trc, mapping_diags

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
      ! ---------------------------------------
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
      ! Diagnostics
      ! 3-D
      ! 2-D


      !
   END SUBROUTINE trc_wri_my_trc


   SUBROUTINE mapping_diags()
      !!---------------------------------------------------------------------
      ! Loop over all BFM variables and check if are needed for xios
      ! Allocate 2D and 3D diagnostics arrays accordingly
      !!---------------------------------------------------------------------
      USE IOM,     ONLY: iom_use
      !
      INTEGER, ALLOCATABLE, DIMENSION(:) :: id_2d, id_3d
      LOGICAL, DIMENSION(stEnd, 2) :: pp_index
      INTEGER :: jn, j_2d, j_3d
      !
      LEVEL1 'Collect 2D & 3D BFM diagnostics'
      !
      pp_index(:,:) = .FALSE.
      DO jn = 1, stEnd
         IF (iom_use(TRIM(var_names(jn)))) THEN
            IF ( (jn >= stPelDiag2dS .AND. jn <= stPelRivE) .OR. &
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
      !
      ! set diagnostics indexes
      jp_dia3d = COUNT(pp_index(:,1))
      jp_dia2d = COUNT(pp_index(:,2))
      ALLOCATE(id_dia3d(jp_dia3d), id_dia2d(jp_dia2d))
      ALLOCATE(trc3d(jpi,jpj,jpk,jp_dia3d), trc2d(jpi,jpj,jp_dia2d))
      !
      j_2d = 1
      j_3d = 1
      !
      DO jn = 1, stEnd
         IF (pp_index(jn,1) == .TRUE.) THEN
            id_dia3d(j_3d) = jn
            j_3d = j_3d + 1
         ENDIF
         IF (pp_index(jn,2) == .TRUE.) THEN
            id_dia2d(j_2d) = jn
            j_2d = j_2d + 1
         ENDIF
      ENDDO
      LEVEL1 'dia3d tot: ', jp_dia3d
      LEVEL1 'ids 3d:', id_dia3d(:)
      LEVEL1 'dia2d tot: ', jp_dia2d
      LEVEL1 'ids 2d:', id_dia2d(:)

   END SUBROUTINE mapping_diags

#else

CONTAINS

   SUBROUTINE trc_wri_my_trc
      !
   END SUBROUTINE trc_wri_my_trc

#endif

END MODULE trcwri_my_trc
