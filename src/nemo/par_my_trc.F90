MODULE par_my_trc
   !!======================================================================
   !!                        ***  par_my_trc  ***
   !! TOP :   set the MY_TRC parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_my_trc.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
   USE par_kind, ONLY: wp, lca

   IMPLICIT NONE

   ! Starting/ending MY_TRC do-loop indices (N.B. no MY_TRC : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC ::   jp_myt0             !: First index of MY_TRC passive tracers
   INTEGER, PUBLIC ::   jp_myt1             !: Last  index of MY_TRC passive tracers
   !
   INTEGER, ALLOCATABLE, DIMENSION(:), PUBLIC :: var_map
#if defined key_xios
   LOGICAL,PUBLIC    :: bfm_iomput=.TRUE. ! use xios in nemo
#else
   LOGICAL,PUBLIC    :: bfm_iomput=.FALSE.
#endif
   !!
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  bottom_level          !: deepest level of water column

   !! boundary forcings
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: fesed   ! Seabed supply of iron

   !! support arrays
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: chl_a
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ph
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: sink_rates

   !! benthic passive tracers  (input and output)
   !! ------------------------------------------
   INTEGER, PUBLIC ::   jp_bgc_b             !: Number of tracers
   INTEGER, PUBLIC ::   jpk_b                !: Number of vertical levels
   CHARACTER(len=lca), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ctrcnm_b   !: tracers short name
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  tr_b          !: tracers concentration

   !! seaice passive tracers  (input and output)
   !! ------------------------------------------
   INTEGER, PUBLIC ::   jp_bgc_i             !: Number of tracers
   INTEGER, PUBLIC ::   jpk_i                !: Number of vertical levels
   CHARACTER(len=lca), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ctrcnm_i   !: tracers short name
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  tr_i          !: tracers concentration
   !!======================================================================
END MODULE par_my_trc
