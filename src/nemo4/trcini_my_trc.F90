MODULE trcini_my_trc
   !!======================================================================
   !!                         ***  MODULE trcini_my_trc  ***
   !! TOP :   initialisation of the MY_TRC tracers
   !!======================================================================
   !! History :        !  2007  (C. Ethe, G. Madec) Original code
   !!                  !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_ini_my_trc   : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_my_trc
   USE trcnam_my_trc     ! MY_TRC SMS namelist
   USE trcsms_my_trc
   !! BFM
   USE mem,        ONLY: D3STATE
   

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_my_trc   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_my_trc.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_my_trc( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_my_trc  ***  
      !!
      !! ** Purpose :   initialization for MY_TRC model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER :: jn
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_ini_my_trc : initialize BFM tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      ! dummy test
      DO jn = jp_myt0, jp_myt1
          tr(:,:,:,jn,Kmm) = D3STATE(jn,1)
      ENDDO

      !
   END SUBROUTINE trc_ini_my_trc

   !!======================================================================
END MODULE trcini_my_trc
