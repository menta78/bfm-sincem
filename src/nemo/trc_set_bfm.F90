#include "domzgr_substitute.h90"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: trc_set_bfm.F90
!
! !DESCRIPTION:
!  Computes additional boundary conditions and transfer sinking velocity
!
! !INTERFACE:
   subroutine trc_set_bfm(kt,m)
!
! !USES:
   ! NEMO
   use oce_trc          ! ocean dynamics and active tracers variables
   use trc              ! ocean passive tracers variables
   ! BFM
   use global_mem, only:RLEN,ZERO,ONE
   use mem_param,  only: AssignAirPelFluxesInBFMFlag,        &
                         AssignPelBenFluxesInBFMFlag
   use mem_PelGlobal, only: p_rR6m, KSINK_rPPY,              &
                            AggregateSink, depth_factor,     &
                            SINKD3STATE 
   use mem
   use mem_Settling
   use constants,    only: SEC_PER_DAY
   use api_bfm
!
!
! !AUTHORS
!   Marcello Vichi (CMCC-INGV)
!
! !REVISION_HISTORY
!  2015: Phyto sinking formulation as in PISCES model.
!  2016: Implement control based on sinking velocites (T. Lovato)
!
! COPYING
!
!   Copyright (C) 2016 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Implicit typing is never allowed
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   implicit none
   !
   ! !INPUT/OUTPUT PARAMETERS:
   integer, intent(IN)     ::  kt  ! ocean time-step index
   integer, intent(IN)     ::  m   ! BFM variable index
   !
   ! !LOCAL VARIABLES:
   ! 3D sinking velocity field
   integer               :: ji, jj, jk, n, iistate
   real(RLEN)            :: zfact,timestep,wsmax
   real(RLEN)            :: wbio(jpi,jpj,jpk)   
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   !  Biological timestep (in days)
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   timestep  = rdt*FLOAT(nn_dttrc)/SEC_PER_DAY

   ! Transfer sinking velocities 
   ! (negative, z-axis is positive upwards)
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   if ( SINKD3STATE(m)%dosink ) then 

      ! Phytoplankton group 
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      if ( SINKD3STATE(m)%group == 1 ) then

         ! Prescribe sinking velocity below depth threshold KSINK_rPPY
         ! (usually below 150m). This accelerate diatoms sinking to balance
         ! the reduction occurring in deeper layers where nutrient concen-
         ! tration is high. It is alternatevly done by AggregationSink
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         if ( KSINK_rPPY > 0 .AND. .NOT. AggregateSink ) &
            where( EPR > KSINK_rPPY ) SINKD3STATE(m)%sedi = p_rR6m

         ! Prescribe burial velocity at the bottom
         SINKD3STATE(m)%sedi(BOTindices) = p_burvel_PI
      endif 

      wbio = -unpack(SINKD3STATE(m)%sedi,SEAmask,ZEROS)

      ! Sinking speeds increase with depth below the turbocline depth due
      ! to aggregation. Velocity is limited according to the depth of the 
      ! layer (80% of it). This modulates the assigned sinking rate
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      if ( AggregateSink ) then
         do jk=1,jpk-1
            do jj=1,jpj
               do ji=1,jpi
                  wsmax=0.8_RLEN*fse3t(ji,jj,jk)/timestep
                  zfact = max(ZERO,exp((fsdepw(ji,jj,jk)-hmld(ji,jj)) &
                                           /depth_factor)-ONE)
                  wbio(ji,jj,jk) = min(wsmax,(ONE+zfact)*wbio(ji,jj,jk))

                end do
            end do
         end do
      endif

      ! Compute vertical sinking with upwind scheme
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      CALL trc_sink_bfm(wbio,m)
   
   endif

   return
   end subroutine trc_set_bfm

!EOC

