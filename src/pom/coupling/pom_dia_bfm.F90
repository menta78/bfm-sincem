#include "INCLUDE.h"
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pom_dia_bfm.F90
!
! !INTERFACE:
   subroutine pom_dia_bfm(ltimed,kt,TT)
!
! !DESCRIPTION:
!
! !USES:
   
   use global_mem, only:RLEN
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use netcdf_bfm, only: save_rst_bfm, ncid_rst
   use mem
   use api_bfm, only: out_delta,D3STATEB
   use constants,  only:SEC_PER_DAY
   use Service, only: nitend,deltat

   implicit none
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in)    :: ltimed !time in days
   integer, intent(in)      :: kt
   real(RLEN)               :: TT
! !LOCAL VARIABLES:
!time in seconds
   real(RLEN)               :: localtime 
   integer                  :: time_to_save
!
!EOP
!-----------------------------------------------------------------------
!BOC
   localtime = ltimed*SEC_PER_DAY      !total seconds
   time_to_save=nint(out_delta*SEC_PER_DAY/deltat)
    if(kt.eq.1) then

    write(6,*) 'IN POM_DIA', nitend,kt, time_to_save,localtime,ltimed, TT, deltat,out_delta
    endif
   !---------------------------------------------
   ! Update means
   !---------------------------------------------
   call calcmean_bfm(ACCUMULATE)

   !---------------------------------------------
   ! Write diagnostic output
   !---------------------------------------------
   !   if ( lwp .AND. MOD( kt, nwritetrc ) == 0 ) then
!   if ( MOD( kt, time_to_save ) == 0 ) then
!     write(6,*) 'TT = ' , TT
!     write(6,*) 'out_delta2 = ', out_delta
!     if (TT .ge. out_delta)  then
   if(TT+(deltat/SEC_PER_DAY).gt.out_delta) then
!      write(6,*) 'pom_dia_bfm : write NetCDF passive tracer concentrations at ', kt, 'time-step'
!      write(6,*) '~~~~~~~~~~~~ '
      call flush(6)
      call calcmean_bfm(MEAN)
      call save_bfm(localtime)
!      TT = TT - int(TT)
       TT = TT - nint(TT)
   end if
   !---------------------------------------------
   ! Close the file (MAV: are we saving the last step?)
   !---------------------------------------------
    if ( kt >= nitend ) then
       if(-TT.le.deltat/SEC_PER_DAY) then
!      call save_rst_bfm(localtime)        !     (G)
      write(6,*) 'D3STATE DOPO WRST ', D3STATE(ppP2c,1), D3STATE(ppP2c,NO_BOXES)
      write(6,*) 'D3STATEB DOPO WRST ', D3STATEB(ppP2c,1), D3STATEB(ppP2c,NO_BOXES)
      call flush(6)
!      call close_ncdf(ncid_rst)
      call close_ncdf(ncid_bfm)
      write (6,*) 'POM_DIA: NETCDF RESTART WRITTEN, TIME--> ', ltimed,kt, nitend
     endif
   endif

!   !---------------------------------------------
!   ! Reset the arrays
!   !---------------------------------------------
!   call ResetFluxes
!
   return
   end subroutine pom_dia_bfm

!EOC


