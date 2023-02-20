!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: external_benthic
!
! DESCRIPTION:
!   Read benthic data, interpolate in time
!   Coupling of benthic processes 
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#if defined BENTHIC_BIO || defined BENTHIC_FULL
!
! INCLUDE
#include "cppdefs.h"
!
! INTERFACE
   subroutine external_benthic
!
! USES
   use global_mem, only: RLEN,ZERO,ONE
   use constants,  only: SEC_PER_DAY
   ! benthic forcings
   use mem,        only: ETW_Ben,ESW_Ben,ERHO_Ben, Depth
   use time,       only: julianday, secondsofday, time_diff, &
                         julian_day,calendar_date
   use envforcing, only: init_forcing_vars, density, &
                         unit_benthic, read_obs, END_OF_FILE, READ_ERROR
   use bfm_error_msg

   IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,parameter                  :: NOBS_Ben=2
   integer                            :: yy,mm,dd,hh,minutes,ss
   real(RLEN)                         :: t,alpha,jday
   real(RLEN), save                   :: dt
   integer, save                      :: data_jul1,data_secs1
   integer, save                      :: data_jul2=0,data_secs2=0
   real(RLEN), save                   :: obs1(NOBS_Ben),obs2(NOBS_Ben)=0.
   integer                            :: ierr,jh,jn
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#ifdef DEBUG
   LEVEL1 'external_forcing (jul,sec): ',julianday,secondsofday
   call  calendar_date(real(julianday,RLEN),yy,mm,dd,jh,jn)
   LEVEL2 'Calendar day:',yy,mm,dd
#endif

   if (init_forcing_vars) then
     data_jul2=0
     data_secs2=0
     obs2(:)=ZERO
     ! check consistency of initial date
     call read_obs(unit_benthic,yy,mm,dd,hh,minutes,ss,NOBS_Ben,obs2,ierr)
     call julian_day(yy,mm,dd,0,0,jday)
     if (jday > julianday )  &
        call bfm_error('external_benthic','Model start date is earlier than start date of sea ice data')
     rewind(unit_benthic)
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .lt. 0) then
      do
         data_jul1 = data_jul2
         data_secs1 = data_secs2
         obs1 = obs2
         call read_obs(unit_benthic,yy,mm,dd,hh,minutes,ss,NOBS_Ben,obs2,ierr)
         select case (ierr)
           case (READ_ERROR)
              call bfm_error('external_benthic','Error reading forcing data')
           case (END_OF_FILE)
              call bfm_error('external_benthic','Model end date is beyond forcing data')
         end select
         call julian_day(yy,mm,dd,0,0,jday)
         data_jul2 = int(jday)
         data_secs2 = hh*3600 + minutes*60 + ss
         if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .gt. 0) EXIT
      end do
      dt = time_diff(data_jul2,data_secs2,data_jul1,data_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(julianday,secondsofday,data_jul1,data_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   ETW_Ben(:) = obs1(1) + t*alpha
   alpha = (obs2(2)-obs1(2))/dt
   ESW_Ben(:) = obs1(2) + t*alpha
!   alpha = (obs2(3)-obs1(3))/dt
   ERHO_Ben(:) = density(ETW_Ben(:),ESW_Ben(:),Depth(:))

#ifdef DEBUG
   LEVEL2 'ETW_Ben=',ETW_Ben
   LEVEL2 'ESW_Ben=',ESW_Ben
   LEVEL2 'ERHO_Ben=',ERHO_Ben
#endif
  return

   end subroutine external_benthic

#endif

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

