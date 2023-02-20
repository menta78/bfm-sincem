!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! MODULE: envforcing
!
! DESCRIPTION
!   This module contains all the ancillary functions for forcing functions.
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
!
! INCLUDE
#include "cppdefs.h"
!
! INTERFACE
   module envforcing
!
! USES
   use global_mem, only:RLEN,ONE,PI
   use constants,  only:SEC_PER_DAY
   use mem, only: iiC,iiN,iiP,iiS

   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public daylength, temperature, salinity, lightAtTime
   public instLight, light, wind, deposition, density
   public read_obs
! !PUBLIC DATA MEMBERS:
   !---------------------------------------------
   ! forcing function parameters
   !---------------------------------------------
   integer,public            :: forcing_method = 1
   !---------------------------------------------
   ! Parameters for analytical forcings
   ! Note: all read from namelist
   !---------------------------------------------
   integer,    public :: ltype
   real(RLEN), public :: tw,ts,tde,sw,ss,lw,ls,ww,ws,CO2inc
   real(RLEN), public :: botdep_c,botdep_n,botdep_p,botdep_si,botox_o
   !---------------------------------------------
   ! arrays for integration routines
   !---------------------------------------------
   logical,public              :: init_forcing_vars=.true.
   logical,public              :: use_external_data=.false.
   logical,public              :: use_event_data=.false.
   logical,public              :: use_seaice_data=.false.
   logical,public              :: use_benthic_data=.false.
   integer,parameter,public    :: unit_forcing=201
   integer,parameter,public    :: unit_data=202
   integer,parameter,public    :: unit_seaice=203
   integer,parameter,public    :: unit_event=204
   integer,parameter,public    :: unit_benthic=205
   character(len=128),public   :: forcing_file,seaice_file, &
                                  data_file, event_file, benthic_file
!
! !PRIVATE DATA MEMBERS:
   real(RLEN),parameter :: RFACTOR=PI/180._RLEN
!  pre-defined parameters
   integer, parameter,public   :: READ_SUCCESS=1
   integer, parameter,public   :: END_OF_FILE=-1
   integer, parameter,public   :: READ_ERROR=-2
   integer, parameter,public   :: NOTHING=0
!

   contains

!-----------------------------------------------------------------------

   FUNCTION daylength(time,latitude,ylength)
!
! DESCRIPTION
!   This function computes the length of the daylight period in hours
!   as a function of time of the year (fractional days) and latitude
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN),intent(in) :: time
   real(RLEN),intent(in) :: latitude
   real(RLEN),intent(in),optional :: ylength
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: daylength

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN)           :: declination,cycle
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   if (present(ylength)) then
      cycle = ylength
   else
      cycle = 360.0_RLEN
   end if
   declination = -0.406*cos(2.*PI*int(time)/cycle)
   daylength = acos(-tan(declination)*tan(latitude*RFACTOR))/PI*24.
   return

   END FUNCTION daylength

!-----------------------------------------------------------------------

   FUNCTION lightAtTime(df,dl)
!
! DESCRIPTION
!   This function determines whether there is light at a certain time
!   of the day. Returns an integer value 0 or 1
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN),intent(in) :: df,dl
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer :: lightAtTime

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: daytime, daylength
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   daytime=df*24. ! time of the day = fraction of the day * 24
   daytime=abs(daytime-12.) ! distance from noon
   daylength=dl/2.
   if(daytime.lt.daylength) then
     lightAtTime=1
   else
     lightAtTime=0
   endif
   return

   END FUNCTION lightAtTime

!-----------------------------------------------------------------------

   FUNCTION instLight(l,dl,df)
!
! DESCRIPTION
!   This function computes the instantaneous light at a certain time of the day
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN),intent(in) :: df,dl,l
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: instLight

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: daylength,daytime
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     daytime=df*24. ! time of the day = fraction of the day * 24
     daytime=abs(daytime-12.) ! distance from noon
     daylength=dl/2.
     if(daytime.lt.daylength) then
       daytime=daytime/daylength*PI
       instLight=l*cos(daytime)+l
     else
       instLight=0.
     endif
     return

   END FUNCTION instLight

!-----------------------------------------------------------------------

   FUNCTION salinity(dy,df)
!
! DESCRIPTION
!   This function provides an articial salinity value given the
!   parameters in the standalone.nml namelist
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: salinity
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     salinity=(ss+sw)/2.-(ss-sw)/2.*cos((dy+(df-.5))*RFACTOR)

   END FUNCTION salinity

!-----------------------------------------------------------------------

   FUNCTION temperature(dy,df)
!
! DESCRIPTION
!   This function provides an articial temperature value given the
!   parameters in the standalone.nml namelist
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: temperature
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     temperature=(ts+tw)/2.-(ts-tw)/2.*cos((dy+(df-.5))*RFACTOR) &
                    -tde*.5*cos(2*Pi*df)

   END FUNCTION temperature

!-----------------------------------------------------------------------

   FUNCTION density(tr,sr,dep)
!
! DESCRIPTION
!   This function computes density in kg/m3 from potential temperature
!   Mellor, 1991, J. Atmos. Oceanic Tech., 609-611
!
! USES
   use global_mem, only:RLEN
   use mem,        only: NO_BOXES

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN),intent(in),dimension(NO_BOXES) :: tr,sr,dep
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: density(NO_BOXES)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN),parameter :: GRAV=9.806E0_RLEN
   real(RLEN),dimension(NO_BOXES) :: TR2,TR3,TR4,TR5,SR2,P,P2,CR
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      TR2=TR  * TR
      TR3=TR2 * TR
      TR4=TR3 * TR
      TR5=TR4 * TR
      SR2=SR  * SR
      ! approximate pressure in units of bars
      P=-GRAV*1.025_RLEN*dep*0.01_RLEN
      p2=p*p
      CR = 1449.1_RLEN+0.0821_RLEN*P+4.55_RLEN*   &
         TR-0.045_RLEN*TR2+1.34_RLEN*(SR-35._RLEN)
      CR=P/(CR*CR)

      density = (999.842594_RLEN          + 6.793952E-2_RLEN*TR &
            -   9.095290E-3_RLEN*TR2   + 1.001685E-4_RLEN*TR3   &
            -   1.120083E-6_RLEN*TR4   + 6.536332E-9_RLEN*TR5   &
            +(  0.824493_RLEN          - 4.0899E-3_RLEN  *TR    &
            +   7.6438E-5_RLEN  *TR2   - 8.2467E-7_RLEN  *TR3   &
            +   5.3875E-9_RLEN  *TR4  )             *SR         &
            +( -5.72466E-3_RLEN        + 1.0227E-4_RLEN  *TR    &
            -   1.6546E-6_RLEN  *TR2  )             *SR**1.5    &
            +   4.8314E-4_RLEN  *SR2                   )        &
            +   1.E5_RLEN       *CR    *  (ONE-(CR+CR))

   END FUNCTION density

!-----------------------------------------------------------------------

   FUNCTION light(dy,df)
!
! DESCRIPTION
!   This function provides an articial light value given the
!   parameters in the standalone.nml namelist
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: light
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     light=(ls+lw)/2.-(ls-lw)/2.*cos(dy*RFACTOR)

   END FUNCTION light

!-----------------------------------------------------------------------

   FUNCTION wind(dy,df)
!
! DESCRIPTION
!   This function provides an artificial wind value given the
!   parameters in the standalone.nml namelist
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: wind
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     wind=(ws+ww)/2.0_RLEN-(ws-ww)/2.0_RLEN*cos((dy+(df-.5_RLEN))*RFACTOR)

   END FUNCTION wind

!-----------------------------------------------------------------------

   FUNCTION deposition(dy,df,botdep_ch,elem)
!
! DESCRIPTION
!   This function provides an articial deposition of detritus given the
!   parameters in the standalone.nml namelist
!   Uses a 3rd order spline with peaks in march and october
!   and minima in july and january 
!
! USES
   use global_mem, only:RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,intent(in)    :: dy,elem
   real(RLEN),intent(in) :: df,botdep_ch
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: deposition

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   real(RLEN) :: botdep_cl,min_factor
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     if (elem.eq.iiN) then
       min_factor=0.1
     else if (elem.eq.iiP) then
       min_factor=0.1
     else if (elem.eq.iiS) then
       min_factor=1
     else
       min_factor=1
     end if

     botdep_cl = 0.01_RLEN*botdep_ch*min_factor

     if (dy.lt.73) then
       deposition=(-0.00000362293541*(dy+df)**3  &
                   +0.00044863654442*(dy+df)**2  &
                   -0.00052603670261*(dy+df)     &
                   +0.057)
     else if (dy.lt.195) then
       deposition=(0.00000180691652*(dy+df-73)**3 &
                  -0.00034478631123*(dy+df-73)**2 &
                  +0.00705503032021*(dy+df-73)    &
                  +1.0)
     else if (dy.lt.287) then
        deposition=(-0.00000259580459*(dy+df-195)**3 &
                    +0.00031654513676*(dy+df-195)**2 &
                    +0.00360960703444*(dy+df-195)    &
                    +0.01*min_factor)
     else
        deposition=(0.00000381561750*(dy+df-287)**3 &
                   -0.00039989693012*(dy+df-287)**2 &
                   -0.00405875795472*(dy+df-287)    &
                   +1.0)
     end if
     deposition=deposition*botdep_ch

     if (deposition.lt.botdep_cl) then
       deposition=botdep_cl
     end if

   END FUNCTION deposition

!-----------------------------------------------------------------------

   subroutine read_obs(unit,yy,mm,dd,hh,min,ss,N,obs,ierr)
!
! DESCRIPTION
!   This routine will read all non-profile observations.
!   The routine allows for reading more than one scalar variable at a time.
!   The number of data to be read is specified by {\tt N}.
!   Data read-in are returned
!   in the 'obs' array. It is up to the calling routine to assign
!   meaningful variables to the individual elements in {\tt obs}.
!
! USES
   use global_mem, only: RLEN

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: N
  ! OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer, intent(out)                :: yy,mm,dd,hh,min,ss
   REAL(RLEN),intent(out)              :: obs(:)
   integer, intent(out)                :: ierr

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   character                 :: c1,c2,c3,c4
   integer                   :: i
   character(len=200)        :: cbuf
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   ierr=0
   read(unit,'(A200)',ERR=100,END=110) cbuf
   read(cbuf,900,ERR=100,END=110) yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   read(cbuf(20:),*,ERR=100,END=110) (obs(i),i=1,N)

   return
100 ierr=READ_ERROR
   return
110 ierr=END_OF_FILE
   return
900 format(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)

   end subroutine read_obs

   END MODULE envforcing

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

