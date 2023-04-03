!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! MODULE: api_bfm
!
! DESCRIPTION:
!   API for the BFM.
!   Storage of variables and diagnostics
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
#include"cppdefs.h"
!
! INTERFACE
  module api_bfm
!
! USES
   use global_mem, only:RLEN,ZERO,bfm_lwp,LOGUNIT,NMLUNIT,bfm_file_FirstUnit
   use mem,        only:NO_D3_BOX_STATES

   implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bfm, find, update_save_delta, GetLun
!
! !PUBLIC DATA MEMBERS:
   logical                            :: bio_calc,bioshade_feedback,bfm_rstctl
   integer                            :: bio_setup  =1
   integer                            :: bfm_init  =0
   integer                            :: surface_flux_method=-1
   integer                            :: bottom_flux_method=-1
   integer                            :: n_surface_fluxes=-1
   character(len=PATH_MAX)            :: out_dir,out_fname,out_title
   integer                            :: out_units
   integer                            :: out_delta,out_secs,save_delta
   real(RLEN)                         :: time_delta
   character(len=PATH_MAX)            :: in_rst_fname, out_rst_fname
   logical                            :: unpad_out

   !---------------------------------------------
   ! parameters for massive parallel computation
   ! the following are the default values if the 
   ! macro BFM_PARALLEL is not defined
   !---------------------------------------------
   logical                            :: parallel = .FALSE.
   logical                            :: parallel_log = .FALSE.
   integer                            :: parallel_rank = 0
   character(LEN=4)                   :: str

   !---------------------------------------------
   ! Dimension lengths for output
   !---------------------------------------------
   integer        :: lon_len
   integer        :: lat_len
   integer        :: depth_len
   integer        :: ocepoint_len,surfpoint_len,botpoint_len

   !---------------------------------------------
   ! BFM variable information for output
   !---------------------------------------------
   integer, dimension(:), allocatable      :: var_ids
   logical, dimension(:), allocatable      :: var_ave
   real(RLEN)                              :: ave_count
   logical                                 :: ave_ctl = .false.
   real(RLEN),allocatable,dimension(:,:)   :: D3ave
   real(RLEN),allocatable,dimension(:,:)   :: D2ave
#if defined INCLUDE_SEAICE
   real(RLEN),allocatable,dimension(:,:)   :: D2ave_ice
#endif
   real(RLEN),allocatable,dimension(:,:)   :: D2ave_ben
   character(len=64), dimension(:), allocatable :: var_names
   character(len=64), dimension(:), allocatable :: var_units
   character(len=64), dimension(:), allocatable :: var_long
   !---------------------------------------------
   ! Indices of the various output variables
   !---------------------------------------------

   integer,public                            :: stStart=0
   integer,public                            :: stEnd=0

   integer,public                            :: stPelStateS=0
   integer,public                            :: stPelDiagS=0
   integer,public                            :: stPelFluxS=0

   integer,public                            :: stPelDiag2dS=0

   integer,public                            :: stPelSurS=0
   integer,public                            :: stPelBotS=0

   integer,public                            :: stPelStateE=0
   integer,public                            :: stPelDiagE=0
   integer,public                            :: stPelFluxE=0

   integer,public                            :: stPelDiag2dE=0

   integer,public                            :: stPelSurE=0
   integer,public                            :: stPelBotE=0

   integer,public                            :: stPelStart=0
   integer,public                            :: stPelEnd=0    

#if defined INCLUDE_SEAICE
   integer,public                            :: stIceStateS=0
   integer,public                            :: stIceDiag2dS=0
   integer,public                            :: stIceFlux2dS=0

   integer,public                            :: stIceStateE=0
   integer,public                            :: stIceDiag2dE=0
   integer,public                            :: stIceFlux2dE=0

   integer,public                            :: stIceStart=0
   integer,public                            :: stIceEnd=0    

#endif

   integer,public                            :: stBenStateS=0
   integer,public                            :: stBenDiag2dS=0
   integer,public                            :: stBenFlux2dS=0

   integer,public                            :: stBenStateE=0
   integer,public                            :: stBenDiag2dE=0
   integer,public                            :: stBenFlux2dE=0

   integer,public                            :: stBenStart=0
   integer,public                            :: stBenEnd=0    


   !---------------------------------------------
   ! Additional output variables
   !---------------------------------------------
   real(RLEN), dimension(:), allocatable   :: c1dim

   !---------------------------------------------
   ! BFM variable information for data input
   ! integer init: select the initialization
   !               0 = homogeneous
   !               1 = analytical
   !               2 = from file
   ! options for init==1
   ! real anv1: value in the surface layer
   ! real anz1: depth of the surface layer
   ! real anv2: value in the bottom layer
   ! real anz2: depth of the bottom layer
   ! options for init==2
   ! char filename: name of the input file
   ! char  varname: name of the var in input file
   ! Options currently used when coupled with NEMO
   ! logical obc: variable has open boundary data
   ! logical sbc: variable has surface boundary data
   ! logical cbc: variable has coastal boundary data
   ! logical rho: initial condition scaled by seawater density
   !---------------------------------------------
   type InputInfo
      integer           :: init
      real(RLEN)        :: unif
      character(LEN=40) :: filename
      character(LEN=40) :: varname
      real(RLEN)        :: anz1
      real(RLEN)        :: anv1
      real(RLEN)        :: anz2
      real(RLEN)        :: anv2
      logical           :: obc
      logical           :: sbc
      logical           :: cbc
      logical           :: rho
   end type InputInfo
   type(InputInfo),dimension(NO_D3_BOX_STATES) :: InitVar

   !---------------------------------------------
   ! Additional 1D arrays
   !---------------------------------------------
   ! progressive indices of bottom and surface points
   integer,allocatable,dimension(:),public     :: BOTindices,SRFindices
   ! real mask of river points at surface
   real(RLEN),allocatable,dimension(:),public  :: RIVmask
   ! Total amount for each variable
   real(RLEN),allocatable,dimension(:),public  :: D3STATE_tot
#ifdef INCLUDE_SEAICE
   real(RLEN),allocatable,dimension(:),public  :: D2STATE_ICE_tot
#endif
   real(RLEN),allocatable,dimension(:),public  :: D2STATE_BEN_tot


#if defined BFM_NEMO || defined BFM_POM
   !---------------------------------------------
   ! Additional 1D arrays
   !---------------------------------------------
   ! absolute indices ocean, surface and bottom points
   integer,allocatable,dimension(:),public     :: ocepoint
   integer,allocatable,dimension(:),public     :: surfpoint, botpoint

   !---------------------------------------------
   ! Additional 3D arrays
   !---------------------------------------------
   real(RLEN),allocatable,dimension(:,:,:),public  :: ZEROS
   ! 3D boolean Land-sea mask
   logical,allocatable,dimension(:,:,:),public     :: SEAmask
   ! 3D boolean sea-bottom mask
   logical,allocatable,dimension(:,:,:),public     :: BOTmask
   ! 3D boolean mask of the surface points
   logical,allocatable,dimension(:,:,:),public     :: SRFmask

   !---------------------------------------------
   ! Additional integration arrays
   ! for leapfrog scheme
   !---------------------------------------------
   real(RLEN),allocatable,dimension(:,:),public  :: D3STATEB
#if defined INCLUDE_SEAICE
   real(RLEN),allocatable,dimension(:,:),public  :: D2STATEB_ICE
#endif
   real(RLEN),allocatable,dimension(:,:),public  :: D2STATEB_BEN

   !---------------------------------------------
   ! Additional allocatable temporary arrays
   !---------------------------------------------
   logical,allocatable,dimension(:),public        :: btmp1D
   logical,allocatable,dimension(:,:),public      :: btmp2D
   logical,allocatable,dimension(:,:,:),public    :: btmp3D
   integer,allocatable,dimension(:),public        :: itmp1D
   integer,allocatable,dimension(:,:),public      :: itmp2D
   integer,allocatable,dimension(:,:,:),public    :: itmp3D
   real(RLEN),allocatable,dimension(:),public     :: rtmp1D
   real(RLEN),allocatable,dimension(:,:),public   :: rtmp2D
   real(RLEN),allocatable,dimension(:,:,:),public :: rtmp3Da
   real(RLEN),allocatable,dimension(:,:,:),public :: rtmp3Db
#endif

contains

   subroutine init_bfm(cpllog)
!
! DESCRIPTION
!   Initialise the bfm module
!
! USES
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,            &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,    &
                  NO_BOXES_XY,                           &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_STATES, Depth, NO_D3_BOX_FLUX
#if defined INCLUDE_SEAICE
   use mem, only: NO_D2_BOX_STATES_ICE,  &
                  NO_D2_BOX_DIAGNOSS_ICE, &
                  NO_D2_BOX_FLUX_ICE, &
                  NO_STATES_ICE, NO_BOXES_ICE, NO_BOXES_Z_ICE 
#endif
   use mem, only: NO_D2_BOX_STATES_BEN,  &
                  NO_D2_BOX_DIAGNOSS_BEN, &
                  NO_D2_BOX_FLUX_BEN, &
                  NO_STATES_BEN, NO_BOXES_BEN, NO_BOXES_Z_BEN 

   use global_mem, only: LOGUNIT
   use time, only: bfmtime, outdeltalab

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   CHARACTER(len=*),INTENT(IN),OPTIONAL   :: cpllog

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer                   :: rc,n
   character(len=PATH_MAX)   :: logfname, thistime, tmpname
   namelist /bfm_nml/ bio_calc,bio_setup,bfm_init,bfm_rstctl,   &
                      out_fname,out_dir,out_units,out_title,    &
                      out_delta,out_secs,bioshade_feedback,     &
                      parallel_log,unpad_out,                   &
                      in_rst_fname
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   !---------------------------------------------
   ! Provide sensible values for namelist parameters
   !---------------------------------------------
   bio_calc      = .TRUE.
   bio_setup     = 1
   bfm_init      = 0
   bfm_rstctl    = .FALSE.
   out_fname     = 'BFM'
   in_rst_fname  = 'in_bfm_restart'
   out_rst_fname = 'out_bfm_restart'
   out_dir       = '.'
   out_title     = 'Another great BFM simulation!'
   out_units     = 0
   out_delta     = 100
   out_secs      = 100
   bioshade_feedback = .FALSE.
   unpad_out     = .FALSE.

   !---------------------------------------------
   !  Open and read the namelist
   !---------------------------------------------
   NMLUNIT=GetLun()
   open(NMLUNIT,file='BFM_General.nml',action='read',status='old',err=99)
   read(NMLUNIT,nml=bfm_nml,err=98)
   close(NMLUNIT)

   LOGUNIT=GetLun()
   tmpname = TRIM(out_fname)
#ifdef BFM_PARALLEL
   parallel = .TRUE.
   ! variable parallel_rank must have been assigned previously
   ! in the coupling with the ocean model 
   ! check if logs have to be produced for each process
   ! and provide a different log file name 
   write(str,'(I4.4)') parallel_rank
    if (parallel_log) then
       if (parallel_rank == 0) then
          bfm_lwp = .TRUE.
          logfname = 'bfm.log'
          if ( PRESENT(cpllog) ) logfname='bfm'//TRIM(cpllog)
       else 
          bfm_lwp = .FALSE.
          logfname = '/dev/null'
       end if
    else 
       ! logs are produced for every process
       bfm_lwp = .TRUE.
       logfname = 'bfm_'//str//'.log'
       if ( PRESENT(cpllog) ) logfname='bfm'//TRIM(cpllog)//'_'//str
    end if
    open(LOGUNIT,file=logfname,action='write',  &
        form='formatted',err=100)

   if (bfm_lwp) WRITE(*,*) "BFM is running in Parallel"

   ! provide different file names for each process domain
   !
   ! restart output file
   out_rst_fname=TRIM(out_fname)
   !
   ! restart input file 
   if (bfm_init == 1 ) &
      in_rst_fname=TRIM(in_rst_fname)//'_'//str

   ! data output file
   thistime=outdeltalab(out_delta)
   out_fname=TRIM(out_fname)//'_'//TRIM(thistime)//'_'//TRIM(bfmtime%date0)//'_'//TRIM(bfmtime%dateEnd)//'_bfm_'//str

#else
   logfname = 'bfm.log'
   bfm_lwp = .TRUE.
   open(LOGUNIT,file=logfname,action='write',  &
        form='formatted',err=100)
#endif
   !
   !-------------------------------------------------------
   ! Write to log bfmtime setting
   !-------------------------------------------------------
   LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
   LEVEL1 ' BIOGEOCHEMICAL FLUX MODEL (BFM) ACTIVITY LOG  '
   LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
   LEVEL1 '               BFM System Team                 '
   LEVEL1 ' '
   LEVEL1 ' '
   LEVEL1 'Step 1 - INITIALIZATION (init_bfm)'
   LEVEL1 ' '
#ifdef BFM_PARALLEL
   LEVEL2 'Producing log for process rank:',parallel_rank
#endif
   if ( .not. bio_calc ) &
    write(LOGUNIT,*) 'WARNING =>  BFM IS NOT ACTIVE (bio_calc=FALSE) '

   LEVEL1 ' '
   LEVEL2 'EXPERIMENT NAME : ', trim(tmpname)
   LEVEL1 ' '
   LEVEL2 'EXPERIMENT SETUP :'
   select case (bio_setup)
      case (0)
        LEVEL2 "Using a SeaIce-Pelagic-Benthic coupled setup (bio_setup=0)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
#ifdef INCLUDE_SEAICE
        LEVEL3 'seaice variables =',NO_D2_BOX_STATES_ICE
        LEVEL3 'seaice diagnostic variables=', NO_D2_BOX_DIAGNOSS_ICE
#endif
        LEVEL3 'benthic variables =',NO_D2_BOX_STATES_BEN
        LEVEL3 'benthic diagnostic variables=', NO_D2_BOX_DIAGNOSS_BEN
      case (1) ! Pelagic only
        LEVEL2 "Using only Pelagic component (bio_setup=1)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
      case (2) ! Benthic only
        LEVEL2 "Using only Benthic component (bio_setup=2)"
        LEVEL3 'benthic variables =',NO_D2_BOX_STATES_BEN
        LEVEL3 'benthic diagnostic variables=', NO_D2_BOX_DIAGNOSS_BEN
      case (3) ! Pelagic-Benthic coupling
        LEVEL2 "Using a Pelagic-Benthic coupled setup (bio_setup=3)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
        LEVEL3 'benthic variables =',NO_D2_BOX_STATES_BEN
        LEVEL3 'benthic diagnostic variables=', NO_D2_BOX_DIAGNOSS_BEN
#ifdef INCLUDE_SEAICE
      case (4) ! SeaIce only
        LEVEL2 "Using only Seaice component (bio_setup=4)"
        LEVEL3 'seaice variables =',NO_D2_BOX_STATES_ICE
        LEVEL3 'seaice diagnostic variables=', NO_D2_BOX_DIAGNOSS_ICE
      case (5) ! Pelagic-SeaIce coupling
        LEVEL2 "Using a Pelagic-Seaice coupled setup (bio_setup=5)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
        LEVEL3 'seaice variables =',NO_D2_BOX_STATES_ICE
        LEVEL3 'seaice diagnostic variables=', NO_D2_BOX_DIAGNOSS_ICE
#endif
   end select

   LEVEL2 'Dimensional informations:'
   LEVEL3 'NO_BOXES_X  =',NO_BOXES_X
   LEVEL3 'NO_BOXES_Y  =',NO_BOXES_Y
   LEVEL3 'NO_BOXES_Z  =',NO_BOXES_Z
   LEVEL3 'NO_BOXES    =',NO_BOXES
   LEVEL3 'NO_BOXES_XY =',NO_BOXES_XY
   LEVEL3 'NO_STATES   =',NO_STATES
#ifdef INCLUDE_SEAICE
   LEVEL2 'Dimensional seaice informations:'
   LEVEL3 'NO_BOXES_Z_ICE =',NO_BOXES_Z_ICE
   LEVEL3 'NO_BOXES_ICE   =',NO_BOXES_ICE
   LEVEL3 'NO_STATES_ICE  =',NO_STATES_ICE
#endif
   LEVEL2 'Dimensional benthic informations:'
   LEVEL3 'NO_BOXES_Z_BEN =',NO_BOXES_Z_BEN
   LEVEL3 'NO_BOXES_BEN   =',NO_BOXES_BEN
   LEVEL3 'NO_STATES_BEN  =',NO_STATES_BEN
   LEVEL1 ' '

#ifndef BFM_NEMO
   LEVEL2 'EXPERIMENT INITIALIZATION :'
   select case (bfm_init)
      case (0)
        LEVEL2 'Initial conditions from user setting (bfm_init=0)'
      case (1)
        LEVEL2 'Multiple restart files (bfm_init=1)'
        LEVEL2 '(Required filename ',trim(in_rst_fname),')'
      case (2)
        LEVEL2 'Single restart file of the global domain (bfm_init=2)'
        LEVEL2 '(Required filename ',trim(in_rst_fname),')'
   end select
#endif

#ifndef BFM_STANDALONE
   LEVEL1 ' '
   LEVEL2 'EXPERIMENT TIME SETTINGS :'
   LEVEL2 ' Start Date (yyyymmdd)  : ', trim(bfmtime%date0)
   LEVEL2 ' End Date   (yyyymmdd)  : ', trim(bfmtime%dateEnd)
   LEVEL2 ' Initial step           : ', bfmtime%step0
   LEVEL2 ' Final step             : ', bfmtime%stepEnd
   LEVEL2 ' Timestep (seconds)     : ', bfmtime%timestep
#endif
   
   !WRITE(LOGUNIT,'(1a)') 'NetCDF date  Start Date  End Date  Julianday0 &
   !                      & JuliandayEnd   step0  stepnow  stepEnd  timestep'
   !WRITE(LOGUNIT,'(a10,4x,a8,3x,a8,2x,f10.1,2x,f10.1,1x,i9,i9,i9,i9)') bfmtime
   !
#ifndef BFM_NEMO
   LEVEL2 ' '
   LEVEL2 "OUTPUT SETTINGS "
   LEVEL3 "Output DATA filename is: ",trim(out_fname)
   WRITE(LOGUNIT,'(12x,a,i8,a)') "Model data saved every ", out_delta, &
                      & " time-steps (if -1 save exact 1 month freq.)"

   LEVEL3 ' '
   LEVEL3 "RESTART filename prefix is: ",trim(out_rst_fname)
   if ( unpad_out ) then
     LEVEL3 "Restart file(s) saved with the same frequency of the output data."
     LEVEL3 ' '
     LEVEL3 "WARNING => NO padding of output step if larger than the final one."
   else
     LEVEL3 "Restart file(s) saved at the end of this experiment."
   endif
#endif
   LEVEL1 ' '
   LEVEL1 ' -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- '
   LEVEL1 ' '

   ! dimension lengths used in the netcdf output
   lon_len = NO_BOXES_X
   lat_len = NO_BOXES_Y
   depth_len = NO_BOXES_Z
   ocepoint_len = NO_BOXES
   surfpoint_len = NO_BOXES_XY
   botpoint_len = NO_BOXES_XY

   !---------------------------------------------
   ! Allocate arrays with attributes of state
   ! variables
   !---------------------------------------------
   ! total number of output states
   n =  NO_D3_BOX_STATES     + NO_D3_BOX_FLUX     + NO_D3_BOX_DIAGNOSS  + &
        NO_D2_BOX_DIAGNOSS
#ifdef INCLUDE_SEAICE
   n = n + NO_D2_BOX_STATES_ICE + NO_D2_BOX_FLUX_ICE + NO_D2_BOX_DIAGNOSS_ICE
#endif
   n = n + NO_D2_BOX_STATES_BEN + NO_D2_BOX_FLUX_BEN + NO_D2_BOX_DIAGNOSS_BEN

   allocate(var_ids(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_ids'
   var_ids=0;
   allocate(var_ave(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_ave'
   var_ave=.false.;
   allocate(var_names(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_names'
   allocate(var_units(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_units'
   allocate(var_long(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_long'
   ! temporary diagnostic variable
   allocate(c1dim(1:NO_BOXES),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating c1dim'

   return
98 FATAL 'I could not read BFM_General.nml'
   stop 'init_bfm'
99 LEVEL2 'I could not open BFM_General.nml'
   LEVEL2 'Simulation starting without the BFM'
   bio_calc = .false.
100 FATAL 'Cannot create log file: ',trim(logfname)
   stop 'init_bfm'

  end subroutine init_bfm

!-----------------------------------------------------------------------

   subroutine update_save_delta(outdelta,savedelta,timedelta)
!
! DESCRIPTION
!   Dynamically set the output stepping for saving data based on outdelta value
!   NOTE: if outdelta is a negative number then outputs the real monthly data
!
! USES
   use time
   use global_mem, only: RLEN,LOGUNIT
   use constants,  only: SEC_PER_DAY
 
   implicit none
   integer,intent(IN)     :: outdelta
   integer,intent(INOUT)  :: savedelta
   real(RLEN),intent(OUT) :: timedelta
   real(RLEN)             :: julian1, julian2  
   integer                :: yyyy,mm,dd,hh,nn,tmptime

   if ( bfmtime%stepnow .eq. bfmtime%stepEnd ) return
   !
   ! if outdelta is finite use it as default output stepping
   if ( outdelta .ge. 0 ) then 
      tmptime = outdelta
      savedelta = savedelta + outdelta
   endif 
   !
   ! if outdelta is negative dinamically set the output to end of the month
   if ( outdelta .lt. 0 ) then
      julian1 = bfmtime%time0 + ( (bfmtime%stepnow - bfmtime%step0) * bfmtime%timestep / SEC_PER_DAY)
      call calendar_date(julian1,yyyy,mm,dd,hh,nn)
      call julian_day(yyyy,mm,eomdays(yyyy,mm),24,0,julian2)
       
      tmptime = int( julian2 - julian1 ) * int(SEC_PER_DAY) / bfmtime%timestep
      savedelta = savedelta + int( julian2 - julian1 ) * int(SEC_PER_DAY) / bfmtime%timestep

      write(LOGUNIT,*)
      write(LOGUNIT,*) 'update_save_delta : Output will be saved for the real monthly value at step ', savedelta 
      write(LOGUNIT,*) 
   endif
   !
   ! check if output is after the end of simulation and adjust the value
   if ( bfmtime%stepEnd .lt. savedelta .and. .NOT. unpad_out ) then
      tmptime = bfmtime%stepEnd - ( savedelta - tmptime )  
      savedelta = bfmtime%stepEnd  
      write(LOGUNIT,*) 'update_save_delta : Last output saving is beyond the end of the simulation.'
      write(LOGUNIT,*) 'update_save_delta : Output saved at the end of the experiment at ', savedelta
      write(LOGUNIT,*)
   endif
   ! set timestep for time output 
   timedelta = real(savedelta,RLEN)
   if (ave_ctl) timedelta = real(savedelta,RLEN) - (real(tmptime,RLEN) / 2.0) 

   end subroutine  update_save_delta

!-----------------------------------------------------------------------

   function find(vector,nt)
!
! DESCRIPTION
!   Finds the location of true elements in logical arrays
!
! USES
   implicit none
!
! !INPUT PARAMETERS:
   logical,intent(IN) :: vector(:)
   integer,intent(IN) :: nt   ! number of true elements in vector
                              ! nt = count(vector)
                              ! enter as an argument for check
! !OUTPUT PARAMETERS:
   integer            :: find(nt)
!
! !LOCAL VARIABLES:
   integer            :: l,m,n
!
!EOP
!-----------------------------------------------------------------------

    if (nt /= count(vector)) stop '#### Error in find: check the input array ####'
    m = size(vector,1)
    n = 1
    do l = 1,m
      if ( vector(l) ) then
        find(n) = l
        n = n + 1
      end if
    end do

   return
   end function find

!-----------------------------------------------------------------------

 integer function GetLun ()
!
! DESCRIPTION:
!   Adapted from Lionel, Shepherd, Clodius, Page, Drummond.
!   and others as posed at comp.lang.fortran on 1997-09-01
      implicit none
      logical :: exs, opn
      integer :: i
      getlun = -1  ! returned if no units are available.
      i = bfm_file_FirstUnit
 L1:  do
        inquire (unit=i,exist=exs,opened=opn)
          if (exs .and. .not. opn) then
            getlun = i
            exit L1
          end if
        i = i + 1
      end do L1
      return
      stop "There are no free Fortran logical units available."
 end function GetLun

!-----------------------------------------------------------------------

  end module api_bfm

