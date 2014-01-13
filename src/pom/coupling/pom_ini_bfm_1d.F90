#include "INCLUDE.h"
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pom_ini_bfm.F90
!
! !INTERFACE:
subroutine pom_ini_bfm_1d
!
! !DESCRIPTION:
!  Initialise the BFM in POM
!  Main communication of array dimensions between
!  BFM and POM
!  Initialisation of variables and netcdf output
! THIS IS THE 1D VERSION
!
!
#ifdef NOPOINTERS
   use mem
#else
!   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
!                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
!                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
!                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
!                  NO_STATES,Depth,Depth_ben, D3STATE, D2STATE
     use mem
#endif
   use POM
   use global_mem, only:RLEN,ZERO,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
   use api_bfm
   use Service
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   use netcdf_bfm, only: init_netcdf_rst_bfm,read_rst_bfm

   IMPLICIT NONE
!
! !LOCAL VARIABLES:

   integer    :: i,j,k,ll,count
   integer    :: status
   integer,parameter    :: namlst=10,unit=11
   integer,allocatable  :: ocepoint(:),surfpoint(:),botpoint(:)
   real(RLEN),allocatable :: tmask(:,:,:)
!EOP
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------
! Model timestep
!----------------------------------------------------------------------
      deltat=dti
!----------------------------------------------
! Iterations needed for an idays simulations
!---------------------------------------------
      nitend=iend
!---------------------------------------------
! out_delta = saving frequency 
!---------------------------------------------
      out_delta = savef

!---------------------------------------------
! Set the BFM dimensions
!---------------------------------------------
      NO_BOXES=KB-1
      NO_BOXES_X=1
      NO_BOXES_Y=1
      NO_BOXES_Z=KB-1
      NO_BOXES_XY=1
!
   allocate(SEAmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))  ! allocate  masks for
   allocate(BOTmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))  ! array packing
   allocate(SRFmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
   allocate(tmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
   SEAmask = .TRUE.
   BOTmask = .TRUE.
   SRFmask = .TRUE.
   tmask=1.
!
   allocate(ZEROS(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))   ! allocate ancillary pack mask
   ZEROS = ZERO
!
#ifdef INCLUDE_BEN
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES_BEN 
#else
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES
#endif
!
   !-------------------------------------------------------
   ! Compressed coordinates for netcdf output
   !-------------------------------------------------------
   lon_len = NO_BOXES_X 
   lat_len = NO_BOXES_Y
   allocate(ocepoint(NO_BOXES))
   ocepoint = 1 
   !-------------------------------------------------------
   ! allocate the arrays for surface and bottom boundary conditions
   ! At the moment only for inorganic components (nutrients and gases)
   !-------------------------------------------------------
!      allocate(bc2D_bio_surf(NO_D3_BOX_STATES,NO_BOXES),stat=status)
!      if (status /= 0)
!     &  stop "# FATAL ERROR in ini_trc_bfm: Error allocating bc2D_bio_surf"
!      bc2D_bio_surf = ZERO
!      allocate(bc2D_bio_bot(NO_D3_BOX_STATES,NO_BOXES),stat=status)
!      if (status /= 0)
!     &  stop "# FATAL ERROR in ini_trc_bfm: Error allocating bc2D_bio_bot"
!      bc2D_bio_surf = ZERO

   !-------------------------------------------------------
   ! Prepares the array containing the indices of the
   ! elements in pelagic BFM 1D arrays that have a benthic layer
   !-------------------------------------------------------
   allocate(BOTindices(NO_BOXES_XY)); BOTindices = 30
!   !-------------------------------------------------------
!   ! Prepares the array containing the indices of the
!   ! elements in pelagic BFM 1D arrays that have a surface layer
!   !-------------------------------------------------------
   allocate(SRFindices(NO_BOXES_XY)); SRFindices = 1
   !---------------------------------------------
   ! Initialise ancillary arrays for output
   !---------------------------------------------
   call init_bfm(namlst)
   !---------------------------------------------
   ! Initialise state variable names and diagnostics
   !---------------------------------------------
   call set_var_info_bfm
   !---------------------------------------------
   ! Allocate memory and give initial values
   !---------------------------------------------
   ! the argument list is mandatory with BFM
   call init_var_bfm(namlst,'bfm.nml',unit,bio_setup)
        do k = 1 , NO_BOXES_Z
                    Depth(k) = dz(k)*h
        end do
   !---------------------------------------------
   ! initialise netcdf output
   ! MAV: CHECK THE MAPPING OF OCEPOINT
   !---------------------------------------------
   call calcmean_bfm(INIT)
   call init_netcdf_bfm(title=out_title,start_time='01-01-0000', &
             time_unit=0,lat=alat,lon=alon,z=z,dz=dz,      &
             oceanpoint=ocepoint,                  &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/))!  &     (G)
!             ,mask3d=tmask )
   !---------------------------------------------
   ! Allocate and initialise additional 
   ! integration arrays
   !---------------------------------------------
   allocate(D3STATEB(NO_D3_BOX_STATES,NO_BOXES))
!   allocate(D3STATEB(NO_BOXES,NO_D3_BOX_STATES))

! #ifdef INCLUDE_BEN
!    allocate(D2STATEB(NO_D2_BOX_STATES_BEN,NO_BOXES_XY))
! #endif   
   !---------------------------------------------
   ! Override homogeneous initialisation 
   ! MAV: crappy stuff, to be improved
   !---------------------------------------------
   call set_initial_conditions
!   call flush(6)

   call init_save_bfm
!---------------------------------------------
! Initialise prior time step for leap-frog
!---------------------------------------------
!    D3STATEB = D3STATE
! #ifdef INCLUDE_BEN
!    D2STATEB = D2STATE_BEN
! #endif
!---------------------------------------------
!  Read restart (Bfm_init = 1 in bfm.nml)
!---------------------------------------------
!  call init_netcdf_rst_bfm(rst_fname)                  (G)
   if (bfm_init == 1) then
      call read_rst_bfm(rst_fname)
   end if
!   
!write(6,*) 'DOPO SET_INI', D3STATE(ppP1c,1), D3STATE(ppP1c,NO_BOXES),maxval(D3STATE(ppP1c,:)),minval(D3STATE(ppP1c,:))
   return
!
end subroutine pom_ini_bfm_1d
!EOC


