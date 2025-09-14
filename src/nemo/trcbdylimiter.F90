module trcbdylimiter
  implicit none
  private
  public :: bdy_limit_tracers

  integer, allocatable :: mbdy(:,:) ! bdy mask
  integer, allocatable :: mbdylim(:,:) ! limiter mask
  logical, save :: initialized = .FALSE.

contains
   
   
   
  subroutine init_mbdy
    ! Build a 2-D open-boundary mask (mbdy) on the T-grid from NEMO BDY lists.
    ! - Marks 1 at all BDY points (optionally up to nn_rimwidth if nbr is provided).
    ! - Assumes BDY lists refer to T-grid (use U/V lists similarly if you need them).
    use par_oce, only: jpi, jpj
    ! The following USE matches typical NEMO-4 BDY exposure.
    ! You may need to adapt names to your branch (idx structure field names vary slightly).
    use bdy_oce, only: nb_bdy, nn_rimwidth         ! number of boundary sets, rim width
    use bdy_oce, only: idx_bdy                      ! type/arrays describing BDY geometry
    implicit none
  
    integer :: nx, ny
    integer :: ibdy, ip, np
    logical :: has_nbr
  
    nx = jpi
    ny = jpj
  
    allocate(mbdy(nx,ny))
    mbdy = 0
  
    ! ---- Loop through boundary sets and mark their T-grid points ----
    do ibdy = 1, nb_bdy
      ! Depending on your NEMO version, the T-grid lists are exposed as one of:
      !   idx_bdy(ibdy)%nblen(1)         : number of T points
      !   idx_bdy(ibdy)%nbi(:)      : i indices
      !   idx_bdy(ibdy)%nbj(:)      : j indices
      !   idx_bdy(ibdy)%nbr(:)      : rim index (1..nn_rimwidth); may not exist on older branches
      !
      ! If your branch uses a single structure (e.g., idx(ibdy)%nbi(:,igrid)) adapt accordingly.
  
      !np = idx_bdy(ibdy)%nblen(1)
      np = size(idx_bdy(ibdy)%nbi, 1)
      has_nbr = associated(idx_bdy(ibdy)%nbr)
  
      if (np > 0) then
        if (has_nbr) then
          do ip = 1, np
            if (idx_bdy(ibdy)%nbr(ip, 1) <= nn_rimwidth(ibdy)) then
              mbdy(idx_bdy(ibdy)%nbi(ip, 1), idx_bdy(ibdy)%nbj(ip, 1)) = 1
            end if
          end do
        else
          do ip = 1, np
            mbdy(idx_bdy(ibdy)%nbi(ip, 1), idx_bdy(ibdy)%nbj(ip, 1)) = 1
          end do
        end if
      end if
    end do
  end subroutine init_mbdy


  !--------------------------------------------------------------------
  ! init_bdylim:
  !   Build/overwrite mbdylim so it has the same shape as mbdy.
  !   mbdylim = 1 for any cell within 'radius' (default 5) grid cells
  !   of at least one mbdy==1 point (Chebyshev distance).
  !--------------------------------------------------------------------
  subroutine init_bdylim(radius)
    implicit none
    integer, intent(in), optional :: radius
    integer :: nx, ny, r
    integer :: i, j, ii, jj
    integer :: i1, i2, j1, j2

    call init_mbdy()

    nx = size(mbdy, 1)
    ny = size(mbdy, 2)
    r  = 10
    if (present(radius)) r = max(0, radius)

    if (.not. allocated(mbdylim)) then
      allocate(mbdylim(nx, ny))
    else
      if (size(mbdylim,1) /= nx .or. size(mbdylim,2) /= ny) then
        deallocate(mbdylim)
        allocate(mbdylim(nx, ny))
      end if
    end if

    mbdylim = 0

    do j = 1, ny
      j1 = max(1, j - r)
      j2 = min(ny, j + r)
      do i = 1, nx
        i1 = max(1, i - r)
        i2 = min(nx, i + r)

        do jj = j1, j2
          do ii = i1, i2
            if (mbdy(ii, jj) /= 0) then
              mbdylim(i, j) = 1
              exit
            end if
          end do
          if (mbdylim(i, j) == 1) exit
        end do

      end do
    end do
    initialized = .true.
  end subroutine init_bdylim

  !--------------------------------------------------------------------
  ! bdy_limit_tracers:
  !   For each tracer k, compute the mean + s*MAD
  !   (nx,ny) field and trim values to that threshold ONLY where
  !   mbdylim==1. Optionally pass a different percentile in (0,1).
  !
  !   Arguments:
  !     tra     [inout] real(:,:,:)  tracer array (nx,ny,ntr)
  !     pfract  [in]    real, optional (default 0.90)
  !
  !   Notes:
  !     - Percentile uses nearest-rank definition: ceil(p * N) after sort.
  !     - No wet mask used here; if you have land points, ensure they are
  !       set to reasonable fill (or prefilter) before calling.
  !--------------------------------------------------------------------
  subroutine bdy_limit_tracers(tra, scoeff)
    USE par_oce,  only: wp
    USE dom_oce, ONLY: tmask
    USE lib_mpp, ONLY: mpp_sum, mpprank
   !USE lib_fortran, ONLY: glob_sum
    USE, INTRINSIC :: ieee_arithmetic
    real,    intent(inout)  :: tra(:,:,:,:)
    real,    intent(in), optional :: scoeff
    logical, allocatable :: wet3d(:,:,:), rim3d(:,:,:), statmask(:,:,:)
    integer :: nwet, idx, k, ntr, nx, ny, nz, trcnt_loc
    real(8),    allocatable :: vals(:), trak(:,:,:)
    real(8):: s, trcount, trsum, trdvsum, trmean, trmad, thrshldhigh, thrshldlow
    ! bdy limiter
    
    if (.not. initialized) call init_bdylim

    ntr = size(tra,4)
    nx = size(tra,1)
    ny = size(tra,2)
    nz = size(tra,3)

    wet3d = (tmask /= 0)                     ! same shape as tra(:,:,:,1)
    rim3d = spread(mbdylim /= 0, 3, nz)      ! replicate mbdylim vertically
    statmask = wet3d .and. .not. rim3d

    ! default to p90 for Normal using mean abs deviation about mean
    ! k = z0.9 / E|Z| = 1.281551565 / 0.797884561 = 1.606186694
    s = 1.606186694
    !s = 2.5 ! corresponds to 97.5% of a normal distribution
    if (present(scoeff)) then
      s = scoeff
    end if

    trcnt_loc = count(statmask)
    trcount = trcnt_loc
   !if (mpprank .eq. 0) WRITE(6,*) 'calling mpp_sum trcbdylim:count'
    call mpp_sum('trcbdylimiter', trcount)
   !if (mpprank .eq. 0) WRITE(6,*) 'called mpp_sum trcbdylim:count, count=', trcount
   !call flush(6)

    nwet = count(wet3d)
    trcnt_loc = count(statmask)
    if (trcnt_loc > 0) then
      allocate(vals(trcnt_loc))
    else
      allocate(vals(1))
    end if
    allocate(trak(nx,ny,nz))

    do k = 1, ntr
      ! ---- gather wet points only
      trak = tra(:,:,:,k)
      vals = 0
      if (trcnt_loc > 0) vals(:) = pack(trak, mask=statmask)
      where (vals /= vals .OR. .NOT. ieee_is_finite(vals))
          vals = 0
      end where

      ! ---- compute mean and mad and threahold
      trsum = sum(vals)
     !if (mpprank .eq. 0) WRITE(6,*) 'calling mpp_sum trcbdylim:sum'
      call mpp_sum('trcbdylimiter', trsum)
     !if (mpprank .eq. 0) WRITE(6,*) 'called mpp_sum trcbdylim:sum, sum=', trsum
     !call flush(6)
      trmean = trsum/trcount
      trdvsum = sum(abs(vals-trmean))
     !if (mpprank .eq. 0) WRITE(6,*) 'calling mpp_sum trcbdylim:mad'
      call mpp_sum('trcbdylimiter', trdvsum)
     !if (mpprank .eq. 0) WRITE(6,*) 'called mpp_sum trcbdylim:mad, trdvsum', trdvsum
     !call flush(6)
      trmad = trdvsum/trcount
      thrshldhigh = trmean + s*trmad
      thrshldlow = max(trmean - s*trmad, 0)

      if (mpprank .eq. 0) 
          WRITE(6,*) ''
          WRITE(6,*) 'tracer ',k
          WRITE(6,*) '  trmean      = ',trmean
          WRITE(6,*) '  thrshldhigh = ',thrshldhigh
          WRITE(6,*) '  thrshldlow = ',thrshldlow
      end if
    
      ! ---- trim only inside boundary buffer & only on wet ocean
      where (rim3d .and. wet3d .and. tra(:,:,:,k) > thrshld)
        tra(:,:,:,k) = thrshld
      end where
    end do
    deallocate(vals)
    deallocate(trak)
  end subroutine bdy_limit_tracers


  subroutine mpp_sum_dp(val)
    ! mpp_sum from nemo does not seem to work here
    USE mpi
    USE lib_mpp, ONLY: mpi_comm_oce
    implicit none
    real(8), intent(inout) :: val
    real(8) :: tmp
    integer :: ierr
    logical :: inited

    call MPI_Initialized(inited, ierr)
    if (.not. inited) return
  
    tmp = val
    call MPI_Allreduce(tmp, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_oce, ierr)
  end subroutine mpp_sum_dp

end module trcbdylimiter
