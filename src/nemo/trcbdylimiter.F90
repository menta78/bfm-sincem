module trcbdylimiter
  implicit none
  private
  public :: bdy_limit_tracers

  integer, allocatable :: mbdy(:,:) ! bdy mask
  integer, allocatable :: mbdylim(:,:) ! limiter mask
  logical, allocatable :: statmask(:,:,:), wet3d(:,:,:), rim3d(:,:,:)
  real(8), save :: trcount
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
      np = size(idx_bdy(ibdy)%nbi, 1)
  
      do ip = 1, np
        mbdy(idx_bdy(ibdy)%nbi(ip, 1), idx_bdy(ibdy)%nbj(ip, 1)) = 1
      end do
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

    allocate(mbdylim(nx, ny))

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
    call init_trcount()
    initialized = .true.
  end subroutine init_bdylim

  subroutine init_trcount
    USE dom_oce, ONLY: tmask
    USE lib_fortran, ONLY: glob_sum
    integer   :: nx, ny, nz
    real(8),    allocatable :: statmask_r(:,:,:)

    nx = size(tmask,1)
    ny = size(tmask,2)
    nz = size(tmask,3)

    wet3d = (tmask /= 0)                     ! same shape as tra(:,:,:,1)
    rim3d = spread(mbdylim /= 0, 3, nz)      ! replicate mbdylim vertically
    statmask = wet3d .and. .not. rim3d
    allocate(statmask_r(nx, ny, nz))
    statmask_r = 0
    where (statmask)
      statmask_r = 1
    end where

    trcount = glob_sum('bdy_limit_tracers', statmask_r)

    deallocate(statmask_r)
  end subroutine init_trcount

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
    USE lib_mpp, ONLY: mpp_sum, mpprank
    USE lib_fortran, ONLY: glob_sum
    USE, INTRINSIC :: ieee_arithmetic
    real,    intent(inout)  :: tra(:,:,:,:)
    real,    intent(in), optional :: scoeff
    integer :: idx, k, ntr, nx, ny, nz
    real(8),    allocatable :: trak(:,:,:)
    real(8):: s, trsum, trdvsum, trmean, trmad, thrshldhigh, thrshldlow
    ! bdy limiter
    
    if (.not. initialized) call init_bdylim

    ntr = size(tra,4)

    ! default to p90 for Normal using mean abs deviation about mean
    ! k = z0.9 / E|Z| = 1.281551565 / 0.797884561 = 1.606186694
    ! s = 1.606186694
    s = 2.5 ! corresponds to 97.5% of a normal distribution
    if (present(scoeff)) then
      s = scoeff
    end if

   !if (mpprank .eq. 0) WRITE(6,*) 'called mpp_sum trcbdylim:count, count=', trcount
   !call flush(6)

    allocate(trak(nx,ny,nz))

    do k = 1, ntr
      ! ---- gather wet points only
      trak = tra(:,:,:,k)
      where (.NOT. statmask .OR. trak /= trak .OR. .NOT. ieee_is_finite(trak))
          trak = 0
      end where

      ! ---- compute mean and mad and threahold
      !trsum = sum(vals)
      !call mpp_sum('trcbdylimiter', trsum)
      trsum = glob_sum('bdy_limit_tracers', trak)
      trmean = trsum/trcount
      trdvsum = glob_sum('bdy_limit_tracers', abs(trak-trmean))
      trmad = trdvsum/trcount
      thrshldhigh = trmean + s*trmad
      thrshldlow = max(trmean - s*trmad, 0.)

     !if (mpprank .eq. 0) then
     !    WRITE(6,*) ''
     !    WRITE(6,*) '  tracer ',k
     !    WRITE(6,*) '  trmean      = ',trmean
     !    WRITE(6,*) '  thrshldhigh = ',thrshldhigh
     !    WRITE(6,*) '  thrshldlow = ',thrshldlow
     !end if
    
      ! ---- trim only inside boundary buffer & only on wet ocean
      where (rim3d .and. wet3d .and. tra(:,:,:,k) > thrshldhigh)
        tra(:,:,:,k) = thrshldhigh
      end where
      where (rim3d .and. wet3d .and. tra(:,:,:,k) < thrshldlow)
        tra(:,:,:,k) = thrshldlow
      end where
    end do
    deallocate(trak)
  end subroutine bdy_limit_tracers

end module trcbdylimiter
