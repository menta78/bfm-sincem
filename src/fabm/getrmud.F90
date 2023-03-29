      subroutine getrmud(sunz,mud)
      USE fabm_types,     ONLY: rk
      IMPLICIT NONE 
!  Computes average cosine for direct irradiance just below
!  sea surface
      real(rk),  intent(in) :: sunz
      real(rk), intent(out) :: mud
      real(rk), parameter   :: refrac_idx = 1.341_rk
      real(rk)              :: rad 
      real(rk)              :: rsza, sinszaw
      real(rk)              :: szaw, rmudl
 
!  Compute average cosine for direct irradiance in the water 
!  column given solar zenith angle (in degrees) at surface.
      rad    = 180.0D0/dacos(-1.0D0) ! radians
      rsza = sunz/rad
      sinszaw = sin(rsza)/refrac_idx
      szaw = asin(sinszaw)
      mud  = cos(szaw)
!     rmudl = 1.0D0/cos(szaw)   !avg cosine direct (1 over)
!     rmud = min(rmudl,1.5D0)
!     rmud = max(rmud,0.0D0)
      
 
      return
      end
