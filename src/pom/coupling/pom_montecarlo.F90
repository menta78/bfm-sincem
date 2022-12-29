MODULE MONTECARLO

USE GLOBAL_MEM, ONLY: RLEN

IMPLICIT NONE

PRIVATE

PUBLIC PROFILE_MONTECARLO_NORMAL


INTEGER, PARAMETER    :: ERFTABSIZE = 10000
REAL(RLEN)            :: ERFTABDELTA = REAL(2, RLEN)/ERFTABSIZE
REAL(RLEN)            :: ERFTAB(0:ERFTABSIZE)
LOGICAL               :: INITIALIZED = .FALSE.

CONTAINS


SUBROUTINE PROFILE_MONTECARLO_NORMAL(MEAN_PROF, VARIANCE_PROF, OUT_PROF)
      REAL(RLEN), INTENT(IN),  DIMENSION(:) :: MEAN_PROF, VARIANCE_PROF
      REAL(RLEN), INTENT(OUT), DIMENSION(:) :: OUT_PROF

      REAL(RLEN) :: PR, ERFI

      IF (.NOT. INITIALIZED) CALL INITIALIZE()

      ! randomizing the probability
      PR = RAND()
      ! computing the inverse normal CDF:
      ! x = \sqrt{2 \sigma^2} \erf^-1(2 p - 1) + \mu
      ! where mu and sigma^2 are the mean and the variance
      ! erf^-1 is the inverse error function
      ERFI = ERFINV(2*PR - 1)
      OUT_PROF = SQRT(2*VARIANCE_PROF)*ERFI + MEAN_PROF
END SUBROUTINE


REAL(RLEN) FUNCTION ERFINV(YY)
      REAL(RLEN), INTENT(IN) :: YY
      INTEGER :: IY

      IY = FLOOR((YY + 1)/ERFTABDELTA)
      ERFINV = ERFTAB(IY)
END FUNCTION



SUBROUTINE INITIALIZE()
      REAL(RLEN) :: YYTAB, YYERF, MINX, MAXX, DELTAX, XX
      INTEGER :: IX, IY, XYFCT

      ! seeding the randomizer
      CALL SRAND(10000)

      ! generating a y->x table for the erf function for rapid inversion
      YYTAB = -1 + ERFTABDELTA
      MINX = -2.7
      MAXX = 2.7
      XYFCT = 10
      DELTAX = (MAXX - MINX)/(ERFTABSIZE*XYFCT)
      ERFTAB(0) = MINX
      IY = 1
      DO IX = 1,(ERFTABSIZE)*XYFCT
        XX = MINX + IX*DELTAX
        YYERF = ERF(XX)
        IF (YYERF .GE. YYTAB) THEN
           ERFTAB(IY) = XX
           YYTAB = YYTAB + ERFTABDELTA
           IY = IY + 1
        END IF
      END DO
      ERFTAB(ERFTABSIZE) = MAXX

      INITIALIZED = .TRUE.
END SUBROUTINE



END MODULE

