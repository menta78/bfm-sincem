! SUBROUTINE LATERAL_FLUX:       
! If the nutrient surface boundary condition is defined with a flux of
! nutrients (NUTSBC_MODE == 1), the same amount of water entering the
! system at the surfcace must leave it laterally. This translates into
! an outgoing flux of constituents.
! For now a constant profile of lateral flux is assumed
      SUBROUTINE SUBTRACT_LATERAL_FLUX(CONC, LAMBDA)
              USE GLOBAL_MEM, ONLY: RLEN, ZERO
              USE POM, ONLY: DTI, KB, LENGTH_SCALE
              IMPLICIT NONE

              ! CONC: concentration of a constituent
              REAL(RLEN), INTENT(INOUT) :: CONC(KB)  
              REAL(RLEN), INTENT(IN)    :: LAMBDA ! relaxation rate (1/s)
              REAL(RLEN) :: CNST_LFLUX
              INTEGER :: IK

              DO IK = 1,KB-1
                 CNST_LFLUX = LAMBDA*CONC(IK) !constituent lateral flux by depth meter. mmol/m2/s

                 ! Applying a forward updstream scheme, assuming CNST_LFLUX is small
                 CONC(IK) = MAX(CONC(IK) - CNST_LFLUX*DTI, ZERO)
              END DO
      END SUBROUTINE

! SUBROUTINE LATERAL_FLUX2:       
! Like SUBTRACT_LATERAL_FLUX, but with a tendency profile to converge to
      SUBROUTINE SUBTRACT_LATERAL_FLUX2(CONC, TND_CONC, LAMBDA)
              USE GLOBAL_MEM, ONLY: RLEN, ZERO
              USE POM, ONLY: DTI, KB, LENGTH_SCALE
              IMPLICIT NONE

              ! CONC: concentration of a constituent
              REAL(RLEN), INTENT(INOUT) :: CONC(KB)  
              ! TND_CONC: tendency concentration of a constituent
              REAL(RLEN), INTENT(IN) :: TND_CONC(KB)
              REAL(RLEN), INTENT(IN)    :: LAMBDA ! relaxation rate (1/s)
              REAL(RLEN) :: CNST_LFLUX
              INTEGER :: IK

              DO IK = 1,KB-1
                 CNST_LFLUX = LAMBDA*(CONC(IK)-TND_CONC(IK)) !constituent lateral flux by depth meter. mmol/m2/s

                 ! Applying a forward updstream scheme, assuming CNST_LFLUX is small
                 CONC(IK) = MAX(CONC(IK) - CNST_LFLUX*DTI, ZERO)
              END DO
      END SUBROUTINE
