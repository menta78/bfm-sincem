! SUBROUTINE LATERAL_FLUX:       
! If the nutrient surface boundary condition is defined with a flux of
! nutrients (NUTSBC_MODE == 1), the same amount of water entering the
! system at the surfcace must leave it laterally. This translates into
! an outgoing flux of constituents.
! For now a constant profile of lateral flux is assumed
      SUBROUTINE SUBTRACT_LATERAL_FLUX(CONC)
              USE GLOBAL_MEM, ONLY: RLEN, ZERO
              USE POM, ONLY: DTI, KB, H
              USE SERVICE, ONLY: DISSURF, CURRENTS_SPEED_PROF
              IMPLICIT NONE

              ! CONC: concentration of a constituent
              REAL(RLEN), INTENT(INOUT) :: CONC(KB)  
              REAL(RLEN) :: WT_LFLUX_BY_METER, CNST_LFLUX_BY_METER
              INTEGER :: IK

              WT_LFLUX_BY_METER = DISSURF/H ! water lateral flux by meter. Kg/m3/s
              DO IK = 1,KB-1
                 CNST_LFLUX_BY_METER = WT_LFLUX_BY_METER/1000*CONC(IK) !constituent lateral flux by depth meter. mmol/m3/s
                 ! adjusting by the vertical profile of currents speed
                 CNST_LFLUX_BY_METER = CNST_LFLUX_BY_METER*CURRENTS_SPEED_PROF(IK) 
                 !Units are divided by m2 as the 1d model is by surface unit
                 ! If we multiply CNST_LFLUX_BY_METER by time, we can compare it immediately with CONC

                 ! Applying a forward updstream scheme, assuming CNST_LFLUX_BY_METER is small
                 CONC(IK) = MAX(CONC(IK) - CNST_LFLUX_BY_METER*DTI, ZERO)
              END DO
      END SUBROUTINE

