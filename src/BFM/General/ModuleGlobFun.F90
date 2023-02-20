!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: ModuleGlobFun
!
! DESCRIPTION
!   List of general model functions
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
!
! INTERFACE
  MODULE mem_globalfun
!
! USES
  USE global_mem, ONLY: RLEN, ZERO, ONE, BASETEMP
  USE mem_Param,  ONLY: p_small
  USE constants,  ONLY: ZERO_KELVIN

  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  ! SHARED GLOBAL FUNCTIONS (must be below contains)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  contains

    ! Convert values in 0 or 1 according to input field
    elemental FUNCTION INSW(x)

        IMPLICIT NONE
        real(RLEN),intent(IN) :: x
        real(RLEN)            :: INSW

        INSW = ZERO
        if (x > ZERO ) INSW=ONE

    end function INSW

    ! Michaelis-Menten saturation curve
    elemental FUNCTION MM(x, m)

        IMPLICIT NONE
        real(RLEN),intent(IN) :: x, m
        real(RLEN)            :: MM

        MM = x / (x + m)
        if ( m == ZERO ) MM = ONE 

    end function MM

    ! Michaelis-Menten saturation curve at power
    elemental FUNCTION MM_POWER(x, m, p)

        IMPLICIT NONE
        real(RLEN),intent(IN) :: x, m
        integer   ,intent(IN) :: p
        real(RLEN)            :: MM_POWER

        MM_POWER = x**p / ( x**p+ m**p)
        if ( m == ZERO ) MM_POWER = ONE 

    end function MM_POWER

    ! Ramp unary real function
    elemental FUNCTION ERAMP(x, m)

        IMPLICIT NONE
        real(RLEN),intent(IN) :: x, m
        real(RLEN)            :: ERAMP

        ERAMP = ZERO
        if ( x > ZERO ) ERAMP = ONE
        if ( x < m )    ERAMP = x / m

    end function

    elemental FUNCTION PartQ(p, a, b, c)

        IMPLICIT NONE
        real(RLEN),intent(IN) :: p, a, b, c
        real(RLEN)            :: PartQ

        real(RLEN) :: a1, b1, c1, norm, r

        c1 = min(p * (-ONE * log(1.0E-20_RLEN)), c);
        b1 = min(b, c1);
        a1 = min(a, b1);
        r  = ZERO

        if ( a == ZERO ) r = ONE

        if (c1 > ZERO .and. p /= ZERO) then
           norm = ONE - exp((- c1) / p)
           PartQ = (exp( -a1 / p) - exp(- b1 / p)) / norm
        else
           PartQ = r
        endif

    end function

    ! temperature dependency for Q10 function
    ! (with optional argument for base temperature)
    elemental function eTq(t, q10, base)

        IMPLICIT NONE
        real(RLEN),intent(IN)          :: t, q10
        real(RLEN),intent(IN),optional :: base
        real(RLEN)                     :: eTq
        real(RLEN)                     :: tref

        tref = BASETEMP
        if ( present(base) ) tref = base

        eTq = exp( log(q10) * (t-tref) / tref)

    end function eTq

    ! Arrhenius equation for temperature dependency 
    ! (with optional argument for base temperature)
    ! activation energy (aex) is in J mol-1
    elemental function eTa(t, aex, base)

        IMPLICIT NONE
        real(RLEN),intent(IN)          :: t, aex
        real(RLEN),intent(IN),optional :: base
        real(RLEN)                     :: eTa
        real(RLEN)                     :: tref, act

        tref = BASETEMP
        if ( present(base) ) tref = base

        act  = ( aex * 1000._RLEN ) / 8.3_RLEN

        eTa = exp( act * ( (ONE / (tref - ZERO_KELVIN)) - (ONE / ( t - ZERO_KELVIN)) ) )  

    end function eTa


    elemental function IntegralExp(alfa,x)

        IMPLICIT NONE
        real(RLEN),intent(IN)  :: alfa, x
        real(RLEN)             :: IntegralExp

        IntegralExp=(exp(alfa * x) -ONE)/alfa

    end function IntegralExp

    elemental subroutine fixratio(fc, fn, fp, qnc, qpc, rc, rn, rp)
    !==========================================================================
    ! Determine release fluxes due to either C, N, or P limitiation
    ! fc,fn,fp     : total fluxes of C,N,P                    [mol/kg]
    ! qnc          : N:C ratio                                [molN/molC]
    ! qpc          : P:C ratio                                [molP/molC]
    ! rc,rn,rp     : release fluxes of C,N,P                  [mol/kg]
    !==========================================================================
       IMPLICIT NONE
       !
       real(RLEN), intent(in)  :: fc, fn, fp, qnc, qpc
       real(RLEN), intent(out) :: rc, rn, rp
       !
       rc = ZERO ; rn = ZERO ; rp = ZERO
       !
       ! carbon
       rc = max( max(fc-fp/qpc,fc-fn/qnc) , ZERO )
       ! 
       ! nitrogen
       rn = max( fn - (fc-rc)*qnc, ZERO )
       ! 
       ! phosphorous
       rp = max( fp - (fc-rc)*qpc, ZERO )
       !
       return
    end subroutine fixratio

    elemental function analytical_ic(z, z1, v1, z2, v2)
    !==========================================================================
    ! Create analytical field depth profile using 2 layer input data
    ! NOTE: if z2,v2 are zeros create only uniform value in the upper layer
    ! z            : Target depth                             [m]
    ! z1, z2       : Depth of upper and lower reference layer [m]
    ! v1, v2       : Field value at upper and lower layers    [field unit]
    !==========================================================================
        IMPLICIT NONE
        !
        real(RLEN),intent(IN) :: z, z1, v1, z2, v2
        real(RLEN)            :: analytical_ic
        real(RLEN)            :: alpha, fout
        !
        fout = ZERO
        !
        ! below lower layer (do it first to solve case z2 = 0)
        if ( z .gt. z2 ) fout = v2
        !
        ! above upper layer
        if ( z .le. z1 ) fout = v1
        !
        ! within upper and lower layer (linear interpolation)
        alpha = (v2-v1) / (z2-z1 +2.E-15)
        if ( (alpha .le. 1.E15) .AND. ( z .gt. z1 .and. z .le. z2 ) )  &
            fout = v1 + alpha * (z-z1)
        !
        analytical_ic = max( fout , p_small ) 
        !
        return
    end function analytical_ic

  end module mem_globalfun

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
