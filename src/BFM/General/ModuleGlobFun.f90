!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleGlobFun
!
! DESCRIPTION
!   List of general model functions

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE mem_globalfun
!
! !USES:
  USE global_mem, ONLY: RLEN, ZERO, ONE, BASETEMP
  USE mem_Param,  ONLY: p_small

!  
!
! !AUTHORS
!   mfstep/ERSEM team
!
! !REVISION_HISTORY
!   --------
!
! COPYING
!   
!   Copyright (C) 2020 BFM System Team (bfm_st@cmcc.it)
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  ! SHARED GLOBAL FUNCTIONS (must be below contains)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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

    end function MM

    ! Michaelis-Menten saturation curve at power
    elemental FUNCTION MM_POWER(x, m, p)

        IMPLICIT NONE
        real(RLEN),intent(IN) :: x, m
        integer   ,intent(IN) :: p
        real(RLEN)            :: MM_POWER

        MM_POWER = x**p / ( x**p+ m**p)

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
    elemental function eTq(t, q10)

        IMPLICIT NONE
        real(RLEN),intent(IN) :: t, q10
        real(RLEN)            :: eTq

        eTq = exp( log(q10) * (t-BASETEMP) / BASETEMP)

    end function eTq

    elemental function IntegralExp(alfa,x)

        IMPLICIT NONE
        real(RLEN),intent(IN)  :: alfa, x
        real(RLEN)             :: IntegralExp
        
        IntegralExp=(exp(alfa * x) -ONE)/alfa

    end function IntegralExp

    elemental function nutlim(fc, fn, fp, qnc, qpc, c, n, p)
    !==========================================================================
    ! Determine whether C, N, or P is the limiting element for a PFT
    ! fc,fn,fp     : total fluxes of C,N,P                    [mol/kg]
    ! qnc          : N:C ratio                                [molN/molC]
    ! qpc          : P:C ratio                                [molP/molC]
    ! c,n,p        : constituents indexes                     [-]
    ! nutlim       : limiting consituent index                [-]
    !
    ! Configurations of fluxes for C, N, P
    !          fc        fn         fp
    !   1      <0        <0         <0
    !   2      <0        <0         >0
    !   3      <0        >0         <0
    !   4      <0        >0         >0
    !   5      >0        <0         <0
    !   6      >0        <0         >0
    !   7      >0        >0         <0
    !   8      >0        >0         >0
    !==========================================================================
       IMPLICIT NONE
       !
       real(RLEN), intent(in) :: fc, fn, fp, qnc, qpc
       integer,    intent(in) :: c, n, p
       integer                :: nutlim
       real(RLEN) :: pu_n, pu_p, pq_n, pq_p
       !
       pu_n = fn / ( p_small + fc ) ; pq_n = pu_n / qnc
       pu_p = fp / ( p_small + fc ) ; pq_p = pu_p / qpc
       !
       nutlim = c
       !
       ! CASE 1
       if ( (fc <0) .AND. (fn <0) .AND. (fp <0) ) then
          if ( pq_p>pq_n .OR. abs(pq_p-pq_n)<p_small ) then
             if (abs(pu_p) > qpc) nutlim = p
          else
             if (abs(pu_n) > qnc) nutlim = n
          endif
       endif
       ! CASE 2
       if ( (fc <0) .AND. (fn <0) .AND. (fp >0) ) then
           if (abs(pu_n) > qnc) nutlim = n
       endif
       ! CASE 3
       if ( (fc <0) .AND. (fn >0) .AND. (fp <0) ) then
          if (abs(pu_p) > qpc) nutlim = p
       endif
       ! CASE 4
       if ( (fc <0) .AND. (fn >0) .AND. (fp >0) ) nutlim = c
       ! CASE 5
       if ( (fc >0) .AND. (fn <0) .AND. (fp <0) ) then
          if ( pq_p<pq_n .OR. abs(pq_p-pq_n)<p_small ) then
             nutlim = p
          else
             nutlim = n
          endif
       endif
       ! CASE 6
       if ( (fc >0) .AND. (fn <0) .AND. (fp >0) ) nutlim = n
       ! CASE 7
       if ( (fc >0) .AND. (fn >0) .AND. (fp <0) ) nutlim = p
       ! CASE 8
       if ( (fc >0) .AND. (fn >0) .AND. (fp >0) ) then
          if ( pq_p<pq_n .OR. abs(pq_p-pq_n)<p_small ) then
             if ( pu_p < qpc ) nutlim = p
          else
             if ( pu_n < qnc ) nutlim = n
          endif
       endif

    end function nutlim

  end module mem_globalfun
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
