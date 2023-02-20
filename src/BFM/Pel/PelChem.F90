!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PelChem
!
! DESCRIPTION
!       This process describes the additional dynamics of dissolved
!       compounds in the watercolumn. Parameterized processes are:
!       - nitrification / denitrification
!       - reoxidation of reduction equivalents
!       - dissolution of biogenic silica
!       This function also calls the carbonate system dynamics (INCLUDE_PELCO2)
!       and iron dynamics (INCLUDE_PELFE) if activated
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
#include "DEBUG.h"
#include "INCLUDE.h"
!
! INTERFACE
  subroutine PelChemDynamics
!
! USES
  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: N4n, N3n, O2o, O4n, N6r, R6s, N5s, P1s
  use mem,  ONLY: ppN4n, ppN3n, ppO2o, ppO4n, ppN6r, ppR6s, ppN5s, &
    flN3O4n, ETW, flPTN6r, NO_BOXES, iiBen, iiPel, flN4N3n, flux_vector
#ifdef INCLUDE_PELFE
  use mem,  ONLY: N7f, R6f, R1f, ppN7f, ppR6f, ppR1f, fscavN7f, R6c
#endif
#endif
  use mem_Param,  ONLY: p_qon_nitri, p_qro, p_qon_dentri, p_small
  use mem_PelChem
  use mem_globalfun,   ONLY: MM, eTq, insw, eTa
  use bfm_error_msg,   ONLY: bfm_error

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, save :: first =0
  integer       :: AllocStatus, DeallocStatus
  real(RLEN),allocatable,save,dimension(:) :: fN4N3n,fN6O2r,eo,     &
                                              er,osat,rPAo,fR6N5s
#ifdef INCLUDE_PELFE
  real(RLEN),allocatable,save,dimension(:) :: fR1N7f, fR6N7f
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     ALLOCATE ( fN6O2r(NO_BOXES),   eo(NO_BOXES),   er(NO_BOXES), &
        &       rPAo(NO_BOXES), fR6N5s(NO_BOXES), osat(NO_BOXES), &    
        &       fN4N3n(NO_BOXES),                                 &
#ifdef INCLUDE_PELFE
        &       fR1N7f(NO_BOXES), fR6N7f(NO_BOXES),               &
#endif 
        &      STAT = AllocStatus )
     IF( AllocStatus /= 0 ) call bfm_error('PelChem','Error allocating arrays')
     first=1
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Regulating factors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  eo  =   MM(  max(p_small,O2o(:)),  p_clO2o)
  er  =   MM(  N6r(:),  p_clN6r)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Nitrification in the water
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  flN4N3n(:) =  max(ZERO,p_sN4N3* N4n(:)* eTq(  ETW(:),  p_q10N4N3)* eo)
  call flux_vector( iiPel, ppN4n,ppN3n, flN4N3n(:) )
  call flux_vector( iiPel, ppO2o,ppO2o,-( flN4N3n(:)* p_qon_nitri) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Denitrification in the water
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  rPAo  =   flPTN6r(:)/ p_qro
  flN3O4n(:) = max(ZERO,p_sN3O4n* eTq( ETW(:), p_q10N4N3)* er* rPAo/ p_rPAo* &
               N3n(:))
  call flux_vector( iiPel, ppN3n,ppO4n, flN3O4n(:) )
  call flux_vector( iiPel, ppN6r,ppN6r,-( p_qro* flN3O4n(:)* p_qon_dentri* &
  insw( -( O2o(:)- N6r(:)/ p_qro))) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Reoxidation of reduction equivalents
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  fN6O2r  =   p_rOS* N6r(:)* eo
  call flux_vector( iiPel, ppN6r,ppN6r,-( fN6O2r) )
  call flux_vector( iiPel, ppO2o,ppO2o,-( fN6O2r/ p_qro) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Dissolution of biogenic silicate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  fR6N5s  =   p_sR6N5* eTa(  ETW(:),  p_aeR6N5)* R6s(:)
  call flux_vector( iiPel, ppR6s,ppN5s, fR6N5s )

#ifdef INCLUDE_PELFE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Dissolved Iron Chemistry (dissolution and scavenging)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Linear regeneration of bioavailable iron
  fR1N7f(:)  =  p_sR1N7* eTq(  ETW(:),  p_q10R6N7)* R1f(:)
  call flux_vector( iiPel, ppR1f, ppN7f, fR1N7f(:) )

  fR6N7f(:)  =  p_sR6N7* eTq(  ETW(:),  p_q10R6N7)* R6f(:)
  call flux_vector( iiPel, ppR6f, ppN7f, fR6N7f(:) )

  ! Scavenging of free dissolved Iron
  ! Adsorption onto Particulate Organic Matter (Eq.11 in Parekh et al.,2015 - GBC )
  fscavN7f(:) = max ( ZERO, p_scavOrg * N7f * (R6c(:)**0.58_RLEN) )
  
  ! Inorganic component ( Linear relaxation to the Iron Ligand concentration )
  fscavN7f(:) = fscavN7f(:) + max(ZERO,p_scavIng*(N7f-p_N7fLigand))
  call flux_vector( iiPel, ppN7f, ppN7f, -fscavN7f(:) )

#endif

#ifdef INCLUDE_PELCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Carbonate chemistry
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call PelagicCSYS()
#endif

  end subroutine PelChemDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
