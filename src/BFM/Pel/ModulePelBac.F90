!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: PelBac
!
! DESCRIPTION
!   Module containing the parameters for bacterioplankton and the 
!   initialization and consistency check
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
#include "cppdefs.h"
!
! INTERFACE
  module mem_PelBac
!
! USES
  use global_mem
  use mem,  ONLY: iiPelBacteria,D3STATETYPE,ppR2c,ppR3c

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  integer,private     :: i

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !NAMELIST PelBacteria_parameters
  !-------------------------------------------------------------------------!
  !  PELAGIC BACTERIA
  !
  ! NAME         [UNIT]/KIND            DESCRIPTION
  ! p_version   integer         Switch for bacteria parameterization
  !                              1 : Baretta-Bekker et al. 1995;
  !                                  Vichi et al., 2007
  !                              2 : Vichi et al., 2004
  !                              3 : Polimene et al., 2006
  ! p_q10                        Q10-value (temperature dependency)
  ! p_chdo      [mmol/m3]        Half-saturation constant for O2 limitation
  ! p_sd        [1/d]            Specific mortality rate
  ! p_sd2       [1/d]            Density dependent specific mortality rate
  ! p_suhR1     [1/d]            Specific potential uptake for nutrient-rich DOM
  ! p_sulR1     [1/d]            Specific potential uptake for nutrient-poor DOM
  ! p_suR2      [1/d]            Specific potential uptake for semi-labile DOC
  ! p_suR3      [1/d]            Specific potential uptake for semi-refractory DOC
  ! p_suR6      [1/d]            Specific potential uptake for POM (1/d)
  ! p_sum       [1/d]            Potential specific growth rate
  ! p_pu_ra     [-]              Activity respiration fraction
  ! p_pu_ra_o   [-]              Additional respiration fraction at low O2 conc
  ! p_srs       [1/d]            Specific rest respiration
  ! p_qncPBA    [mmolN/mgC]      Optimal N/C ratio 
  ! p_qpcPBA    [mmolP/mgC]      Optimal P/C ratio 
  ! p_qlnc      [mmolN/mgC]      Minimal N/C ratio 
  ! p_qlpc      [mmolP/mgC]      Minimal P/C ratio 
  ! p_qun       [mmolN/mgC/day]  Membrane affinity for N 
  ! p_qup       [mmolP/mgC/day]  Membrane affinity for N 
  ! p_chn       [mmolN/m3]       Half saturation ammonium conc. for uptake
  ! p_chp       [mmolP/m3]       Half saturation phosphate conc. for uptake
  ! p_ruen      [1/d]            Relaxation timescale for N uptake/remin.
  ! p_ruep      [1/d]            Relaxation timescale for P uptake/remin.
  ! p_rec       [1/d]            Relaxation timescale for semi-labile excretion
  ! p_pu_ea_R3  [-]              Excretion of semi-refractory DOC
  ! p_chuc      [mgC/m3]         Half saturation total carbon uptake from substrate
  ! p_chuc_lim  [-]              Scaling factor of p_chuc to set minimum value
  integer     :: p_version(iiPelBacteria), itrp
  integer, parameter ::       BACT1=1,BACT2=2,BACT3=3
  real(RLEN)  :: p_q10(iiPelBacteria)
  real(RLEN)  :: p_chdo(iiPelBacteria)
  real(RLEN)  :: p_sd(iiPelBacteria)
  real(RLEN)  :: p_sd2(iiPelBacteria)
  real(RLEN)  :: p_suhR1(iiPelBacteria)
  real(RLEN)  :: p_sulR1(iiPelBacteria)
  real(RLEN)  :: p_suR2(iiPelBacteria)
  real(RLEN)  :: p_suR6(iiPelBacteria)
  real(RLEN)  :: p_suR3(iiPelBacteria)
  real(RLEN)  :: p_sum(iiPelBacteria)
  real(RLEN)  :: p_pu_ra(iiPelBacteria)
  real(RLEN)  :: p_pu_ra_o(iiPelBacteria)
  real(RLEN)  :: p_srs(iiPelBacteria)
  real(RLEN)  :: p_qncPBA(iiPelBacteria)
  real(RLEN)  :: p_qpcPBA(iiPelBacteria)
  real(RLEN)  :: p_qlnc(iiPelBacteria)
  real(RLEN)  :: p_qlpc(iiPelBacteria)
  real(RLEN)  :: p_qun(iiPelBacteria)
  real(RLEN)  :: p_qup(iiPelBacteria)
  real(RLEN)  :: p_chn(iiPelBacteria)
  real(RLEN)  :: p_chp(iiPelBacteria)
  real(RLEN)  :: p_ruen(iiPelBacteria)
  real(RLEN)  :: p_ruep(iiPelBacteria)
  real(RLEN)  :: p_rec(iiPelBacteria)
  real(RLEN)  :: p_pu_ea_R3(iiPelBacteria)
  real(RLEN)  :: p_chuc(iiPelBacteria)
  real(RLEN)  :: p_dep(iiPelBacteria)
  real(RLEN)  :: p_dep_exp(iiPelBacteria)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitPelBac

  contains

  subroutine InitPelBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /PelBacteria_parameters/ p_version, p_q10, p_chdo, p_sd, p_sd2, p_suhR1, &
    p_sulR1, p_suR2, p_suR6, p_sum, p_pu_ra, p_pu_ra_o, p_pu_ea_R3, p_srs, &
    p_suR3, p_qpcPBA, p_qlpc, p_qncPBA, p_qlnc, p_qun, p_qup, p_chn, p_chp, &
    p_ruen, p_ruep, p_rec, p_chuc, p_dep, p_dep_exp

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  LEVEL1 "#  Reading PelBac parameters.."
  open(NMLUNIT,file='Pelagic_Ecology.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=PelBacteria_parameters,err=101)
  close(NMLUNIT)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Check consistency of parameters according to the parametrization
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do i=1,iiPelBacteria
     LEVEL1 "#  Checking PelBacteria parameters for group:",i
     LEVEL1 "#   Use formulation (p_version) -> ",p_version
     select case ( p_version(i) )
       case ( BACT3 ) ! Polimene et al. (2006)
         p_sulR1(i) = ZERO
         LEVEL1 "#   forcing p_sulR1=0"
         if (p_pu_ea_R3(i) + p_pu_ra(i) .GT. 0.3_RLEN) then
           LEVEL1 "#  Warning: Bacterial growth efficiency is lower than 0.3!"
           LEVEL1 "#  The release of capsular material is possibly larger than p_pu_ra/4"
         end if
       case ( BACT1 ) ! Vichi et al. 2007
         p_sulR1(i) = ZERO
         p_suR2(i) = ZERO
         p_suR3(i) = ZERO
         LEVEL1 "#   forcing p_sulR1, p_suR2, p_suR3=0"
       case ( BACT2 ) ! Vichi et al. 2004
         p_suR3(i) = ZERO
         LEVEL1 "#   forcing p_suR3=0"
     end select
  end do
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Check across bacterial groups if R2 and R3 must be kept in transport
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  itrp=maxval(p_version)
  select case (itrp)
     case ( 1 )
        D3STATETYPE(ppR3c)=NOTRANSPORT
        D3STATETYPE(ppR2c)=NOTRANSPORT
        LEVEL1 " Disable R2c & R3c transport as no bacterial group use it "
     case ( 4 )
        D3STATETYPE(ppR2c)=NOTRANSPORT
        LEVEL1 " Disable R2c transport as no bacterial group use it "
     case default
        LEVEL1 " Active R2c & R3c transport"
  end select
  LEVEL1 ""
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Write parameter list to the log
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  LEVEL1 "#  Namelist is:"
  if (bfm_lwp) write(LOGUNIT,nml=PelBacteria_parameters)

  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPelBac.f90","Pelagic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitPelBac.f90","PelBacteria_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end subroutine InitPelBac

  end module mem_PelBac

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
