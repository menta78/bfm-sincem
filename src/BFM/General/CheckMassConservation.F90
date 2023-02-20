!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! ROUTINE: CheckMassConservation
!
! DESCRIPTION
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
  subroutine CheckMassConservationDynamics
!
! USES
  use global_mem,   ONLY: RLEN,ZERO,ONE,LOGUNIT
  use constants,    ONLY: MW_C, MW_P, MW_N, MW_SI
  use mem
  use mem_Param,    ONLY: CalcBenthicFlag,p_d_tot
  use time, only: bfmtime

  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)  :: s
  integer         :: i,j
  real(RLEN),save :: prevsysc,prevsysn,prevsysp,prevsyss
  real(RLEN),save :: initialc,initialn,initialp,initials
  integer         :: prec
  logical,save    :: first=.TRUE.
  logical,save    :: flag=.FALSE.
  real(RLEN),parameter :: p_prec=1.e-12_RLEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  totpelc(:)=ZERO
  totpelp(:)=ZERO
  totpeln(:)=ZERO
  totpels(:)=ZERO

  ! Pelagic Bacteria
  do i=1, iiPelBacteria
     totpelc(:)=totpelc(:) + PelBacteria(i,iiC)
     totpeln(:)=totpeln(:) + PelBacteria(i,iiC)*qncPBA(:,i)
     totpelp(:)=totpelp(:) + PelBacteria(i,iiC)*qpcPBA(:,i)
  end do
  ! PhytoPlankton
  do i=1, iiPhytoPlankton
     totpelc(:)=totpelc(:) + PhytoPlankton(i,iiC)
     totpeln(:)=totpeln(:) + PhytoPlankton(i,iiC)*qncPPY(:,i)
     totpelp(:)=totpelp(:) + PhytoPlankton(i,iiC)*qpcPPY(:,i)
     totpels(:)=totpels(:) + PhytoPlankton(i,iiC)*qscPPY(:,i)
  end do
  ! MicroZooplankton
  do i=1, iiMicroZooplankton
     totpelc(:)=totpelc(:) + MicroZooplankton(i,iiC)
     totpeln(:)=totpeln(:) + MicroZooplankton(i,iiC)*qncMIZ(:,i)
     totpelp(:)=totpelp(:) + MicroZooplankton(i,iiC)*qpcMIZ(:,i)
  end do
  ! MesoZooPlankton
  do i=1, iiMesoZooPlankton
     totpelc(:)=totpelc(:) + MesoZooPlankton(i,iiC)
     totpeln(:)=totpeln(:) + MesoZooPlankton(i,iiC)*qncMEZ(:,i)
     totpelp(:)=totpelp(:) + MesoZooPlankton(i,iiC)*qpcMEZ(:,i)
   end do
  ! Pelagic Detritus
  do i=1, iiPelDetritus
     if ( ppPelDetritus(i,iiC)/=0) then
        s=PelDetritus(i,iiC)
        totpelc(:)=totpelc(:) + s
     end if
     if ( ppPelDetritus(i,iiN)/=0) then
        s=PelDetritus(i,iiN)
        totpeln(:)=totpeln(:) + s
     end if
     if ( ppPelDetritus(i,iiP)/=0) then
        s=PelDetritus(i,iiP)
        totpelp(:)=totpelp(:) + s
     end if
     if ( ppPelDetritus(i,iiS)/=0) then
        s=PelDetritus(i,iiS)
        totpels(:)=totpels(:) + s
     end if
  end do

#ifdef INCLUDE_PELCO2
  totpelc(:) = totpelc(:)+ O3c(:)
#endif
  ! Convert from default units to g and multiply for the water volume
  totpelc(:) = totpelc(:)*Volume(:)/1000.0_RLEN
  ! Convert from default units to g and multiply for the water volume
  totpeln(:) = (totpeln(:)+ ( N3n(:) + N4n(:) + O4n(:))) &
               *Volume(:)*MW_N/1000.0_RLEN
  totpelp(:) = (totpelp(:)+ N1p(:)) &
               *Volume(:)*MW_P/1000.0_RLEN
  totpels(:) = (totpels(:)+ N5s(:)) &
               *Volume(:)*MW_SI/1000.0_RLEN

  totsysc(:) = sum(totpelc(:))
  totsysn(:) = sum(totpeln(:))
  totsysp(:) = sum(totpelp(:))
  totsyss(:) = sum(totpels(:))

  ! Mass conservation variables
  totbenc(:)  =  ZERO
  totbenp(:)  =  ZERO
  totbenn(:)  =  ZERO
  totbens(:)  =  ZERO
  
  if ( CalcBenthicFlag ) then

#if defined BENTHIC_BIO
     totbenc(:) = ( Y1c(:)+ Y2c(:)+ Y3c(:)+ &
                    Y4c(:)+ Y5c(:)+ H1c(:)+ &
                    H2c(:)+ Q1c(:)+ Q6c(:)+ &
                    Q16c(:)+ Q11c(:) )

     totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ &
                    Y4p(:)+ Y5p(:)+ H1p(:)+ &
                    H2p(:)+ Q1p(:)+ Q6p(:)+ &
                    Q16p(:)+ Q11p(:)+ K1p(:)+ K11p(:))
     totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ &
                    Y4n(:)+ Y5n(:)+ H1n(:)+ &
                    H2n(:)+ Q1n(:)+ Q6n(:)+ &
                    Q16n(:)+ Q11n(:)+ K4n(:)+ K14n(:))
     totbens(:)  =  Q6s(:)+ Q16s(:)

#elif defined BENTHIC_FULL
     totbenc(:) = ( Y1c(:)+ Y2c(:)+ Y3c(:)+ &
                    Y4c(:)+ Y5c(:)+ H1c(:)+ &
                    H2c(:)+ Q1c(:)+ Q6c(:)+ &
                    Q11c(:) )

     totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ &
                    Y4p(:)+ Y5p(:)+ H1p(:)+ &
                    H2p(:)+ Q1p(:)+ Q6p(:)+ &
                    Q11p(:)+ K1p(:)+ K11p(:)+ K21p(:))
     totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ &
                    Y4n(:)+ Y5n(:)+ H1n(:)+ &
                    H2n(:)+ Q1n(:)+ Q6n(:)+ &
                    Q11n(:)+ G4n(:)+ K3n(:)+&
                    K4n(:)+ K14n(:)+ K24n(:))
     totbens(:) =   K5s(:)+ Q6s(:)

#else
      totbenc(:)  =  ( Q1c(:)+ Q11c(:)+ Q6c(:)+ Q16c(:))
      totbenp(:)  =  ( Q1p(:)+ Q11n(:)+ Q6p(:)+ Q16p(:))
      totbenn(:)  =  ( Q1n(:)+ Q11p(:)+ Q6n(:)+ Q16n(:))
      totbens(:)  =    Q6s(:)+ Q16s(:)
#endif

#ifdef INCLUDE_BENCO2
     totbenc(:) = totbenc(:) + G3c(:)
#endif
     ! Convert from default units to g and multiply for the sediment volume
     totbenc(:) = totbenc(:)/1000.0_RLEN*Area2d(:)
     totbenn(:) = totbenn(:)*MW_N/1000.0_RLEN*Area2d(:)
     totbenp(:) = totbenp(:)*MW_P/1000.0_RLEN*Area2d(:)
     totbens(:) = totbens(:)*MW_Si/1000.0_RLEN*Area2d(:)

  endif

  ! Add benthic mass to the total
  totsysc(:) = totsysc(:)+sum(totbenc(:))
  totsysn(:) = totsysn(:)+sum(totbenn(:))
  totsysp(:) = totsysp(:)+sum(totbenp(:))
  totsyss(:) = totsyss(:)+sum(totbens(:))

  ! Store and check previous value
  if (first) then
     write(LOGUNIT,*) "Initializing Mass Conservation"
     first = .FALSE.
     flag  = .FALSE.
     initialc = totsysc(1)
     initialn = totsysn(1)
     initialp = totsysp(1)
     initials = totsyss(1)
  else
     prec = precision(prevsysc)
     write(LOGUNIT,*) ""
     write(LOGUNIT,*) "Check Mass Conservation at step ", bfmtime%stepnow
     if (bfmtime%stepnow ==  bfmtime%step0) & 
     write(LOGUNIT,"(a,i6,a,1D15.8)") "---> Using precision digits ",prec, ", specified precision threshold ", p_prec
     write(LOGUNIT,"(15x,A,15x,A)")  "Current", "Previous" 
     write(LOGUNIT,"(A,2D22.15)") "---> C :",totsysc(1),prevsysc
     write(LOGUNIT,"(A,2D22.15)") "---> N :",totsysn(1),prevsysn
     if (abs(totsysn(1)/initialn-ONE)>p_prec) then
        flag = .TRUE.
        write(LOGUNIT,"(A,1D22.15)") "------> Change in N larger than specified precision : ",totsysn(1)/initialn-ONE
     end if
     write(LOGUNIT,"(A,2D22.15)") "---> P :",totsysp(1),prevsysp
     if (abs(totsysp(1)/initialp-ONE)>p_prec) then
        flag = .TRUE.
        write(LOGUNIT,"(A,1D22.15)") "------> Change in P larger than specified precision : ",totsysp(1)/initialp-ONE
     end if
     write(LOGUNIT,"(A,2D22.15)") "---> Si:",totsyss(1),prevsyss
     if (abs(totsyss(1)/initials-ONE)>p_prec) then
        flag = .TRUE.
        write(LOGUNIT,"(A,1D22.15)") "------> Change in Si larger than specified precision : ",totsyss(1)/initials-ONE
     end if
     if (flag)  then
        write(LOGUNIT,*) ""
        write(LOGUNIT,*) "Check also BFM_General.nml settings!"
        call flush(LOGUNIT)
        stop "Mass conservation violation in BFM! Check log file."
     end if
  end if
  prevsysc = totsysc(1) !first element is sufficient
  prevsysn = totsysn(1)
  prevsysp = totsysp(1)
  prevsyss = totsyss(1)

  end subroutine CheckMassConservationDynamics

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
