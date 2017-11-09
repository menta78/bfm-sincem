!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: PROFTS
!
!DESCRIPTION    
!
! !INTERFACE
!
      Subroutine PROFTS(FF,WFSURF,WFBOT,SWRAD,FSURF,NBC,DT2,NTP,UMOL)
!
! USES:
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use global_mem,ONLY: RLEN, ZERO,ONE
!
      use POM,ONLY: H,KB,A,C,KH,DZ,DZZ,VH,VHP,Z,ilong
!
!-------------------------------------------------------------------------!
!
!BOC
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      IMPLICIT NONE
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  Scalar Arguments
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!     -----TWICE THE TIME STEP
!
      REAL(RLEN)                          :: DT2
!
!     -----SURFACE TEMPERATURE / SALINITY / TRACER-----
!
      REAL(RLEN)                          :: FSURF
!   
!     -----SURFACE INCIDENT SHORT WAVE RADIATION-----
!
      REAL(RLEN)                          :: SWRAD
!
!      -----SURFACE/bottom HEAT FLUX LOSS TERM OR SALINITY / TRACER FLUX-----
!
      REAL(RLEN)                          :: WFSURF, WFBOT
!
!     -----FLAG FOR BOUNDARY CONDITION DEFINITION-----
!
!     ************************************************
!     ************************************************
!     **                                            **
!     ** NBC=1: SURF. B.C. IS WFSURF+SWRAD. NO      **
!     **        RADIATIVE PENETRATION.              **
!     ** NBC=2; SURF. B.C. IS WFSURF. SWRAD         **
!     **        PENETRATES WATER COLUMN             **
!     ** NBC=3; SURF. B.C. IS TSURF. NO SWRAD       **
!     **        RADIATIVE PENETRATION               **
!     ** NBC=4; SURF. B.C. IS TSURF. SWRAD          **
!     **        PENETRATES WATER COLUMN             **
!     **                                            **
!     ** NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE   **
!     ** NEGATIVE VALUES WHEN WATER COLUMN IS       **
!     ** WARMING.                                   **
!     **                                            **
!     ************************************************
!     ************************************************
!
      INTEGER(ilong)                      :: NBC
!
!
!     -----FLAG FOR JERLOV WATER TYPE CHOICE-----
!
!     ************************************************
!     ************************************************
!     **                                            **
!     ** JERLOV WATER TYPE CHOICE IS RELEVANT ONLY  **
!     ** WHEN NBC = 2 OR NBC = 4                    **
!     **                                            **
!     ************************************************
!     ************************************************
!
      INTEGER(ilong)                      :: NTP
!
!     -----BACKGROUND DIFFUSIVITY------
!
      REAL(RLEN)                          :: UMOL
!
!     -----TEMPERATURE/SALINITY/TRACER-----
!
      REAL(RLEN)                          :: FF(KB)
!
!     -----COUNTERS-----
!
      INTEGER(ilong)                      :: K,KI
!     
!     -----SW PROFILE----
!
      REAL(RLEN)                          :: RAD(KB)
!
!     -----IRRADIANCE PARAMETERS AFTER PAULSON & SIMPSON JPO 1977, 952-956-----
!
      REAL(RLEN)                          :: R(5), AD1(5), AD2(5)
!
!     -----JERLOV WATER TYPES-----
!
!     NTP        = 1      2        3      4      5
!    JERLOV TYPE = I      IA       IB     II     III
!
     DATA R   /      .58,  .8,       .67,   .77,  .78/
     DATA AD1 /      .35,  .2,      1.00,  1.50, 1.40/
     DATA AD2 /    23.00, 5.88235, 17.00, 14.00, 7.90/
!
!     -----INTRINSIC FUNCTION-----
!
      INTRINSIC EXP
!
!     -----START COMPUTATION OF VERTICAL PROFILE----
!
      DO K = 2,KB - 1
!
         A(K-1) = -DT2 * (KH(K)+UMOL)/ (DZ(K-1) * DZZ(K-1) * H * H)
         C(K)   = -DT2 * (KH(K)+UMOL)/ (DZ(K)   * DZZ(K-1) * H * H)
!
      END DO
!
      RAD(:)=ZERO
!
      select case (NBC)
!
         case (2, 4)
!
!        ***********************************************************
!        ***********************************************************
!        **                                                       **
!        ** PENETRATIVE RADIATION CALCULATION.                    **
!        ** AT THE BOTTOM ANY                                     **
!        ** UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER         **
!        **                                                       **
!        ***********************************************************
!        ***********************************************************
!
         RAD(:) = SWRAD * ((    R(NTP) * EXP(Z(:)*H/AD1(NTP))               &  
                        +  (ONE-R(NTP) * EXP(Z(:)*H/AD2(NTP)))))
!        
         RAD(KB)=ZERO
!
      end select

      select case (NBC)

         case (1) 
!
              VH(1)  =  A(1)/(A(1)-ONE)
              VHP(1) = -DT2*(WFSURF+SWRAD)/(-DZ(1)*H) - FF(1)
              VHP(1) =  VHP(1)/(A(1)-ONE)
!
              RAD(:) = ZERO
!
         case (2)
!
              VH(1)  = A(1)/(A(1)-ONE)
              VHP(1) = DT2*(WFSURF+RAD(1)-RAD(2))/(DZ(1)*H) -FF(1)
              VHP(1) = VHP(1)/(A(1)-ONE)
!
         case (3, 4)

              VH(1)  = ZERO
              VHP(1) = FSURF
!
      end select
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  The following section solves the equation                      
!  DT2*(KH*FF')' -FF = FB                                        
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      DO K = 2,KB - 2
!
           VHP(K) = ONE/(A(K)+C(K)*(ONE-VH(K-1))-ONE)
           VH(K)  = A(K)*VHP(K)
           VHP(K) = (C(K)*VHP(K-1) - FF(K))*VHP(K)
!
      END DO
!
      WFBOT=ZERO  ! Giulia
      FF(KB-1) = (C(KB-1)*VHP(KB-2)-FF(KB-1)+(WFBOT*DT2/(DZ(KB-1)*H)))/ &
                 (C(KB-1)* (ONE-VH(KB-2))-ONE)
!
      DO K = 2,KB - 1
!
           KI = KB - K
           FF(KI) = VH(KI)*FF(KI+1) + VHP(KI)
!
      END DO
!
!     -----ZEROING-----
!
      VH(:)  = ZERO
      VHP(:) = ZERO
      A(:)   = ZERO
      C(:)   = ZERO
!
      RETURN
!
     end subroutine PROFTS

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
