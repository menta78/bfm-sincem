! This file is included in BFMcoupler_fields_load.F and BFMcoupler_ini_forcing.F
! File BFMcoupler_LOAD.h


C     BFMcoupler_ldRec     :: time-record currently loaded (in temp arrays *[1])

      COMMON /BFMcoupler_LOAD_I/ BFMcoupler_ldRec
      INTEGER BFMcoupler_ldRec(nSx,nSy)

      COMMON /BFMcouplerLOAD_RS/
     &    AtmosPCO20,AtmosPCO21,AtmosP0,AtmosP1,
     &    AtmosWIND0,AtmosWIND1,
     &    N1p_surfF0,N1p_surfF1,N3n_surfF0,N3n_surfF1,
     &    N5s_surfF0,N5s_surfF1,O3c_surfF0,O3c_surfF1,
     &    O3h_surfF0,O3h_surfF1,N1p_botF0,N1p_botF1,
     &    N3n_botF0,N3n_botF1,N4n_botF0,N4n_botF1,
     &    O2o_botF0,O2o_botF1,O3c_botF0,O3c_botF1,
     &    O3h_botF0,O3h_botF1
#ifdef READ_xESP
     &    ,xESP0,xESP1
#endif
#ifdef READ_PAR
     &    ,spar0,spar1
#endif
      _RS AtmosWIND0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS AtmosWIND1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS AtmosPCO20  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS AtmosPCO21  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS AtmosP0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS AtmosP1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N1p_surfF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N1p_surfF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N3n_surfF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N3n_surfF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N5s_surfF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N5s_surfF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3c_surfF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3c_surfF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3h_surfF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3h_surfF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N1p_botF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N1p_botF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N3n_botF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N3n_botF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N4n_botF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS N4n_botF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O2o_botF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O2o_botF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3c_botF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3c_botF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3h_botF0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS O3h_botF1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef READ_xESP
      _RS xESP0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS xESP1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef READ_PAR
      _RS spar0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS spar1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
