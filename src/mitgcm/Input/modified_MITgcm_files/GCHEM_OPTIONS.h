C $Header: /u/gcmpack/MITgcm/pkg/gchem/GCHEM_OPTIONS.h,v 1.9 2011/12/24 01:04:47 jmc Exp $
C $Name: checkpoint64v $


#ifndef GCHEM_OPTIONS_H
#define GCHEM_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_GCHEM

CBOP
C    !ROUTINE: GCHEM_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C options for biogeochemistry package
CEOP

#undef GCHEM_SEPARATE_FORCING
c undefining gchem_separate_forcing actives BFMcoupler_calc_tendency and add_tendency
c #define GCHEM_SEPARATE_FORCING
c defining gchem_separate_forcing actives GCHEM_FORCING_SEP and there BFMcoupler_calc_tendency is used

#endif /* ALLOW_GCHEM */
#endif /* GCHEM_OPTIONS_H */
