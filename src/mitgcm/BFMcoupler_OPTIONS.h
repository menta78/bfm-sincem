! file BFMcoupler_OPTIONS.h

#ifndef BFMcoupler_OPTIONS_H
#define BFMcoupler_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_BFMCOUPLER
C     Package-specific Options & Macros go here
c  eventual specific options for the BFMcoupler go here

#define USE_QSW
#undef READ_PAR
#define USE_SINK
c #undef USE_SINK
#define USE_SHADE
c #undef READ_xESP
#define READ_xESP
c #define BFMcoupler_DEBUG
#undef BFMcoupler_DEBUG

#endif /* ALLOW_BFMCOUPLER */
#endif /* BFMcoupler_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
