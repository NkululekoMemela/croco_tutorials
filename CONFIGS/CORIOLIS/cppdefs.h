/* This cppdefs.h file is based on the example described in 
** coriolis_cppdfes.mov downloaded from
** http://oces.us/eas-ocean-modeling/ROMS-Tutorial/tutorials.html
*/

# define CORIOLIS

# if defined CORIOLIS

# undef OPENMP
# undef MPI

/* Options for Internal Oscillations Example */

/* Model Physics */
#define SOLVE3D
#define UV_COR

/* Computational grid and initial conditions */
#define ANA_GRID
#define ANA_INITIAL

/* Surface boundary conditions */
#define ANA_SMFLUX
#define ANA_STFLUX

/* Bottom boiundary conditions */
/* #define UV_LDRAG - commented out as no longer a cpp option in CROCO */
#define ANA_BTFLUX

/* Lateral boundary conditions */
#define EW_PERIODIC
#define NS_PERIODIC

/* Other options */
#define FLOATS

# endif

/*
**------------------------------------------------------
** Include other internal CPP definitions
**------------------------------------------------------
*/

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
