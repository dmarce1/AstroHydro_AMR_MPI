#ifndef OPTIONS_____H
#define OPTIONS_____H

//#define ROTATING_DISC
//#define BLAST_WAVE
//#define SOD
//#define POISSON_TEST
#define SINGLE_STAR
//#define BINARY_STAR
//#define FMM_TEST
//#define RADIATION_TEST

#define USE_HYDRO_GRID




#ifdef FMM_TEST
#define STELLAR_STATE
#define HYDRO_GRAV_GRID
#define NFRAC 2
#include "./hydro_FMM_grid/defs.h"
#include "stellar_state/state.h"
class HydroFMMGrid;
typedef HydroFMMGrid ProblemGrid;
#endif

#ifdef BINARY_STAR
#define STELLAR_STATE
#define HYDRO_GRAV_GRID
#define NFRAC 2
#include "./binary_star/defs.h"
#include "stellar_state/state.h"
class BinaryStar;
typedef BinaryStar ProblemGrid;
#endif


#ifdef SINGLE_STAR
#define STELLAR_STATE
#define NFRAC 2
#define HYDRO_GRAV_GRID
#define USE_FMM
#include "./single_star/defs.h"
#include "stellar_state/state.h"
class SingleStar;
typedef SingleStar ProblemGrid;
#endif

#ifdef RADIATION_TEST
#define USE_RADIATION
#define EULER_STATE
#include "./radiation_test/defs.h"
#include "euler_state/state.h"
class RadiationTest;
typedef RadiationTest ProblemGrid;
#endif

#ifdef POISSON_TEST
#undef USE_HYDRO_GRID
#include "./poisson_test/poisson_test_defs.h"
class PoissonTest;
typedef PoissonTest ProblemGrid;
#endif

#ifdef BLAST_WAVE
#define EULER_STATE
#include "./blast_wave/euler_defs.h"
#include "euler_state/state.h"
class GridBlastWave;
typedef GridBlastWave ProblemGrid;
#endif

#ifdef SOD
#define EULER_STATE
#include "./sod/euler_defs.h"
#include "euler_state/state.h"
class GridSod;
typedef GridSod ProblemGrid;
#endif

#ifdef ROTATING_DISC
#undef EULER_STATE
#include "rotating_disc/euler_defs.h"
#include "rotating_disc/state.h"
class RotatingDisc;
typedef RotatingDisc ProblemGrid;
#endif

#endif


