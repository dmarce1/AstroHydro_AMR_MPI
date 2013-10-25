#ifndef OPTIONS_SCF_
#define OPTIONS_SCF_
#include "../real.h"

#define GCON 1.0
#define DIAGNOSTIC_FREQ 1
#define INX                 12
#define GRID_DIM            (1.0)
#define MINMOD_THETA        1.3
#define PPM
#define BW 					(3)
#define EULER_GAMMA         (5.0/3.0)
#define GNX 				(INX+2*BW)
#define PNX 				(INX+2)
#define TIME_MAX           	(1.0e+99)
#define OUTPUT_TIME_FREQ   	(0.025/100.0)
#define MAXDTINC            (1.25)
#define GRID_CFL_FACTOR     (1.0/3.0)
#define MAXINITDT           (1.0e-5)
#define ZTWD
#define DRIVING			    0.01
#define DRIVING_TIME		0.0
//#define NGRID_LIMIT         1000

#endif
