#ifndef OPTIONS_SCF_1
#define OPTIONS_SCF_1
#include "../real.h"

//#define HYDRO_OFF

#define INX                 12
#define GRID_DIM            (0.75e+10)
#define MINMOD_THETA        1.3
#define PPM
#define BW 					(3)

#define EULER_GAMMA         (5.0/3.0)
#define GNX 				(INX+2*BW)
#define PNX 				(INX+2)
#define TIME_MAX           	(60.0)
#define OUTPUT_TIME_FREQ   	(1.0e-1)
#define MAXDTINC            (1.25)
#define GRID_CFL_FACTOR     0.125
#define MAXINITDT           (1.0e-3)
#define CHECKPT_FREQ        32
#define GCON (6.67e-8)
#endif
