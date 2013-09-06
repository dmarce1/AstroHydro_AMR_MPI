#ifndef OPTIONS_SCF_45
#define OPTIONS_SCF_45

#define INX 8

#define GRID_DIM            (1.0)
#define MINMOD_THETA        1.3
#define EULER_GAMMA         (5.0/3.0)
#define PNX 				(INX+2)
#define TIME_MAX           	(1e-2)
#define OUTPUT_TIME_FREQ   	(1.0)
#define MAXDTINC           (1.25)
#define GRID_CFL_FACTOR    0.4
#define MAXINITDT          (1.0e-2)
#define BW 3
#define GNX (INX+2*BW)
#define PPM

#define USE_FMM

#endif
