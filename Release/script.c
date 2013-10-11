
#include <stdio.h>
#include <stdlib.h>


int main( void ) {
	double interval = 0.0025;
	char str[1200];
	double x, y;
	for( x = 0.0; x <= 0.5; x += interval ) {
		y = 0.49;
//		for( y = 0.0; y <= 0.95; y += interval ) {
			sprintf( str, "/usr/bin/mpirun -np 4 ./AstroHydro_AMR_MPI 5 %f %f %f\n", x, y, 0.0 );
			system( str );
		}
//	}


}
