/*
 * taylor.cpp
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#include "taylor.h"
#include <stdlib.h>

#ifdef USE_FMM
static int map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
static int map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } }, { { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8, 9 } } };
static int map4[3][3][3][3];

__attribute__((constructor))
static void init_map4(void) {
	int m = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					map4[i][j][k][l] = -1;
				}
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = i; j < 3; j++) {
			for (int k = j; k < 3; k++) {
				for (int l = k; l < 3; l++) {
					map4[i][j][k][l] = m;
					map4[i][j][l][k] = m;
					map4[i][k][j][l] = m;
					map4[i][k][l][j] = m;
					map4[i][l][j][k] = m;
					map4[i][l][k][j] = m;

					map4[j][i][k][l] = m;
					map4[j][i][l][k] = m;
					map4[j][k][i][l] = m;
					map4[j][k][l][i] = m;
					map4[j][l][i][k] = m;
					map4[j][l][k][i] = m;

					map4[k][i][j][l] = m;
					map4[k][i][l][j] = m;
					map4[k][j][i][l] = m;
					map4[k][j][l][i] = m;
					map4[k][l][i][j] = m;
					map4[k][l][j][i] = m;

					map4[l][i][j][k] = m;
					map4[l][i][k][j] = m;
					map4[l][j][i][k] = m;
					map4[l][j][k][i] = m;
					map4[l][k][i][j] = m;
					map4[l][k][j][i] = m;
					m++;
				}
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					if(map4[i][j][k][l] == -1 ) {
						printf( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11\n");
						abort();
					}
				}
			}
		}
	}	if (m != 15) {
		printf("m = %i\n", m);
		exit(1);
	}
}

Taylor::Taylor() {
// TODO Auto-generated constructor stub

}

Real Taylor::operator ()() const {
	return l0;
}

Real Taylor::operator ()(int i) const {
	return l1[i];
}

Real Taylor::operator ()(int i, int j) const {
	return l2[map2[i][j]];
}

Real Taylor::operator ()(int i, int j, int k) const {
	return l3[map3[i][j][k]];
}

Real& Taylor::operator ()() {
	return l0;
}

Real& Taylor::operator ()(int i) {
	return l1[i];
}

Real& Taylor::operator ()(int i, int j) {
	return l2[map2[i][j]];
}

Real& Taylor::operator ()(int i, int j, int k) {
	return l3[map3[i][j][k]];
}


void Taylor::zero() {
	l0 = 0.0;
	for (int i = 0; i < 3; i++) {
		l1[i] = 0.0;
	}
	for (int i = 0; i < 6; i++) {
		l2[i] = 0.0;
	}
	for (int i = 0; i < 10; i++) {
		l3[i] = 0.0;
	}
}

Taylor& Taylor::operator +=(const Taylor& t) {
	l0 += t.l0;
	for (int i = 0; i < 3; i++) {
		l1[i] += t.l1[i];
	}
	for (int i = 0; i < 6; i++) {
		l2[i] += t.l2[i];
	}
	for (int i = 0; i < 10; i++) {
		l3[i] += t.l3[i];
	}
	return *this;
}

Taylor& Taylor::operator =(const Taylor& t) {
	l0 = t.l0;
	for (int i = 0; i < 3; i++) {
		l1[i] = t.l1[i];
	}
	for (int i = 0; i < 6; i++) {
		l2[i] = t.l2[i];
	}
	for (int i = 0; i < 10; i++) {
		l3[i] = t.l3[i];
	}
	return *this;
}

Taylor::~Taylor() {
// TODO Auto-generated destructor stub
}

#endif
