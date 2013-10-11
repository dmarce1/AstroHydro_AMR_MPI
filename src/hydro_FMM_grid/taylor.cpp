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

