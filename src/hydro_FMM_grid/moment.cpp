/*
 * Moment.cpp
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#include "moment.h"


static const int m2_map[3][3] = { { 0, 3, 4 }, { 3, 1, 5 }, { 4, 5, 2 } };
static int map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } }, { { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8, 9 } } };



Moment::Moment() {
	// TODO Auto-generated constructor stub

}

Real Moment::M() const {
	return m;
}

Real Moment::M2(int i, int j) const {
	return m2[m2_map[i][j]];
}


Real& Moment::M() {
	return m;
}

Real& Moment::M2(int i, int j) {
	return m2[m2_map[i][j]];
}

Moment& Moment::operator=(const Moment& cp) {
	m = cp.m;
	is_leaf = cp.is_leaf;
	for( int i = 0; i < 6; i++) {
		m2[i] = cp.m2[i];
	}
	X = cp.X;
	return *this;
}

Moment::~Moment() {
	// TODO Auto-generated destructor stub
}

