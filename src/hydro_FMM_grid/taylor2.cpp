#include "taylor2.h"
#include <stdlib.h>

#ifdef USE_FMM
static int map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };

Taylor2::Taylor2() {
// TODO Auto-generated constructor stub

}

Real Taylor2::operator ()() const {
	return l0;
}

Real Taylor2::operator ()(int i) const {
	return l1[i];
}

Real Taylor2::operator ()(int i, int j) const {
	return l2[map2[i][j]];
}

Real& Taylor2::operator ()() {
	return l0;
}

Real& Taylor2::operator ()(int i) {
	return l1[i];
}

Real& Taylor2::operator ()(int i, int j) {
	return l2[map2[i][j]];
}

void Taylor2::zero() {
	l0 = 0.0;
	for (int i = 0; i < 3; i++) {
		l1[i] = 0.0;
	}
	for (int i = 0; i < 6; i++) {
		l2[i] = 0.0;
	}
}

Taylor2& Taylor2::operator +=(const Taylor2& t) {
	l0 += t.l0;
	for (int i = 0; i < 3; i++) {
		l1[i] += t.l1[i];
	}
	for (int i = 0; i < 6; i++) {
		l2[i] += t.l2[i];
	}
	return *this;
}

Taylor2& Taylor2::operator =(const Taylor2& t) {
	l0 = t.l0;
	for (int i = 0; i < 3; i++) {
		l1[i] = t.l1[i];
	}
	for (int i = 0; i < 6; i++) {
		l2[i] = t.l2[i];
	}
	return *this;
}

Taylor2::~Taylor2() {
// TODO Auto-generated destructor stub
}

#endif

