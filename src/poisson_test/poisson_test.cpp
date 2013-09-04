/*
 * poisson_test.cpp
 *
 *  Created on: Mar 6, 2013
 *      Author: dmarce1
 */

#include "poisson_test.h"

#ifdef POISSON_TEST

Real PoissonTest::mass0;
Real PoissonTest::x0 = -0.5;
Real PoissonTest::y0 = 0.0;
Real PoissonTest::z0 = 0.0;
Real PoissonTest::r0 = 0.1;

PoissonTest* PoissonTest::new_octnode() const {
	return new PoissonTest;
}

void PoissonTest::initialize() {
	Real r1;
	Real x, y, z;
	const Real dm = 1.0 / (0.1 * 0.1 * 0.1 * 4.0 / 3.0 * M_PI);
	int M = 16;
	Real r0, a;
	for (int k = 0; k < PNX; k++) {
		for (int j = 0; j < PNX; j++) {
			for (int i = 0; i < PNX; i++) {
				x = this->xc(i);
				y = this->yc(j);
				z = this->zc(k);
				r0 = max(0.1, get_dx());
				r1 = sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
				if (r1 < r0 - sqrt(3.0 / 4.0) * get_dx()) {
					a = dm;
				} else if (r1 > r0 + sqrt(3.0 / 4.0) * get_dx()) {
					a = 1.0e-20;
				} else {
					a = 0.0;
					for (int i1 = 0; i1 < M; i1++) {
						for (int j1 = 0; j1 < M; j1++) {
							for (int k1 = 0; k1 < M; k1++) {
								x = xc(i) + get_dx() * (Real(i1) + 0.5) / Real(M) - 0.5 * get_dx();
								y = yc(j) + get_dx() * (Real(j1) + 0.5) / Real(M) - 0.5 * get_dx();
								z = zc(k) + get_dx() * (Real(k1) + 0.5) / Real(M) - 0.5 * get_dx();
								r1 = sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
								if (r1 < r0) {
									a += dm / M / M / M;
								} else {
									a += 1.0e-20 / M / M / M;
								}
							}
						}
					}
				}
				set_source(i, j, k, a);
			}
		}
	}

}

Real PoissonTest::get_output_point(int i, int j, int k, int l) const {
	switch (l) {
	case 0:
		return get_source(i, j, k);
	case 1:
		return get_phi(i, j, k);
	case 2:
		return get_residual(i, j, k);
	case 3:
		return get_fx(i, j, k);
	case 4:
		return get_fy(i, j, k);
	case 5:
		return get_fz(i, j, k);
	case 6:
		return analytic_phi(xc(i), yc(j), zc(k));
	}
	return 0.0;
}

int PoissonTest::nvar_output() const {
	return 7;
}

const char* PoissonTest::output_field_names(int i) const {
	switch (i) {
	case 0:
		return "rho";
	case 1:
		return "numerical_phi";
	case 2:
		return "residual";
	case 3:
		return "get_fx";
	case 4:
		return "get_fy";
	case 5:
		return "get_fz";
	case 6:
		return "analytic_phi";
	}
	return "dummy";
}

Real PoissonTest::analytic_phi(Real x, Real y, Real z) const {
	Real r = sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
	Real a = mass0 / 4.0 / M_PI;
	if (r > r0) {
		return -a / r;
	} else {
		return 0.5 * a * r * r / (r0 * r0 * r0) - 1.5 * a / r0;
	}
}

Vector<Real, 3> PoissonTest::momentum_error() {
	Vector<Real, 3> sum = 0.0;
	Real dv = pow(get_dx(), 3.0);
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				if (!this->poisson_zone_is_refined(i, j, k)) {
					sum[0] += 0.5 * (get_fx(i, j, k) + get_fx(i + 1, j, k)) * get_source(i, j, k) * dv;
					sum[1] += 0.5 * (get_fy(i, j, k) + get_fy(i, j + 1, k)) * get_source(i, j, k) * dv;
					sum[2] += 0.5 * (get_fz(i, j, k) + get_fz(i, j, k + 1)) * get_source(i, j, k) * dv;
				}
			}
		}
	}
	return sum;
}

Real PoissonTest::phi_error() {
	Real sum = 0.0;
	Real dv = pow(get_dx(), 3.0);
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				if (!this->poisson_zone_is_refined(i, j, k)) {
					sum += pow(get_phi(i, j, k) - analytic_phi(Poisson::xc(i), Poisson::yc(j), Poisson::zc(k)), 2) * dv;
				}
			}
		}
	}
	return sum;
}

PoissonTest::PoissonTest() {
// TODO Auto-generated constructor stub

}

PoissonTest::~PoissonTest() {
// TODO Auto-generated destructor stub
}

#endif
