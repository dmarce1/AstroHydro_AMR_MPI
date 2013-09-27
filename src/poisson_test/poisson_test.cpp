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
#ifdef USE_FMM
	for (int k = BW; k < GNX - BW; k++) {
		for (int j = BW; j < GNX - BW; j++) {
			for (int i = BW; i < GNX - BW; i++) {
#else
	for (int k = 0; k < PNX; k++) {
		for (int j = 0; j < PNX; j++) {
			for (int i = 0; i < PNX; i++) {
#endif
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
#ifdef USE_FMM
				(*this)(i, j, k) = Vector<Real, STATE_NF>(0.0);
				(*this)(i, j, k).set_rho(a / 4.0 / M_PI);
#else
				set_source(i, j, k, a);
#endif
			}
		}
	}

}

Real PoissonTest::get_output_point(int i, int j, int k, int l) const {
	switch (l) {
#ifdef USE_FMM
	case 0:
	return (*this)(i, j, k).rho();
	case 1:
	return get_phi(i, j, k);
	case 2:
	return gx(i, j, k);
	case 3:
	return gy(i, j, k);
	case 4:
	return gz(i, j, k);
#else
	case 0:
		return get_source(i, j, k);
	case 1:
		return get_phi(i, j, k);
	case 2:
		return get_fx(i, j, k);
	case 3:
		return get_fy(i, j, k);
	case 4:
		return get_fz(i, j, k);
#endif
	case 5:
		return analytic_phi(xc(i), yc(j), zc(k));
	case 6:
		return analytic_force(xc(i), yc(j), zc(k))[0];
	case 7:
		return analytic_force(xc(i), yc(j), zc(k))[1];
	case 8:
		return analytic_force(xc(i), yc(j), zc(k))[2];
	}
	return 0.0;
}

int PoissonTest::nvar_output() const {
	return 9;
}

const char* PoissonTest::output_field_names(int i) const {
	switch (i) {
	case 0:
		return "rho";
	case 1:
		return "phi_n";
	case 2:
		return "fx_n";
	case 3:
		return "fy_n";
	case 4:
		return "fz_n";
	case 5:
		return "phi_a";
	case 6:
		return "fx_a";
	case 7:
		return "fy_a";
	case 8:
		return "fz_a";
	}
	return "dummy";
}

Real PoissonTest::analytic_phi(Real x, Real y, Real z) const {
	Real r = sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
#ifdef USE_FMM
	Real a = mass0;
#else
	Real a = mass0 / 4.0 / M_PI;
#endif
	if (r > r0) {
		return -a / r;
	} else {
		return 0.5 * a * r * r / (r0 * r0 * r0) - 1.5 * a / r0;
	}
}

_3Vec PoissonTest::analytic_force(Real x, Real y, Real z) const {
	Real r = sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
	_3Vec force;
	Real d;
#ifdef USE_FMM
	Real a = mass0;
#else
	Real a = mass0 / 4.0 / M_PI;
#endif
	if (r > r0) {
		d = -a / (r * r * r);
	} else {
		d = -a / (r0 * r0 * r0);
	}
	force[0] = (x - x0) * d;
	force[1] = (y - y0) * d;
	force[2] = (z - z0) * d;
	return force;
}

Vector<Real, 3> PoissonTest::momentum_error() {
	Vector<Real, 3> sum = 0.0;
	Real dv = pow(get_dx(), 3.0);
#ifdef USE_FMM
	for (int k = BW; k < GNX - BW; k++) {
		for (int j = BW; j < GNX - BW; j++) {
			for (int i = BW; i < GNX - BW; i++) {
				if (!this->zone_is_refined(i, j, k)) {
					sum[0] += gx(i, j, k) * (*this)(i, j, k).rho() * dv;
					sum[1] += gy(i, j, k) * (*this)(i, j, k).rho() * dv;
					sum[2] += gz(i, j, k) * (*this)(i, j, k).rho() * dv;
#else
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				if (!this->poisson_zone_is_refined(i, j, k)) {
					sum[0] += 0.5 * (get_fx(i, j, k) + get_fx(i + 1, j, k)) * get_source(i, j, k) * dv;
					sum[1] += 0.5 * (get_fy(i, j, k) + get_fy(i, j + 1, k)) * get_source(i, j, k) * dv;
					sum[2] += 0.5 * (get_fz(i, j, k) + get_fz(i, j, k + 1)) * get_source(i, j, k) * dv;
#endif
				}
			}
		}
	}
	return sum;
}

Real PoissonTest::phi_error() {
	Real sum = 0.0;
	Real dv = pow(get_dx(), 3.0);
#ifdef USE_FMM
	for (int k = BW; k < GNX - BW; k++) {
		for (int j = BW; j < GNX - BW; j++) {
			for (int i = BW; i < GNX - BW; i++) {
				if (!this->zone_is_refined(i, j, k)) {
					sum += pow(get_phi(i, j, k) - analytic_phi(xc(i), yc(j), zc(k)), 2) * dv;
#else
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				if (!this->poisson_zone_is_refined(i, j, k)) {
					sum += pow(get_phi(i, j, k) - analytic_phi(Poisson::xc(i), Poisson::yc(j), Poisson::zc(k)), 2) * dv;
#endif
				}
			}
		}
	}
	return sum;
}

_3Vec PoissonTest::force_error() {
	_3Vec sum = 0.0;
	Real dv = pow(get_dx(), 3.0);
#ifdef USE_FMM
	for (int k = BW; k < GNX - BW; k++) {
		for (int j = BW; j < GNX - BW; j++) {
			for (int i = BW; i < GNX - BW; i++) {
				if (!this->zone_is_refined(i, j, k)) {
					_3Vec tmp = analytic_force(xc(i), yc(j), zc(k));
					sum[0] += pow(gx(i, j, k) - tmp[0], 2) * dv;
					sum[1] += pow(gy(i, j, k) - tmp[1], 2) * dv;
					sum[2] += pow(gz(i, j, k) - tmp[2], 2) * dv;
#else
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				if (!this->poisson_zone_is_refined(i, j, k)) {
					_3Vec tmp = analytic_force(xc(i), yc(j), zc(k));
					sum[0] += pow((get_fx(i, j, k) + get_fx(i + 1, j, k)) / 2.0 - tmp[0], 2) * dv;
					sum[1] += pow((get_fy(i, j, k) + get_fy(i, j + 1, k)) / 2.0 - tmp[1], 2) * dv;
					sum[2] += pow((get_fz(i, j, k) + get_fz(i, j, k + 1)) / 2.0 - tmp[2], 2) * dv;
#endif
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
