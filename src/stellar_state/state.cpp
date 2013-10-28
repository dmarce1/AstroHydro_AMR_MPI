#include "../defs.h"

#ifdef STELLAR_STATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../oct_node/oct_face.h"
#include "state.h"
const Real StellarState::gamma = 5.0 / 3.0;
Real StellarState::rho_floor = 1.0e-13;
Real StellarState::ei_floor = 1.0e-20;
const int StellarState::d_index = 0;
const int StellarState::sx_index = 1;
const int StellarState::sy_index = 2;
const int StellarState::sz_index = 3;
const int StellarState::et_index = 4;
const int StellarState::tau_index = 5;
const int StellarState::frac_index = 6;
bool StellarState::cylindrical = true;
_3Vec StellarState::drift_vel = 0.0;

Vector<Real, STATE_NF> StellarState::x_scalar_flux_coeff(const _3Vec& X) {
	Vector<Real, STATE_NF> c = 0.0;
	if (cylindrical) {
		c[sx_index] = X[0] / sqrt(X[0] * X[0] + X[1] * X[1]);
		c[sy_index] = 0.0;
	} else {
		c[sx_index] = 0.0;
	}
	return c;
}

Vector<Real, STATE_NF> StellarState::y_scalar_flux_coeff(const _3Vec& X) {
	Vector<Real, STATE_NF> c = 0.0;
	if (cylindrical) {
		c[sx_index] = X[1] / sqrt(X[0] * X[0] + X[1] * X[1]);
		c[sy_index] = 0.0;
	} else {
		c[sy_index] = 0.0;
	}
	return c;
}

Vector<Real, STATE_NF> StellarState::z_scalar_flux_coeff(const _3Vec& X) {
	Vector<Real, STATE_NF> c = 0.0;
	c[sz_index] = 0.0;
	return c;
}

Real StellarState::scalar_flux(const _3Vec& X) const {
	return pg(X);
}

Vector<Real, STATE_NF> StellarState::x_vector_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	const Real v = vx(X);
	flux = (*this) * (v - drift_vel[0]);
	if (State::cylindrical) {
		flux[sy_index] -= X[1] * pg(X);
	} else {
		flux[sx_index] += pg(X);
	}
	flux[et_index] += v * pg(X);
	return flux;
}

Vector<Real, STATE_NF> StellarState::y_vector_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	const Real v = vy(X);
	flux = (*this) * (v - drift_vel[1]);
	if (State::cylindrical) {
		flux[sy_index] += X[0] * pg(X);
	} else {
		flux[sy_index] += pg(X);
	}
	flux[et_index] += v * pg(X);
	return flux;
}

Vector<Real, STATE_NF> StellarState::z_vector_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	const Real v = vz();
	flux = (*this) * (v - drift_vel[2]);
	flux[sz_index] += pg(X);
	flux[et_index] += v * pg(X);
	return flux;
}

Vector<Real, STATE_NF> StellarState::source(const _3Vec& X, Real t) const {
	Vector<Real, STATE_NF> s = 0.0;
	if (!cylindrical) {
		s[sx_index] += omega * (*this)[sy_index];
		s[sy_index] -= omega * (*this)[sx_index];
	} else {
		Real R = sqrt(X[0] * X[0] + X[1] * X[1]);
		s[sx_index] += pow((*this)[sy_index], 2) / (rho() * pow(R, 3));
	}
#ifdef DRIVING
	const Real period = 2.0 * M_PI / omega;
	if (t < DRIVING_TIME * period) {
		if (cylindrical) {
			s[sy_index] -= (*this)[sy_index] * DRIVING / period;
		} else {
			printf("Error in state - cartesian not yet supported for drivingt\n");
			abort();
		}
	}
#endif
	return s;
}

const char* StellarState::field_name(int i) {
	assert(i >= 0);assert(i < STATE_NF);
	switch (i) {
	case State::d_index:
		return "d";
	case et_index:
		return "et";
	case sx_index:
		if (cylindrical) {
			return "sr";
		} else {
			return "sx";
		}
	case sy_index:
		if (cylindrical) {
			return "lz";
		} else {
			return "sy";
		}
	case tau_index:
		return "tau";
	case sz_index:
		return "sz";
	case pot_index:
		return "pot";
	default:
		static char* str = NULL;
		if (str != NULL) {
			free(str);
		}
		asprintf(&str, "frac%i", i - frac_index + 1);
		return str;
	}
}

void StellarState::floor(const _3Vec& X) {
	Real rho1, rho2;
	Real de = pot();
	(*this)[d_index] = max((*this)[d_index], rho_floor);
	const Real ei0 = ei(X);
	if (ei0 > 0.1 * et()) {
		(*this)[tau_index] = pow(max(ei0, ei_floor), 1.0 / gamma);
	}
	(*this)[tau_index] = max(pow(ei_floor, 1.0 / gamma), (*this)[tau_index]);
	if (NFRAC > 1) {
		Real tot;
		for (int i = 0; i < NFRAC; i++) {
			(*this)[frac_index + i] /= rho();
		}
		tot = 0.0;
		for (int i = 0; i < NFRAC; i++) {
			(*this)[frac_index + i] = min(1.0, (*this)[frac_index + i]);
			(*this)[frac_index + i] = max(0.0, (*this)[frac_index + i]);
			tot += (*this)[frac_index + i];
		}
		if (tot <= 0.0) {
			for (int i = 0; i < NFRAC; i++) {
				(*this)[frac_index + i] = (1.0 / Real(NFRAC)) * rho();
			}
		} else {
			for (int i = 0; i < NFRAC; i++) {
				(*this)[frac_index + i] *= rho() / tot;
			}
		}
	}
	de -= pot();
	(*this)[et_index] += de;
}

StellarState::StellarState() :
		Vector<Real, STATE_NF>() {
}

StellarState::StellarState(const Vector<Real, STATE_NF>& v) :
		Vector<Real, STATE_NF>(v) {
	return;
}

Real StellarState::rho() const {
	return (*this)[d_index];
}

Real StellarState::sx() const {
	return (*this)[sx_index];
}

Real StellarState::sy() const {
	return (*this)[sy_index];
}

Real StellarState::sz() const {
	return (*this)[sz_index];
}

Real StellarState::et() const {
	return (*this)[et_index];
}

Real StellarState::vx(const _3Vec& X) const {
	Real v_x;
	if (cylindrical) {
		Real sr, st, R;
		R = sqrt(X[0] * X[0] + X[1] * X[1]);
		sr = sx();
		st = sy() / R;
		v_x = (X[0] * sr - X[1] * st) / R / rho();
	} else {
		v_x = sx() / rho();
	}
	return v_x + X[1] * omega;
}

Real StellarState::vy(const _3Vec& X) const {
	Real v_y;
	if (cylindrical) {
		Real sr, st, R;
		R = sqrt(X[0] * X[0] + X[1] * X[1]);
		sr = sx();
		st = sy() / R;
		v_y = (X[1] * sr + X[0] * st) / R / rho();
	} else {
		v_y = sy() / rho();
	}
	return v_y - X[0] * omega;
}

Real StellarState::vz() const {
	return sz() / rho();
}

Real StellarState::ek(const _3Vec& X) const {
	const Real x = vx(X);
	const Real y = vy(X);
	const Real z = vz();
	return 0.5 * rho() * (x * x + y * y + z * z);
}

Real StellarState::hd() const {
	const Real x = pow(rho() / PhysicalConstants::B, 1.0 / 3.0);
	Real h;
	if (x < 0.01) {
		h = 4.0 * PhysicalConstants::A / PhysicalConstants::B * x * x;
	} else {
		h = 8.0 * PhysicalConstants::A / PhysicalConstants::B * (sqrt(x * x + 1.0) - 1.0);
	}
	return max(h, 0.0);
}

Real StellarState::ed() const {
	return max(hd() * rho() - pd(), 0.0);
}

Real StellarState::pd() const {
	const Real x = pow(rho() / PhysicalConstants::B, 1.0 / 3.0);
	Real p;
	if (x < 0.01) {
		p = 1.6 * PhysicalConstants::A * pow(x, 5);
	} else {
		p = PhysicalConstants::A * (x * (2.0 * x * x - 3.0) * sqrt(x * x + 1.0) + 3.0 * asinh(x));
	}
	return max(p, 0.0);
}

Real StellarState::ei(const _3Vec& X) const {
	const Real ei0 = max(et() - ek(X) - ed(), ei_floor);
	if (ei0 > 0.001 * et()) {
		return ei0;
	} else {
		return pow(max((*this)[tau_index], 0.0), gamma);
	}
}

Real StellarState::pg(const _3Vec& X) const {
//	return (gamma - 1.0) * ei(X);
	return (gamma - 1.0) * ei(X) + pd();
}

Real StellarState::cs(const _3Vec& X) const {
	Real x, dp_depsilon, dp_drho, cs2;
//	return sqrt(gamma * pg(X) / rho());
	x = pow(rho() / PhysicalConstants::B, 1.0 / 3.0);
	dp_drho = ((8.0 * PhysicalConstants::A) / (3.0 * PhysicalConstants::B)) * x * x / sqrt(x * x + 1.0) + gamma * (gamma - 1.0) * ei(X) / rho();
	dp_depsilon = (gamma - 1.0) * rho();
	cs2 = (pg(X) / (rho() * rho())) * dp_depsilon + dp_drho;
	return sqrt(cs2);
}

Real StellarState::max_abs_x_eigen(const _3Vec& X) const {
	Real s = fabs(vx(X)) + cs(X);
	return s;
}

Real StellarState::max_abs_y_eigen(const _3Vec& X) const {
	Real s = fabs(vy(X)) + cs(X);
	return s;
}

Real StellarState::max_abs_z_eigen(const _3Vec& X) const {
	Real s = fabs(vz()) + cs(X);
	return s;
}

void StellarState::set_rho(Real a, int frac) {
	(*this)[d_index + frac] = a;
}

void StellarState::set_sx(Real a) {
	(*this)[sx_index] = a;
}

void StellarState::set_sy(Real a) {
	(*this)[sy_index] = a;
}

void StellarState::set_sz(Real a) {
	(*this)[sz_index] = a;
}

void StellarState::set_et(Real a) {
	(*this)[et_index] = a;
}

void StellarState::set_tau(Real a) {
	(*this)[tau_index] = a;
}

void StellarState::to_prim(const _3Vec& X) {
//	(*this).to_cartesian(X);
	for (int i = 1; i < STATE_NF; i++) {
		(*this)[i] /= rho();
	}
	if (cylindrical) {
		(*this)[sy_index] -= omega * (X[0] * X[0] + X[1] * X[1]);
	} else {
		(*this)[sx_index] += X[1] * omega;
		(*this)[sy_index] -= X[0] * omega;
	}
}

void StellarState::from_prim(const _3Vec& X) {
	if (cylindrical) {
		(*this)[sy_index] += omega * (X[0] * X[0] + X[1] * X[1]);
	} else {
		(*this)[sx_index] -= X[1] * omega;
		(*this)[sy_index] += X[0] * omega;
	}
	for (int i = 1; i < STATE_NF; i++) {
		(*this)[i] *= rho();
	}
//	(*this).from_cartesian(X);
}

Real StellarState::omega = 0.0;

#endif
