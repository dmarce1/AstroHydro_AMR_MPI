#include "../defs.h"

#ifdef EULER_STATE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../oct_node/oct_face.h"
#include "state.h"
const Real EulerState::gamma = EULER_GAMMA;
const Real EulerState::rho_floor = 1.0e-13;
const Real EulerState::ei_floor = 1.0e-20;
const int EulerState::d_index = 0;
const int EulerState::sx_index = 1;
const int EulerState::sy_index = 2;
const int EulerState::sz_index = 3;
const int EulerState::et_index = 4;
const int EulerState::tau_index = 5;
#ifdef USE_RADIATION
const int EulerState::er_index = 6;
const int EulerState::frac_index = 7;
#else
const int EulerState::frac_index = 6;
#endif

Vector<Real, STATE_NF> EulerState::source(const _3Vec& X, Real t) const {
	Vector<Real, STATE_NF> s = 0.0;
	return s;
}

const char* EulerState::field_name(int i) {
	assert(i >= 0);assert(i < STATE_NF);
	switch (i) {
	case d_index:
		return "d";
	case et_index:
		return "et";
	case sx_index:
		return "sx";
	case sy_index:
		return "sy";
	case sz_index:
		return "sz";
	default:
		//case tau_index:
		return "tau";
	}
}

void EulerState::floor(const _3Vec& X) {
	(*this)[d_index] = max((*this)[d_index], rho_floor);
	const Real ei0 = ei();
	if (ei0 > 0.1 * et()) {
		(*this)[tau_index] = pow(max(ei0, ei_floor), 1.0 / gamma);
	}
	(*this)[tau_index] = max(pow(ei_floor, 1.0 / gamma), (*this)[tau_index]);
}

EulerState::EulerState() :
		Vector<Real, STATE_NF>() {
}

EulerState::EulerState(const Vector<Real, STATE_NF>& v) :
		Vector<Real, STATE_NF>(v) {
	return;
}

Real EulerState::rho() const {
	return (*this)[d_index];
}

Real EulerState::sx() const {
	return (*this)[sx_index];
}

Real EulerState::sy() const {
	return (*this)[sy_index];
}

Real EulerState::sz() const {
	return (*this)[sz_index];
}

Real EulerState::et() const {
	return (*this)[et_index];
}

Real EulerState::ei() const {
	const Real ek = 0.5 * (sx() * sx() + sy() * sy() + sz() * sz()) / rho();
	const Real ei0 = max(et() - ek, 0.0);
	Real this_ei;
	if (ei0 > 0.001 * et()) {
		this_ei = ei0;
	} else {
		this_ei = pow(tau(), gamma);
	}
	return this_ei;
}

Real EulerState::tau() const {
	return (*this)[tau_index];
}

Real EulerState::pg() const {
	return (gamma - 1.0) * ei();
}

Real EulerState::max_abs_x_eigen(const _3Vec& X, Real rad_f) const {
#ifdef USE_RADIATION
	const Real cs2 = (gamma * pg() + rad_f * (4.0 / 9.0) * er()) / rho();
#else
	const Real cs2 = gamma * pg() / rho();
#endif
	return fabs(sx() / rho()) + sqrt(cs2);
}

Real EulerState::max_abs_y_eigen(const _3Vec& X, Real rad_f) const {
#ifdef USE_RADIATION
	const Real cs2 = (gamma * pg() + rad_f * (4.0 / 9.0) * er()) / rho();
#else
	const Real cs2 = gamma * pg() / rho();
#endif
	return fabs(sy() / rho()) + sqrt(cs2);
}

Real EulerState::max_abs_z_eigen(const _3Vec& X, Real rad_f) const {
#ifdef USE_RADIATION
	const Real cs2 = (gamma * pg() + rad_f * (4.0 / 9.0) * er()) / rho();
#else
	const Real cs2 = gamma * pg() / rho();
#endif
	return fabs(sz() / rho()) + sqrt(cs2);
}

Vector<Real, STATE_NF> EulerState::x_vector_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	flux = *this * (sx() / rho());
	flux[et_index] += pg() * (sx() / rho());
	return flux;
}

Vector<Real, STATE_NF> EulerState::y_vector_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	flux = *this * (sy() / rho());
	flux[et_index] += pg() * (sy() / rho());
	return flux;
}

Vector<Real, STATE_NF> EulerState::z_vector_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	flux = *this * (sz() / rho());
	flux[et_index] += pg() * (sz() / rho());
	return flux;
}

void EulerState::set_rho(Real a) {
	(*this)[d_index] = a;
}

void EulerState::set_sx(Real a) {
	(*this)[sx_index] = a;
}

void EulerState::set_sy(Real a) {
	(*this)[sy_index] = a;
}

void EulerState::set_sz(Real a) {
	(*this)[sz_index] = a;
}

void EulerState::set_et(Real a) {
	(*this)[et_index] = a;
}

void EulerState::set_tau(Real a) {
	(*this)[tau_index] = a;
}

void EulerState::to_prim(const _3Vec& X) {
	for (int i = 1; i < STATE_NF; i++) {
		(*this)[i] /= rho();
	}
}

void EulerState::from_prim(const _3Vec& X) {
	for (int i = 1; i < STATE_NF; i++) {
		(*this)[i] *= rho();
	}
}

void EulerState::to_con(const _3Vec& X) {
	return;
}

void EulerState::from_con(const _3Vec& X) {
	return;
}

#endif
