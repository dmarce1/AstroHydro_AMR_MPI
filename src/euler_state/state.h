#ifndef EULER0_STATE_H_
#define EULER0_STATE_H_

#include "../defs.h"
#ifdef EULER_STATE
#include "../real.h"
#include "state.h"
#include "../vector.h"
#include "../oct_node/oct_face.h"

//#define TRACE

#ifdef USE_RADIATION
#define STATE_NF 7
#else
#define STATE_NF 6
#endif

class EulerState: public Vector<Real, STATE_NF> {
public:
	static const Real gamma;
	static const Real ei_floor;
	static const Real rho_floor;
	static const int d_index;
	static const int sx_index;
	static const int sy_index;
	static const int sz_index;
	static const int et_index;
	static const int tau_index;
	static const int frac_index;
#ifdef USE_RADIATION
	static const int er_index;
	void set_er( Real a ) {
		(*this)[er_index] = a;
	}
	Real er() const {
		return (*this)[er_index];
	}
	Real cv() const {
		return 3.0 / 2.0;
	}
#endif
	EulerState();
	EulerState(const Vector<Real, STATE_NF>&);
	~EulerState() {
		return;
	}
	static bool low_order_variable(int i) {
		return false;
	}
	static bool smooth_variable(int i) {
		return false;
	}
	Vector<Real, STATE_NF> source(const _3Vec& X, Real t) const;
	Real max_abs_x_eigen(const _3Vec& X, Real frad = 0.0) const;
	Real max_abs_y_eigen(const _3Vec& X, Real frad = 0.0) const;
	Real max_abs_z_eigen(const _3Vec& X, Real frad = 0.0) const;
	Vector<Real, STATE_NF> x_vector_flux(const _3Vec& X) const;
	Vector<Real, STATE_NF> y_vector_flux(const _3Vec& X) const;
	Vector<Real, STATE_NF> z_vector_flux(const _3Vec& X) const;
	Real scalar_flux(const _3Vec& X) const {
		return pg();
	}
	Real vx() const {
		return sx() / rho();
	}
	Real vy() const {
		return sy() / rho();
	}
	Real vz() const {
		return sz() / rho();
	}
	static Vector<Real, STATE_NF> x_scalar_flux_coeff(const _3Vec& X) {
		Vector<Real, STATE_NF> c = 0.0;
		c[sx_index] = 1.0;
		return c;
	}
	static Vector<Real, STATE_NF> y_scalar_flux_coeff(const _3Vec& X) {
		Vector<Real, STATE_NF> c = 0.0;
		c[sy_index] = 1.0;
		return c;
	}
	static Vector<Real, STATE_NF> z_scalar_flux_coeff(const _3Vec& X) {
		Vector<Real, STATE_NF> c = 0.0;
		c[sz_index] = 1.0;
		return c;
	}
	void to_prim(const _3Vec& X);
	void from_prim(const _3Vec& X);
	void to_con(const _3Vec& X);
	void from_con(const _3Vec& X);
	static const char* field_name(int i);
	void floor(const _3Vec&);
	Real rho() const;
	Real tau() const;
	Real sz() const;
	Real sx() const;
	Real sy() const;
	Real pg() const;
	Real et() const;
	Real ei() const;
	void set_rho(Real);
	void set_sx(Real);
	void set_sy(Real);
	void set_sz(Real);
	void set_et(Real);
	void set_tau(Real);
	Real speed() {
		return sqrt(sx() * sx() + sy() * sy() + sz() * sz()) / rho();
	}
};

typedef EulerState State;

#endif

#endif /* EULER_STATE_H_ */
