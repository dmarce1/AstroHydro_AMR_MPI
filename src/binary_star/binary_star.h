/*
 * binary_star.h
 *
 *  Created on: Mar 20, 2013
 *      Author: dmarce1
 */

#ifndef BINARY_STAR_H_
#define BINARY_STAR_H_

#include "../hydro_grav_grid/hydro_grav_grid.h"
#include "../hydro_FMM_grid/hydro_FMM_grid.h"
#include "dwd.h"

#ifdef HYDRO_GRAV_GRID
#ifdef USE_FMM
class BinaryStar: public HydroFMMGrid {
public:
	static binary_parameters_t bparam;
	static bool bparam_init;
#else
	class BinaryStar: public HydroGravGrid {
#endif
private:
	static int preferred_node_count;
	static Real refine_floor;
	static Real code_to_cm, code_to_s, code_to_K, code_to_g;
	static Real lz_t0;
	static _3Vec com_vel_correction;
	static Real M1, M2, Ka, Kd;
	static bool scf_code;
	static Real q_ratio, polyK, rho_max_a, rho_max_d;
	static void read_from_file(const char*, int*, int*);
	static void write_to_file(int, int, const char*);
	static _3Vec a0, d0;
	virtual void set_refine_flags();
	Real radius(int i, int j, int k);
	static void assign_fracs(Real donor_frac, Real min_phi, Real max_phi);
	static void find_l(Real m1_x, Real m2_x, Real* l1_phi, Real* l1_x, int lnum);
	static void next_omega(Real phi0);
	static Real virial_error();
	static Real find_acc_com();
	static void find_phimins(Real* phi_min_a, Real* a_x, Real* phi_min_d, Real* d_x);
	static Real Ax, Bx, Cx, Aphi, Bphi, Cphi;
	static void find_rho_max(Real* rho1, Real* rho2);
	static void find_phi_min(Real* phi1, Real* phi2);
	static Real find_K(int frac, Real phi0, Real center, Real l1_x, Real* span);
	static void next_rho(Real, Real, Real, Real, Real, Real, Real);
	static void find_mass(int, Real*, Real*);
	virtual void compute_flow_off();
public:
	Real grad(int i, int j, int k) {
		_3Vec x = X(i, j, k);
		Real R = sqrt(x[0] * x[0] + x[1] * x[1]);
		return (gx(i, j, k) * x[0] + gy(i, j, k) * x[1]) / R;
	}
	_3Vec geff(int i, int j, int k) const {
		_3Vec g;
		g[0] = gx(i, j, k);
		g[1] = gy(i, j, k);
		g[2] = gz(i, j, k);
		Real O2 = pow(State::get_omega(), 2);
		g[0] += HydroGrid::xc(i) * O2;
		g[1] += HydroGrid::yc(j) * O2;
		return g;
	}
	Real glz(int i, int j, int k) const {
		return HydroGrid::xc(i) * gy(i, j, k) - HydroGrid::yc(j) * gx(i, j, k);
	}
	static void integrate_conserved_variables(Vector<Real, STATE_NF>*);
	static void diagnostics(Real);
	virtual BinaryStar* new_octnode() const {
		return new BinaryStar;
	}
	virtual void initialize();
	static void run(int argc, char* argv[]);
	static void scf_run(int argc, char* argv[]);
	BinaryStar();
	virtual ~BinaryStar();
	virtual Real get_output_point(int i, int j, int k, int l) const {
		State U = (*this)(i, j, k);
		_3Vec x = X(i, j, k);
		switch (l) {
		case 0:
			return U.rho();
		case 1:
			return U.frac(1);
		case 2:
			return U.frac(0);
		case 3:
			return U.vx(x);
		case 4:
			return U.vy(x);
		case 5:
			return U.vz();
		case 6:
			return U.temp(x);
#ifndef USE_FMM
			case 7:
			return get_phi(i - BW + 1, j - BW + 1, k - BW + 1);
			case 8:
			return get_fx(i - BW + 1, j - BW + 1, k - BW + 1);
			case 9:
			return get_fy(i - BW + 1, j - BW + 1, k - BW + 1);
			case 10:
			return get_fz(i - BW + 1, j - BW + 1, k - BW + 1);
		}
#else
		case 7:
			return get_phi(i, j, k);
		case 8:
			return gx(i, j, k);
		case 9:
			return gy(i, j, k);
		case 10:
			return gz(i, j, k);
		case 11:
			return glz(i, j, k);
		case 12:
			return (*this)(i, j, k).pg(this->X(i, j, k));
		case 13:
			return (*this)(i, j, k).ei(this->X(i, j, k)) + (*this)(i, j, k).ed();
		}
#endif
		assert(false);
		return 0.0;
	}
	virtual const char* output_field_names(int i) const {
		switch (i) {
		case 0:
			return "rho";
		case 1:
			return "He";
		case 2:
			return "CO";
		case 3:
			return "vx";
		case 4:
			return "vy";
		case 5:
			return "vz";
		case 6:
			return "T";
		case 7:
			return "phi";
		case 8:
			return "gx";
		case 9:
			return "gy";
		case 10:
			return "gz";
		case 11:
			return "glz";
		case 12:
			return "p";
		case 13:
			return "e_gas";
		} assert(false);
		return "";
	}
	virtual int nvar_output() const {
		return 14;
	}
};
#endif

#endif /* BINARY_STAR_H_ */
