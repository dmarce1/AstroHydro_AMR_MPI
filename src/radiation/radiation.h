/*
 * radiation.h
 *
 *  Created on: May 17, 2013
 *      Author: dmarce1
 */

#ifndef RADIATION_H_
#define RADIATION_H_

#include "../multigrid/multigrid.h"
#include "../physical_constants.h"

#ifdef USE_RADIATION

#define RAD_INS_CNT 7

class Radiation: public MultiGrid {
private:
	int ip;
	void null();
	typedef void (Radiation::*ifunc_t)();
	static void boundary_communicate();
	static ifunc_t cs[RAD_INS_CNT + 1];
	_3Vec* mpi_3vec_buffer[24];
	static MPI_Datatype MPI_send_bnd_t[6];
	static MPI_Datatype MPI_recv_bnd_t[6];
	int cx, ax;
	static void MPI_datatypes_init();
	void gradE_real_boundary_communicate();
	void gradE_real_boundary_begin_loop();
	void gradE_real_boundary_end_loop();
	void gradE_real_boundary_wait();
	void gradE_amr_boundary_communicate();
	void gradE_amr_boundary_wait();
	void compute_gradE();
	static Real change_max();
	static Real dt;
	static bool initialized;
	static Real const sigma_T = 0.2;
	Array3d<_3Vec, PNX, PNX, PNX> gradE;
	Array3d<Real, PNX, PNX, PNX> erad0;
	Array3d<Real, PNX, PNX, PNX> erad1;
	Array3d<Real, PNX, PNX, PNX> eint1;
	Array3d<Real, PNX, PNX, PNX> rho;
	Array3d<Real, PNX, PNX, PNX> D1;
	Array3d<Real, PNX, PNX, PNX> D2;
	Array3d<Real, PNX, PNX, PNX> D3;
	Array3d<Real, PNX, PNX, PNX> eint0;
	Array3d<Real, PNX, PNX, PNX> eint;
	Array3d<Real, PNX, PNX, PNX> Bp;
	Array3d<Real, PNX, PNX, PNX> cv;
	static Real flux_limiter(Real);
	void compute_eint();
	void compute_RHS();
	virtual void residual_error_compute(int);
	virtual void relax_compute(int);
	virtual void vdown_init_compute(int);
	virtual void vdown_init_adjust_send(int);
	virtual void vdown_init_adjust_recv_wait(int);
	void jacobian_compute(Real);
	virtual void forces_compute(int);
	Real erad(int i, int j, int k) const;
	Real& erad(int i, int j, int k);
	Real erad(const _3Vec& x) const {
		return erad(x[0], x[1], x[2]);
	}
	void compute_coupling();
	Real& erad(const _3Vec& x) {
		return erad(x[0], x[1], x[2]);
	}
	Real* Jacobian(int, int, int) const;
	void store0();
	void store1();
protected:
	_3Vec get_rad_momentum_term(int i, int j, int k) const {
		Real R;
		_3Vec g;
		g[0] = (erad(i + 1, j, k) - erad(i - 1, j, k)) / get_dx() / 2.0;
		g[1] = (erad(i, j + 1, k) - erad(i, j - 1, k)) / get_dx() / 2.0;
		g[2] = (erad(i, j, k + 1) - erad(i, j, k - 1)) / get_dx() / 2.0;
		R = sqrt(g.dot(g)) / erad(i, j, k) / chi(rho(i, j, k));
		return -g * flux_limiter(R);
	}
	static Real chi(Real rho);
	static Real kappa_E(Real rho);
	static Real kappa_p(Real rho);
	void compute_difco();
public:
	virtual void allocate_arrays();
	virtual void deallocate_arrays();
	Real get_eint(int i, int j, int k) const {
		return eint(i, j, k);
	}
	Real get_lambda(int i, int j, int k) const {
		Real gx, gy, gz, R;
		gx = (erad(i + 1, j, k) - erad(i - 1, j, k)) / get_dx() / 2.0;
		gy = (erad(i, j + 1, k) - erad(i, j - 1, k)) / get_dx() / 2.0;
		gz = (erad(i, j, k + 1) - erad(i, j, k - 1)) / get_dx() / 2.0;
		R = sqrt(gx * gx + gy * gy + gz * gz) / erad(i, j, k) / chi(rho(i, j, k));
		return flux_limiter(R);
	}
	Real get_Trad(int i, int j, int k) const {
		static const Real a = 4.0 * PhysicalConstants::sigma / PhysicalConstants::c;
		return pow(erad(i, j, k) / a, 0.25);
	}
	Real get_Tgas(int i, int j, int k) const {
		return eint(i, j, k) / cv(i, j, k);
	}
	void set_erad(int i, int j, int k, Real e) {
		erad(i, j, k) = e;
	}
	void set_eint(int i, int j, int k, Real e) {
		eint(i, j, k) = e;
	}
	void set_cv(int i, int j, int k, Real e) {
		cv(i, j, k) = e;
	}
	void set_rho(int i, int j, int k, Real e) {
		rho(i, j, k) = e;
	}
	Real get_erad(int i, int j, int k) const {
		return erad(i, j, k);
	}
	Real get_chi(int i, int j, int k) const {
		return chi(rho(i, j, k));
	}
	Real get_kappa_E(int i, int j, int k) const {
		return kappa_E(rho(i, j, k));
	}
	Vector<Vector<Real, 3>, 3> get_eddington_tensor(int i, int j, int k) const {
		Vector<Vector<Real, 3>, 3> F;
		_3Vec g;
		Real R, f, lambda, f1, f2, abs_g;
		g[0] = (erad(i + 1, j, k) - erad(i - 1, j, k)) / get_dx() / 2.0;
		g[1] = (erad(i, j + 1, k) - erad(i, j - 1, k)) / get_dx() / 2.0;
		g[2] = (erad(i, j, k + 1) - erad(i, j, k - 1)) / get_dx() / 2.0;
		abs_g = sqrt(g.dot(g));
//		printf("%e %e\n", erad(i, j, k), rho(i, j, k));
		R = abs_g / erad(i, j, k) / chi(rho(i, j, k));
		lambda = flux_limiter(R);
		f = lambda * (1.0 + lambda * R * R);
		f1 = (1.0 - f) / 2.0;
		f2 = (3.0 * f - 1.0) / 2.0;
		if (abs_g != 0.0) {
			g /= abs_g;
		} else {
			g = 0.0;
		}
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				F[i][j] = f2 * g[i] * g[j];
				if (i == j) {
					F[i][j] += f1;
				}
			}
		}
		return F;
	}
	static Real implicit_solve(Real dt);
	Radiation();
	virtual ~Radiation();
};
#endif
#endif /* RADIATION_H_ */
