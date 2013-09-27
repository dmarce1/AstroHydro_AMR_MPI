#ifndef POISSON_TEST_H_
#define POISSON_TEST_H_

#include "../poisson/poisson.h"
#include "../hydro_FMM_grid/hydro_FMM_grid.h"
#ifdef POISSON_TEST
#include "poisson_test_defs.h"
#ifdef USE_FMM
class PoissonTest: public HydroFMMGrid {
#else
	class PoissonTest: public Poisson {
#endif
public:
	static Real mass0;
	static Real x0, y0, z0;
	static Real r0;
	static void compute_error(Real* phi_err, _3Vec* force_error, Real* mx_err, Real* my_err, Real* mz_err);
	virtual PoissonTest* new_octnode() const;
	virtual Real get_output_point(int i, int j, int k, int l) const;
	virtual int nvar_output() const;
	virtual const char* output_field_names(int i) const;
	Real analytic_phi(Real x, Real y, Real z) const;
	_3Vec analytic_force(Real x, Real y, Real z) const;
	void initialize();
	Real phi_error();
	_3Vec force_error();
	Vector<Real, 3> momentum_error();
public:

	virtual void set_refine_flags() {
		ChildIndex c;
		if (get_level() < 1) {
			for (int i = 0; i < OCT_NCHILD; i++) {
				set_refine_flag(i, true);
			}
		} else if (get_level() < OctNode::get_max_level_allowed()) {
#ifdef USE_FMM
			for (int k = BW; k < GNX - BW; k++) {
				c.set_z(2 * k / GNX);
				for (int j = BW; j < GNX - BW; j++) {
					c.set_y(2 * j / GNX);
					for (int i = BW; i < GNX - BW; i++) {
						c.set_x(2 * i / GNX);
						if (!get_refine_flag(c)) {
							if ((*this)(i, j, k).rho() > 1.0e-3) {
#else
								for (int k = 1; k < PNX - 1; k++) {
									c.set_z(2 * k / PNX);
									for (int j = 1; j < PNX - 1; j++) {
										c.set_y(2 * j / PNX);
										for (int i = 1; i < PNX - 1; i++) {
											c.set_x(2 * i / PNX);
											if (!get_refine_flag(c)) {
												if (S(i, j, k) > 1.0e-3) {
#endif
								set_refine_flag(c, true);
							}
						}
					}
				}
			}
		}
	}

	static void run(int argc, char* argv[]);
	virtual ~PoissonTest();
	PoissonTest();
};

#endif /* POISSON_TEST_H_ */
#endif
