#ifndef POISSON_TEST_H_
#define POISSON_TEST_H_

#include "../poisson/poisson.h"
#include "poisson_test_defs.h"

#ifdef POISSON_TEST
class PoissonTest: public Poisson {
public:
	static Real mass0;
	static Real x0, y0, z0;
	static Real r0;
	static void compute_error(Real* phi_err, Real* mx_err, Real* my_err, Real* mz_err);
	virtual PoissonTest* new_octnode() const;
	virtual Real get_output_point(int i, int j, int k, int l) const;
	virtual int nvar_output() const;
	virtual const char* output_field_names(int i) const;
	Real analytic_phi(Real x, Real y, Real z) const;
	void initialize();
	Real phi_error();
	Vector<Real, 3> momentum_error();
public:
	static void run(int argc, char* argv[]);
	virtual ~PoissonTest();
	PoissonTest();
};

#endif /* POISSON_TEST_H_ */
#endif
