#ifndef RADPOISSON_TEST_H_
#define RADPOISSON_TEST_H_

#include "../radiation/radiation.h"
#include "../hydro_rad_grid/hydro_rad_grid.h"
#include "defs.h"

#ifdef RADIATION_TEST

class RadiationTest: public Radiation {
private:
	Real radius(int i, int j, int k) const {
		return sqrt(xc(i) * xc(i) + yc(j) * yc(j) + zc(k) * zc(k));
	}
	virtual RadiationTest* new_octnode() const;
	virtual Real get_output_point(int i, int j, int k, int l) const;
		virtual int nvar_output() const;
		virtual const char* output_field_names(int i) const;
		void initialize();
public:
	static void run(int argc, char* argv[]);
	virtual ~RadiationTest();
	RadiationTest();
	virtual void set_refine_flags();
};

#endif /* POISSON_TEST_H_ */

#endif
