/*
 * poisson_test.cpp
 *
 *  Created on: Mar 6, 2013
 *      Author: dmarce1
 */

#include "radiation_test.h"

#ifdef RADIATION_TEST

RadiationTest* RadiationTest::new_octnode() const {
	return new RadiationTest;
}

void RadiationTest::set_refine_flags() {
	ChildIndex c;
	if (get_level() < 1) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			set_refine_flag(i, true);
		}
	} else if (get_level() < OctNode::get_max_level_allowed()) {
		for (int k = 1; k < PNX - 1; k++) {
			c.set_z(2 * k / PNX);
			for (int j = 1; j < PNX - 1; j++) {
				c.set_y(2 * j / PNX);
				for (int i = 1; i < PNX - 1; i++) {
					c.set_x(2 * i / PNX);
					if (radius(i, j, k) < 0.25)
						set_refine_flag(c, true);
				}
			}
		}
	}
}

void RadiationTest::initialize() {
	for (int k = 0; k < PNX; k++) {
		for (int j = 0; j < PNX; j++) {
			for (int i = 0; i < PNX; i++) {
	//			set_erad(i, j, k, 1.0);
//				if (radius(i, j, k) < 0.25) {
//					set_rho(i, j, k, 1.0e+2);
//				} else {
//					set_rho(i, j, k, 1.0e-3);

			//	}
	//			if (radius(i, j, k) < 0.015) {
	//				set_eint(i, j, k, 1.0);
	//			} else {
	//				set_eint(i, j, k, 1.0e-7);
	//			}
				set_cv(i, j, k, 1.0);
				set_erad(i, j, k, Real((rand() & 0xFF) + 1) / Real(0xFF + 1));
				set_eint(i, j, k, Real((rand() & 0xFF) + 1) / Real(0xFF + 1));
				set_rho(i, j, k, Real((rand() & 0xFF) + 1) / Real(0xFF + 1));
			}
		}
	}
}

Real RadiationTest::get_output_point(int i, int j, int k, int l) const {
	switch (l) {
	case 0:
		return get_erad(i, j, k);
	case 1:
		return get_eint(i, j, k);
	case 2:
		return get_Trad(i, j, k);
	case 3:
		return get_Tgas(i, j, k);
	case 4:
		return get_chi(i, j, k);
	case 5:
		return get_kappa_E(i, j, k);
	}
	return 1.0 / 0.0;
}

int RadiationTest::nvar_output() const {
	return 6;
}

const char* RadiationTest::output_field_names(int i) const {
	switch (i) {
	case 0:
		return "erad";
	case 1:
		return "eint";
	case 2:
		return "Trad";
	case 3:
		return "Tgas";
	case 4:
		return "chi";
	case 5:
		return "kappa_E";
	}
	return "dummy";
}

RadiationTest::RadiationTest() {
// TODO Auto-generated constructor stub

}

RadiationTest::~RadiationTest() {
// TODO Auto-generated destructor stub
}

#endif
