#include "radiation_test.h"
#include "../physical_constants.h"

#ifdef RADIATION_TEST

void RadiationTest::run(int argc, char* argv[]) {
	PhysicalConstants::c = 1.0;
	PhysicalConstants::sigma = 1.0;
	Real dtrad;
	for (int iters = 0; iters <= OctNode::get_max_level_allowed(); iters++) {
		OctNode::initialize_grids();
		check_for_refine();
		if (MPI_rank() == 0) {
			printf("maxlevel = %i nodes = %i\n", get_max_level(), get_node_cnt());
		}
	}
	OctNode::initialize_grids();
	get_root()->output("X", 0, PNX, 1);
	int iter = 1;
	const Real dxmin = dynamic_cast<RadiationTest*>(get_root())->get_dx() / Real(1 << OctNode::get_max_level_allowed());
	dynamic_cast<RadiationTest*>(get_root())->phi_calc_amr_bounds();
	const Real dt0 = dxmin / PhysicalConstants::c;
	Real dt = dt0 * 1.0e-4;
	for (Real t = 0.0; t < 1.0e+3; t += dt) {
		if (MPI_rank() == 0) {
			printf("t=%e dt=%e dtrad=%e\n", t, dt, dtrad);
		}
		dtrad = implicit_solve(dt);
		if (dt < dt0) {
			dt *= pow(10.0, 0.1);
		}
	//	dt = min(dtrad, dt);
		get_root()->output("X", iter, PNX, 1);
		iter++;
	}

}
#endif
