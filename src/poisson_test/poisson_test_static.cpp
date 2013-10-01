#include "poisson_test.h"

#ifdef POISSON_TEST

void PoissonTest::run(int argc, char* argv[]) {
	Real phi_err, x_err, y_err, z_err;
	int iters = 0;
	iters = 0;
	Real boundary_time, solve_time, start_time;
	get_root()->find_local_nodes();
	HydroFMMGrid::solve_on = false;
	for (int iters = 0; iters <= OctNode::get_max_level_allowed(); iters++) {
		OctNode::initialize_grids();
		check_for_refine();
		if (MPI_rank() == 0) {
			printf("maxlevel = %i nodes = %i\n", get_max_level(), get_node_cnt());
		}
	}
	OctNode::initialize_grids();

	start_time = MPI_Wtime();
#ifdef USE_FMM
	PoissonTest::mass0 = dynamic_cast<HydroFMMGrid*>(get_root())->com_sum()[3];
#else
	compute_physical_boundaries();
	PoissonTest::mass0 = 4.0 * M_PI * source_sum();
#endif
	printf("mass=%e\n", mass0);
	boundary_time = MPI_Wtime() - start_time;

	start_time = MPI_Wtime();
#ifdef USE_FMM
	HydroFMMGrid::solve_on = true;
	FMM_solve();
	HydroFMMGrid::solve_on = false;
#else
	Real err = vcycle();
	for (iters = 0; err > 1.0e-6; iters++) {
		err = vcycle();
		if (MPI_rank() == 0) {
			printf("%i %e\n", iters, err);
		}
	}
#endif
	solve_time = MPI_Wtime() - start_time;

	if (MPI_rank() == 0) {
		printf("Boundary time = %e : Solve time = %e\n", boundary_time, solve_time);
	}
#ifdef USE_FMM
	get_root()->output("X", 0, GNX, BW);
#else
	get_root()->output("X", 0, PNX, 1);
#endif
	_3Vec ferr;
	compute_error(&phi_err, &ferr, &x_err, &y_err, &z_err);
	if (MPI_rank() == 0) {
		printf("%e %e %e %e %e %e %e\n", PoissonTest::x0, PoissonTest::y0, PoissonTest::z0, phi_err, x_err, y_err, z_err);
		printf("%e %e %e \n", ferr[0], ferr[1], ferr[2]);
	}
}

void PoissonTest::compute_error(Real* phi_err, _3Vec* fee, Real* mx_err, Real* my_err, Real* mz_err) {
	int count;
	Real phi_sum = 0.0;
	Vector<Real, 3> m_sum = 0.0;
	Real v[4], v0[4];
	Real f[3];
	Real f0[3];
	_3Vec f_sum = 0.0;
	count = OctNode::get_local_node_cnt();
	for (int i = 0; i < count; i++) {
		phi_sum += dynamic_cast<PoissonTest*>(OctNode::get_local_node(i))->phi_error();
		m_sum += dynamic_cast<PoissonTest*>(OctNode::get_local_node(i))->momentum_error();
		f_sum += dynamic_cast<PoissonTest*>(OctNode::get_local_node(i))->force_error();
	}
	v[0] = m_sum[0];
	v[1] = m_sum[1];
	v[2] = m_sum[2];
	v[3] = phi_sum;
	MPI_Reduce(v, v0, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD );
	*mx_err = v0[0];
	*my_err = v0[1];
	*mz_err = v0[2];
	*phi_err = sqrt(v0[3] / 8.0);
	f0[0] = f_sum[0];
	f0[1] = f_sum[1];
	f0[2] = f_sum[2];
	MPI_Reduce(f0, f, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD );
	(*fee)[0] = f[0];
	(*fee)[1] = f[1];
	(*fee)[2] = f[2];

}

#endif
