#include "poisson_test.h"

#ifdef POISSON_TEST

void PoissonTest::run(int argc, char* argv[]) {
	Real phi_err, x_err, y_err, z_err;
	int iters = 0;
	iters = 0;
	Real boundary_time, solve_time, start_time;
	get_root()->find_local_nodes();
	for (int iters = 0; iters <= OctNode::get_max_level_allowed(); iters++) {
		OctNode::initialize_grids();
		check_for_refine();
		if (MPI_rank() == 0) {
			printf("maxlevel = %i nodes = %i\n", get_max_level(), get_node_cnt());
		}
	}
	OctNode::initialize_grids();

	start_time = MPI_Wtime();
	compute_physical_boundaries();
	PoissonTest::mass0 = 4.0 * M_PI * source_sum();
	printf("mass=%e\n",mass0);
	boundary_time = MPI_Wtime() - start_time;

	start_time = MPI_Wtime();
	Real err = vcycle();
	for (iters = 0; err > 1.0e-6; iters++) {
		err = vcycle();
		if (MPI_rank() == 0) {
			printf("%i %e\n", iters, err);
		}
	}
	solve_time = MPI_Wtime() - start_time;

	if (MPI_rank() == 0) {
		printf("Boundary time = %e : Solve time = %e\n", boundary_time, solve_time);
	}

	get_root()->output("X", 0, PNX, 1);
	compute_error(&phi_err, &x_err, &y_err, &z_err);
	if (MPI_rank() == 0) {
		printf("%e %e %e %e %e %e %e\n", PoissonTest::x0, PoissonTest::y0, PoissonTest::z0, phi_err, x_err, y_err, z_err);
	}
}

void PoissonTest::compute_error(Real* phi_err, Real* mx_err, Real* my_err, Real* mz_err) {
	int count;
	Real phi_sum = 0.0;
	Vector<Real, 3> m_sum = 0.0;
	Real v[4], v0[4];
	count = OctNode::get_local_node_cnt();
	for (int i = 0; i < count; i++) {
		phi_sum += dynamic_cast<PoissonTest*>(OctNode::get_local_node(i))->phi_error();
		m_sum += dynamic_cast<PoissonTest*>(OctNode::get_local_node(i))->momentum_error();
	}
	v[0] = m_sum[0];
	v[1] = m_sum[1];
	v[2] = m_sum[2];
	v[3] = phi_sum;
	MPI_Reduce(&v, &v0, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
	*mx_err = v0[0];
	*my_err = v0[1];
	*mz_err = v0[2];
	*phi_err = sqrt(v0[3] / 8.0);
}

#endif
