/*
 * hydro_rad_grid.cpp
 *
 *  Created on: May 22, 2013
 *      Author: dmarce1
 */

#include "hydro_rad_grid.h"
#ifdef USE_RADITION

HydroRadGrid::HydroRadGrid() {
	// TODO Auto-generated constructor stub

}

HydroRadGrid::~HydroRadGrid() {
	// TODO Auto-generated destructor stub
}




void HydroRadGrid::run(int argc, char* argv[]) {
	Real start_time, dt, tm;
	FILE* fp;
	bool do_output, last_step;
	int step_cnt = 0;
	int ostep_cnt = 0;
	int nnodes;
	Real avg_node_cnt;
	int max_node_cnt, min_node_cnt;

	setup_grid_structure();
	for (int i = 0; i < OctNode::get_local_node_cnt(); i++) {
		dynamic_cast<MultiGrid*>(OctNode::get_local_node(i))->phi_calc_amr_bounds();
	}

	if (OUTPUT_TIME_FREQ <= TIME_MAX) {
		get_root()->output("X", 0.0, GNX, BW);
	}

	avg_node_cnt = 0.0;
	max_node_cnt = min_node_cnt = OctNode::get_node_cnt();

	start_time = MPI_Wtime();
	do {
		dt = next_dt(&do_output, &last_step, &ostep_cnt);
		step(dt);
		if (step_cnt % (GNX - 2 * BW)== 0){
			//	check_for_refine();
			if (MultiGrid::check_for_refine()) {
				HydroGrid::redistribute_grids();
			}
		}
		step_cnt++;
		if (MPI_rank() == 0) {
			nnodes = OctNode::get_node_cnt();
			avg_node_cnt = (avg_node_cnt * Real(step_cnt - 1) + Real(nnodes)) / Real(step_cnt);
			max_node_cnt = max(max_node_cnt, nnodes);
			min_node_cnt = min(min_node_cnt, nnodes);
			printf("step=%i t=%e dt=%e lmax=%i ngrids=%i avg ngrids=%i", step_cnt, HydroGrid::get_time(), dt, OctNode::get_max_level(), nnodes,
					(int) avg_node_cnt);
		}
		if (do_output) {
			if (MPI_rank() == 0) {
				printf("*");
			}
			get_root()->output("X", nint(HydroGrid::get_time() / OUTPUT_TIME_FREQ), GNX, BW);
		}
		if (MPI_rank() == 0) {
			printf("\n");
		}
	} while (!last_step);

	if (MPI_rank() == 0) {
		char str[81];
		sprintf(str, "scaling.%i.txt", OctNode::get_max_level_allowed());
		tm = (MPI_Wtime() - start_time);
		fp = fopen(str, "at");
		fprintf(fp, "%i %e %i %i %i\n", MPI_size(), tm, (int) avg_node_cnt, min_node_cnt, max_node_cnt);
		fclose(fp);
	}
}
#endif
