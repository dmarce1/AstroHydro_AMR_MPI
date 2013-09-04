#include "single_star.h"

#ifdef HYDRO_GRAV_GRID

Vector<Real, 2>* SingleStar::radial_avg_tmp;
Vector<Real, 2>* SingleStar::radial_avg;
int* SingleStar::radial_bin_cnt;
int SingleStar::radial_N;
int* SingleStar::radial_bin_cnt_tmp;

#ifdef SINGLE_STAR
void SingleStar::read_from_file(const char* str, int* i1, int* i2) {
	if (MPI_rank() == 0) {
		printf("Reading checkpoint file\n");
	}
	PhysicalConstants::set_cgs();
	char* fname;
	Real omega;
	FILE* fp;
	_3Vec O;
	asprintf(&fname, "checkpoint.%s.%i.bin", str, MPI_rank());
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		printf("Error - Checkpoint file not found!");
		abort();
	}
	free(fname);
//	fread(&refine_floor, sizeof(Real), 1, fp);
	fread(&State::ei_floor, sizeof(Real), 1, fp);
	fread(&State::rho_floor, sizeof(Real), 1, fp);
//	fread(&com_vel_correction, sizeof(_3Vec), 1, fp);
	fread(&O, sizeof(_3Vec), 1, fp);
	fread(&omega, sizeof(Real), 1, fp);
//	set_origin(O);
	fread(&last_dt, sizeof(Real), 1, fp);
	fread(i1, sizeof(int), 1, fp);
	fread(i2, sizeof(int), 1, fp);
	get_root()->read_checkpoint(fp);
	fclose(fp);
	State::set_omega(omega);
	int count = OctNode::get_local_node_cnt();
	for (int i = 0; i < count; i++) {
		dynamic_cast<MultiGrid*>(OctNode::get_local_node(i))->phi_calc_amr_bounds();
	}
	solve_poisson();
}

void SingleStar::write_to_file(int i1, int i2, const char* idname) {
	char* fname;
	Real omega = State::get_omega();
	FILE* fp;
	_3Vec O = get_origin();

//	asprintf(&fname, "mv -f checkpoint.hello.%i.bin checkpoint.goodbye.%i.bin 2> /dev/null\n", MPI_rank(), MPI_rank());
//	system(fname);
//	free(fname);

	asprintf(&fname, "checkpoint.%s.%i.bin", idname, MPI_rank());
	fp = fopen(fname, "wb");
	free(fname);
//	fwrite(&refine_floor, sizeof(Real), 1, fp);
	fwrite(&State::ei_floor, sizeof(Real), 1, fp);
	fwrite(&State::rho_floor, sizeof(Real), 1, fp);
//	fwrite(&com_vel_correction, sizeof(_3Vec), 1, fp);
	fwrite(&O, sizeof(_3Vec), 1, fp);
	fwrite(&omega, sizeof(Real), 1, fp);
	fwrite(&last_dt, sizeof(Real), 1, fp);
	fwrite(&i1, sizeof(int), 1, fp);
	fwrite(&i2, sizeof(int), 1, fp);
	get_root()->write_checkpoint(fp);
	fclose(fp);
}
void SingleStar::run(int argc, char* argv[]) {
	//State::turn_off_total_energy();
	Real dxmin = dynamic_cast<SingleStar*>(get_root())->get_dx() / Real(1 << get_max_level_allowed());
	radial_N = int(2.5e+9 / dxmin + 1.5);
	PhysicalConstants::set_cgs();
	shadow_off();

	Real dt;
	bool do_output, last_step;
	int step_cnt = 0;
	int ostep_cnt = 0;
	int nnodes;


	if (argc == 2) {
		setup_grid_structure();
		if (OUTPUT_TIME_FREQ <= TIME_MAX ) {
			get_root()->output("X", 0.0, GNX, BW);
		}
	} else {
		read_from_file(argv[2], &step_cnt, &ostep_cnt);
	}


	hydro_time = poisson_boundary_time = poisson_interior_time = 0.0;

	do {
		if (step_cnt % CHECKPT_FREQ == 0) {
			if (step_cnt % (CHECKPT_FREQ * 2) == 0) {
				write_to_file(step_cnt, ostep_cnt, "hello");
			} else {
				write_to_file(step_cnt, ostep_cnt, "goodbye");
			}
		}
		dt = next_dt(&do_output, &last_step, &ostep_cnt);
		step(dt);
		if (step_cnt % (GNX - 2 * BW)== 0){
			check_for_refine();
		}
		step_cnt++;
		if (MPI_rank() == 0) {
			nnodes = OctNode::get_node_cnt();
			printf("step=%i t=%e dt=%e lmax=%i ngrids=%i", step_cnt, HydroGrid::get_time(), dt, OctNode::get_max_level(), nnodes);
		}
		if (do_output) {
			if (MPI_rank() == 0) {
				printf("*\n");
			}
			get_root()->output("X", nint(HydroGrid::get_time() / OUTPUT_TIME_FREQ), GNX, BW);
		}
		if (MPI_rank() == 0) {
			printf("\n");
		}
		Vector<State, 4> vec = dynamic_cast<HydroGrid*>(get_root())->state_sum();
		State sum = vec[0];
		_3Vec com = get_center_of_mass();
		if (MPI_rank() == 0) {
			FILE* fp = fopen("sums.dat", "at");
			fprintf(fp, "%e %e %e %e %e %e %e %e\n", get_time(), sum[0], sum[1], sum[2], sum[3], sum[4] + 0.5 * sum[6], sum[State::frac_index + 0],
					sum[State::frac_index + 1]);
			fclose(fp);
			fp = fopen("com.dat", "at");
			fprintf(fp, "%e %e %e %e\n", get_time(), com[0], com[1], com[2]);
			fclose(fp);
		}
	} while (!last_step);
	if (MPI_rank() == 0) {
		char* str;
		if (asprintf(&str, "time.%i.txt", get_max_level_allowed()) == 0) {
			printf("Unable to create filename\n");
		} else {
			FILE* fp = fopen(str, "at");
			fprintf(fp, "%i %e %e %e %e\n", MPI_size(), hydro_time + poisson_boundary_time + poisson_interior_time, hydro_time, poisson_boundary_time,
					poisson_interior_time);
			fclose(fp);
			free(str);
		}
	}
}
#endif
#endif
