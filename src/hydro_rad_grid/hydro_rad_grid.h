/*
 * hydro_rad_grid.h
 *
 *  Created on: May 22, 2013
 *      Author: dmarce1
 */

#ifndef HYDRO_RAD_GRID_H_
#define HYDRO_RAD_GRID_H_

#include "../radiation/radiation.h"
#include "../hydro_grid/hydro_grid.h"
#include "../indexer3d.h"
#ifdef USE_RADIATION

class HydroRadGrid: public HydroGrid, public Radiation {
private:
	virtual void max_dt_compute(int dir) {
		Real this_dt;
		Real dtinv = 0.0;
		_3Vec x;
		State q0[GNX], ql[GNX], qr[GNX];
		int k, j, i;
		Real this_dtinv, dx, f;
		dx = HydroGrid::get_dx();
		for (k = BW; k < GNX - BW; k++) {
			for (j = BW; j < GNX - BW; j++) {
				this_dtinv = 0.0;
				if (dir == 0) {
					for (i = 0; i < GNX; i++) {
						q0[i] = (*this)(i, j, k);
						q0[i].to_prim(HydroGrid::X(i, j, k));
					}
					reconstruct(q0, ql, qr);
					for (i = BW; i < GNX - BW + 1; i++) {
						x = HydroGrid::Xfx(i, j, k);
						ql[i].from_prim(x);
						qr[i].from_prim(x);
						f = 1.0 - exp(-dx * chi(0.5 * (ql[i].rho() + qr[i].rho())));
						if( f < 0.0 ) {
							printf( "%e\n", f);
						}
						this_dtinv = max(this_dtinv, ql[i].max_abs_x_eigen(x, f), qr[i].max_abs_x_eigen(x, f));
					}
				} else if (dir == 1) {
					for (i = 0; i < GNX; i++) {
						q0[i] = (*this)(j, i, k);
						q0[i].to_prim(HydroGrid::X(j, i, k));
					}
					reconstruct(q0, ql, qr);
					for (i = BW; i < GNX - BW + 1; i++) {
						x = HydroGrid::Xfy(j, i, k);
						ql[i].from_prim(x);
						qr[i].from_prim(x);
						f = 1.0 - exp(-dx * chi(0.5 * (ql[i].rho() + qr[i].rho())));
						this_dtinv = max(this_dtinv, ql[i].max_abs_y_eigen(x, f), qr[i].max_abs_y_eigen(x, f));
					}
				} else {
					for (i = 0; i < GNX; i++) {
						q0[i] = (*this)(j, k, i);
						q0[i].to_prim(HydroGrid::X(j, k, i));
					}
					reconstruct(q0, ql, qr);
					for (i = BW; i < GNX - BW + 1; i++) {
						x = HydroGrid::Xfz(j, k, i);
						ql[i].from_prim(x);
						qr[i].from_prim(x);
						f = 1.0 - exp(-dx * chi(0.5 * (ql[i].rho() + qr[i].rho())));
						this_dtinv = max(this_dtinv, ql[i].max_abs_z_eigen(x, f), qr[i].max_abs_z_eigen(x, f));
					}
				}
				dtinv = max(this_dtinv, dtinv);
			}
		}
		dtinv /= HydroGrid::get_dx();
		if (dtinv == 0.0) {
			dtinv = 1.0E-10;
		}
		this_dt = 1.0 / dtinv;
		this_dt *= GRID_CFL_FACTOR;
		eax = min(this_dt, eax);
		HydroGrid::inc_instruction_pointer(dir);
	}

	void compute_dudt(int dir) {
		const int o = BW - 1;
		const Real hinv = HydroGrid::get_dx() * 0.5;
		Vector<Vector<Real, 3>, 3> fedd, gradv;
		HydroGrid::compute_dudt(dir);
		_3Vec radmom;
		if (dir == 0) {
			for (int k = BW; k < GNX - BW; k++) {
				for (int j = BW; j < GNX - BW; j++) {
					for (int i = BW; i < GNX - BW; i++) {
						fedd = get_eddington_tensor(i - o, j - o, k - o);
						gradv[0][0] = ((*this)(i + 1, j, k).vx() - (*this)(i - 1, j, k).vx()) * hinv;
						gradv[0][1] = ((*this)(i, j + 1, k).vx() - (*this)(i, j - 1, k).vx()) * hinv;
						gradv[0][2] = ((*this)(i, j, k + 1).vx() - (*this)(i, j, k - 1).vx()) * hinv;
						gradv[1][0] = ((*this)(i + 1, j, k).vy() - (*this)(i - 1, j, k).vy()) * hinv;
						gradv[1][1] = ((*this)(i, j + 1, k).vy() - (*this)(i, j - 1, k).vy()) * hinv;
						gradv[1][2] = ((*this)(i, j, k + 1).vy() - (*this)(i, j, k - 1).vy()) * hinv;
						gradv[2][0] = ((*this)(i + 1, j, k).vz() - (*this)(i - 1, j, k).vz()) * hinv;
						gradv[2][1] = ((*this)(i, j + 1, k).vz() - (*this)(i, j - 1, k).vz()) * hinv;
						gradv[2][2] = ((*this)(i, j, k + 1).vz() - (*this)(i, j, k - 1).vz()) * hinv;
						radmom = get_rad_momentum_term(i - o, j - o, k - o);
						Vector<Real, STATE_NF> d = 0.0;
						d[State::sx_index] += radmom[0];
						d[State::sy_index] += radmom[1];
						d[State::sz_index] += radmom[2];
						d[State::et_index] += (*this)(i, j, k).vx() * radmom[0];
						d[State::et_index] += (*this)(i, j, k).vy() * radmom[1];
						d[State::et_index] += (*this)(i, j, k).vz() * radmom[2];
						d[State::er_index] -= gradv[0].dot(fedd[0]) * (*this)(i, j, k).er();
						d[State::er_index] -= gradv[1].dot(fedd[1]) * (*this)(i, j, k).er();
						d[State::er_index] -= gradv[2].dot(fedd[2]) * (*this)(i, j, k).er();
						add_to_dudt(i, j, k, d);
					}
				}
			}
		}
	}
	void load_radiation_grid() {
		//	printf( "loading\n");
		const int o = BW - 1;
		for (int k = BW - 1; k < GNX - BW + 1; k++) {
			for (int j = BW - 1; j < GNX - BW + 1; j++) {
				for (int i = BW - 1; i < GNX - BW + 1; i++) {
					set_erad(i - o, j - o, k - o, (*this)(i, j, k).er());
				}
			}
		}
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					set_eint(i - o, j - o, k - o, (*this)(i, j, k).ei());
					set_cv(i - o, j - o, k - o, (*this)(i, j, k).cv());
					set_rho(i - o, j - o, k - o, (*this)(i, j, k).rho());
				}
			}
		}
	}
	static void setup_grid_structure() {
		Real dt;
		HydroGrid::set_time(0.0);
		OctNode::get_root()->find_local_nodes();
		OctNode::initialize_grids();
		dt = max_dt_driver();
		dt *= min(1.0, MAXINITDT);
		if (MPI_rank() == 0) {
			printf("\t       ");
		}
		for (int l = 0; l <= OctNode::get_max_level_allowed(); l++) {
			//	HydroRadGrid::step(dt);
			HydroGrid::set_time(0.0);
			check_for_refine();
			if (MPI_rank() == 0) {
				printf("%i %i\n", OctNode::get_max_level(), OctNode::get_node_cnt());
			}
			OctNode::initialize_grids();
			dt = min(max_dt_driver(), MAXINITDT);
		}
	}
	void unload_radiation_grid() {
		//	printf( "loading\n");
		const int o = BW - 1;
		Real deint, ei;
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					(*this)(i, j, k).set_er(get_erad(i - o, j - o, k - o));
					ei = (*this)(i, j, k).ei();
					deint = get_eint(i - o, j - o, k - o) - ei;
					(*this)(i, j, k)[State::et_index] += deint;
					(*this)(i, j, k)[State::tau_index] += (1.0 / State::gamma) * pow(ei, 1.0 / State::gamma - 1.0) * deint;
					(*this)(i, j, k)[State::tau_index] = max((*this)(i, j, k)[State::tau_index], pow(State::ei_floor, 1.0 / State::gamma));
				}
			}
		}
	}
	static void step(Real dt) {
		for (int l = 0; l < get_local_node_cnt(); l++) {
			HydroRadGrid* g = dynamic_cast<HydroRadGrid*>(get_local_node(l));
			g->load_radiation_grid();
		}
		HydroGrid::step(dt);
		for (int l = 0; l < get_local_node_cnt(); l++) {
			HydroRadGrid* g = dynamic_cast<HydroRadGrid*>(get_local_node(l));
			g->load_radiation_grid();
		}
		implicit_solve(dt);
		for (int l = 0; l < get_local_node_cnt(); l++) {
			HydroRadGrid* g = dynamic_cast<HydroRadGrid*>(get_local_node(l));
			g->unload_radiation_grid();
		}
	}
	virtual void initialize() {
	}
	virtual void read(FILE* fp) {
		HydroGrid::read(fp);

	}
	virtual void write(FILE* fp) const {
		HydroGrid::write(fp);
	}
	virtual void set_refine_flags() {
		HydroGrid::set_refine_flags();
	}
	virtual Real xf(int i) const {
		return HydroGrid::xf(i);
	}
	virtual Real yf(int i) const {
		return HydroGrid::yf(i);
	}
	virtual Real zf(int i) const {
		return HydroGrid::zf(i);
	}
	virtual void allocate_arrays() {
		Radiation::allocate_arrays();
		HydroGrid::allocate_arrays();
	}
	virtual void deallocate_arrays() {
		Radiation::deallocate_arrays();
		HydroGrid::deallocate_arrays();
	}
	virtual void create_child(const ChildIndex& c) {
		HydroGrid::create_child(c);
		Radiation::create_multigrid_child(c);
	}
	virtual HydroRadGrid* new_octnode() const {
		return new HydroRadGrid;
	}
	virtual Real get_output_point(int i, int j, int k, int l) const {
		return HydroGrid::get_output_point(i, j, k, l);
	}
	virtual int nvar_output() const {
		return HydroGrid::nvar_output();

	}
	virtual const char* output_field_names(int i) const {
		return HydroGrid::output_field_names(i);
	}
protected:
	static void run(int argc, char* argv[]);
public:
	virtual void init() {
		HydroGrid::init();
		MultiGrid::init();
	}
	HydroRadGrid();
	virtual ~HydroRadGrid();
};
#endif
#endif /* HYDRO_RAD_GRID_H_ */
