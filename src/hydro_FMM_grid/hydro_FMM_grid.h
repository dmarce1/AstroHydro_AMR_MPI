/*
 * hydro_FMM_grid.h
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#ifndef HYDRO_FMM_GRID_H_
#define HYDRO_FMM_GRID_H_

#include "../defs.h"

#ifdef USE_FMM

#define USE_FMM_ANGULAR

#include "../hydro_grid/hydro_grid.h"

#include "moment.h"
#include "taylor.h"
#ifdef USE_FMM_ANGULAR
#include "taylor2.h"
#endif

#define FORDER 1
#define FNX (INX+2*FBW)

#define FBW (2*FORDER)

#define FSTAGE 15

class HydroFMMGrid: public HydroGrid {
private:
#ifdef USE_FMM_ANGULAR
	static MPI_Datatype MPI_taylor2_t;
	static MPI_Datatype MPI_comm_taylor2_t[8];
	Taylor2* taylor2_buffer;
	Array3d<Taylor2, FNX, FNX, FNX> RxL;
#endif
	static Real d0_array[INX][INX][INX];
	static Real d1_array[2 * INX + 1][2 * INX + 1][2 * INX + 1][3];
	HydroFMMGrid* neighbors[26];
	typedef void (HydroFMMGrid::*ifunc_t)(int);
	static ifunc_t cs[FSTAGE + 1];
	static MPI_Datatype MPI_send_bnd_t[26];
	static MPI_Datatype MPI_recv_bnd_t[26];
	static MPI_Datatype MPI_comm_child1_t[8];
	static MPI_Datatype MPI_comm_taylor_t[8];
	static MPI_Datatype MPI_comm_child3_t[8];
	static MPI_Datatype MPI_moment_t;
	static MPI_Datatype MPI_taylor_t;
	static MPI_Datatype MPI_4force_t;
	MPI_Request send_request[26], recv_request[26];
	Vector<Real,4>* _4force_buffer;
	Moment* moment_buffer;
	Taylor* taylor_buffer;
	Array3d<Moment, FNX, FNX, FNX> m;
	Array3d<Taylor, FNX, FNX, FNX> L;
	Array3d<Vector<Real, 4>, FNX, FNX, FNX> g;
	Vector<Real, 3> mom_sum;
	bool is_leaf(int i, int j, int k) const;

public:
	static bool solve_on;
	static Vector<Real,6> momentum_sum();
	static Vector<Real, 4> com_sum();
	Real gx(int i, int j, int k) const;
#ifdef USE_FMM_ANGULAR
	Real dlz(int i, int j, int k) const;
#endif
	Real gy(int i, int j, int k) const;
	Real gz(int i, int j, int k) const;
	virtual void compute_dudt(int dir);
	static void step(Real dt);
	Real get_phi(int i, int j, int k) const;
	static void FMM_solve();
	static void MPI_datatypes_init();
	virtual void allocate_arrays();
	virtual void deallocate_arrays();
	virtual HydroFMMGrid* new_octnode() const;
	virtual void initialize();
	HydroFMMGrid();
	virtual ~HydroFMMGrid();
	void find_neighbors();
	void moments_recv(int);
	void moments_recv_wait(int);
	void _4force_recv(int);
	void _4force_recv_wait(int);
	void moments_send(int);
	void moments_send_wait(int);
	void _4force_send(int);
	void _4force_send_wait(int);
	void moments_communicate_all(int);
	void moments_communicate_wait_all(int);
	void compute_interactions(int);
	void expansion_recv(int);
	void expansion_recv_wait(int);
	void expansion_send(int);
	void expansion_send_wait(int);
	void null(int);
	static bool check_for_refine();
	static void to_conserved_energy();
	static void from_conserved_energy();
	static void pot_to_hydro_grid();

	static Real get_phi_at(Real x, Real y, Real z) {
		Real p, tmp;
		HydroFMMGrid* g;
		p = 0.0;
		for (int l = 0; l < get_local_node_cnt(); l++) {
			g = dynamic_cast<HydroFMMGrid*>(get_local_node(l));
			if (MPI_rank() == g->proc()) {
				for (int k = FBW; k < FNX - FBW; k++) {
					if (z >= g->zf(k + BW - FBW) && z < g->zf(k + 1 + BW - FBW)) {
						for (int j = FBW; j < FNX - FBW; j++) {
							if (y >= g->yf(j + BW - FBW) && y < g->yf(j + 1 + BW - FBW)) {
								for (int i = FBW; i < FNX - FBW; i++) {
									if (x >= g->xf(i + BW - FBW) && x < g->xf(i + 1 + BW - FBW)) {
										if (!g->zone_is_refined(i + BW - FBW, j + BW - FBW, k + BW - FBW)) {
											p = g->g(i, j, k)[3];
										}
									}
								}
							}
						}
					}
				}
			}
		}
		tmp = p;
		MPI_Allreduce(&tmp, &p, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
		return p;
	}

};

#endif /* HYDRO_FMM_GRID_H_ */
#endif
