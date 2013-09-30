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

#include "../hydro_grid/hydro_grid.h"

#include "moment.h"
#include "taylor.h"


#define FORDER 2
#define FNX (INX+2*FBW)

#define FBW (2*FORDER)

#define FSTAGE 11

class HydroFMMGrid: public HydroGrid {
private:
	HydroFMMGrid* neighbors[26];
	typedef void (HydroFMMGrid::*ifunc_t)(int);
	static ifunc_t cs[FSTAGE+1];
	static MPI_Datatype MPI_send_bnd_t[26];
	static MPI_Datatype MPI_recv_bnd_t[26];
	static MPI_Datatype MPI_comm_child1_t[8];
	static MPI_Datatype MPI_comm_child2_t[8];
	static MPI_Datatype MPI_moment_t;
	static MPI_Datatype MPI_taylor_t;
	MPI_Request send_request[26], recv_request[26];
	Moment* moment_buffer;
	Taylor* taylor_buffer;
	Array3d<Moment, FNX, FNX, FNX> m;
	Array3d<Taylor, FNX, FNX, FNX> L;
	Array3d<Real, FNX, FNX, FNX> phi;
	Array3d<_3Vec, FNX, FNX, FNX> g;
	Vector<Real,3> mom_sum;
	bool is_leaf(int i, int j, int k) const;

public:
	static void momentum_sum();
	static Vector<Real,4> com_sum();
	Real gx(int i, int j, int k) const;
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
	void moments_send(int);
	void moments_send_wait(int);
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
};

#endif /* HYDRO_FMM_GRID_H_ */
#endif
