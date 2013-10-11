/*
 * hydro_FMM_grid.cpp
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#include "hydro_FMM_grid.h"
#ifdef USE_FMM
#include "../tag.h"
#include <mpi.h>

#define TAG_MOMENT 1
#define TAG_4FORCE 5
#define TAG_BOUND 2
#define TAG_EXPANSION 3
#define TAG_EXPANSION2 4

int l1count = 0;
int l2count = 0;

HydroFMMGrid::ifunc_t HydroFMMGrid::cs[FSTAGE + 1] = { &HydroFMMGrid::moments_recv, &HydroFMMGrid::moments_recv_wait, &HydroFMMGrid::moments_send,
		&HydroFMMGrid::moments_send_wait, &HydroFMMGrid::moments_communicate_all, &HydroFMMGrid::moments_communicate_wait_all, &HydroFMMGrid::compute_interactions,
		&HydroFMMGrid::expansion_recv, &HydroFMMGrid::expansion_recv_wait, &HydroFMMGrid::expansion_send, &HydroFMMGrid::expansion_send_wait,
		&HydroFMMGrid::_4force_recv, &HydroFMMGrid::_4force_recv_wait, &HydroFMMGrid::_4force_send, &HydroFMMGrid::_4force_send_wait, &HydroFMMGrid::null };

MPI_Datatype HydroFMMGrid::MPI_send_bnd_t[26];
MPI_Datatype HydroFMMGrid::MPI_recv_bnd_t[26];
MPI_Datatype HydroFMMGrid::MPI_comm_child1_t[8];
MPI_Datatype HydroFMMGrid::MPI_comm_taylor_t[8];
MPI_Datatype HydroFMMGrid::MPI_comm_child3_t[8];
MPI_Datatype HydroFMMGrid::MPI_4force_t;
MPI_Datatype HydroFMMGrid::MPI_moment_t;
MPI_Datatype HydroFMMGrid::MPI_taylor_t;
#ifdef USE_FMM_ANGULAR
MPI_Datatype HydroFMMGrid::MPI_taylor2_t;
MPI_Datatype HydroFMMGrid::MPI_comm_taylor2_t[8];
#endif
Real HydroFMMGrid::d0_array[INX][INX][INX];
Real HydroFMMGrid::d1_array[2 * INX + 1][2 * INX + 1][2 * INX + 1][3];
bool HydroFMMGrid::solve_on = true;

int face_id[26];
int face_opp_id[26];

void HydroFMMGrid::FMM_solve() {
//	printf("!\n");
	MPI_datatypes_init();
	HydroFMMGrid** list = new HydroFMMGrid*[get_local_node_cnt()];
	for (int i = 0; i < get_local_node_cnt(); i++) {
		list[i] = dynamic_cast<HydroFMMGrid*>(get_local_node(i));
		list[i]->ip = new int[1];
		list[i]->ip[0] = 0;
	}
	int last_ip;
	bool done;
	const int nprocs = get_local_node_cnt();

	do {
		done = true;
		for (int i = 0; i < nprocs; i++) {
			if (list[i]->ip[0] < FSTAGE) {
				do {
					//	printf("%i %i %i\n", list[i]->get_id(), list[i]->get_level(), list[i]->ip[0]);
					last_ip = list[i]->ip[0];
					(list[i]->*cs[list[i]->ip[0]])(0);
				} while (last_ip != list[i]->ip[0] && list[i]->ip[0] <= FSTAGE);
				if (list[i]->ip[0] <= FSTAGE) {
					done = false;
				}
			}
		}
	} while (!done);
	for (int i = 0; i < nprocs; i++) {
		delete[] list[i]->ip;
	}

	delete[] list;
	pot_to_hydro_grid();
//	printf("%i %i\n", l1count, l2count);
//	momentum_sum();
}

HydroFMMGrid::HydroFMMGrid() {
	// TODO Auto-generated constructor stub

}

void HydroFMMGrid::initialize() {
}

HydroFMMGrid* HydroFMMGrid::new_octnode() const {
	return new HydroFMMGrid;
}

void HydroFMMGrid::allocate_arrays() {
	HydroGrid::allocate_arrays();
	m.allocate();
	L.allocate();
	g.allocate();
#ifdef USE_FMM_ANGULAR
	RxL.allocate();
#endif
	for (int k = 0; k < FNX; k++) {
		for (int j = 0; j < FNX; j++) {
			for (int i = 0; i < FNX; i++) {
				L(i, j, k).zero();
#ifdef USE_FMM_ANGULAR
				RxL(i, j, k).zero();
#endif
				m(i, j, k).M() = 0.0;
				m(i, j, k).X = 0.0;
				g(i, j, k) = 0.0;
				m(i, j, k).is_leaf = false;
				for (int a = 0; a < 3; a++) {
					for (int b = a; b < 3; b++) {
						m(i, j, k).M2(a, b) = 0.0;
#ifdef USE_HIGH_ORDER_POT
						for (int c = b; c < 3; c++) {
							m(i, j, k).M3(a, b, c) = 0.0;
						}
#endif
					}
				}
			}
		}
	}
}

void HydroFMMGrid::deallocate_arrays() {
#ifdef USE_FMM_ANGULAR
	RxL.deallocate();
#endif
	HydroGrid::deallocate_arrays();
	m.deallocate();
	L.deallocate();
	g.deallocate();
}

HydroFMMGrid::~HydroFMMGrid() {
// TODO Auto-generated destructor stub
}

void HydroFMMGrid::moments_recv(int) {
	int tag;
	HydroFMMGrid* child;
	const Real dv = pow(get_dx(), 3);
	for (int k = 0; k < FNX; k++) {
		for (int j = 0; j < FNX; j++) {
			for (int i = 0; i < FNX; i++) {
#ifdef USE_FMM_ANGULAR
				RxL(i, j, k).zero();
#endif
				L(i, j, k).zero();
				m(i, j, k).M() = 0.0;
				m(i, j, k).X = 0.0;
				m(i, j, k).is_leaf = false;
				for (int a = 0; a < 3; a++) {
					for (int b = a; b < 3; b++) {
						m(i, j, k).M2(a, b) = 0.0;
#ifdef USE_HIGH_ORDER_POT
						for (int c = b; c < 3; c++) {
							m(i, j, k).M3(a, b, c) = 0.0;
						}
#endif
					}
				}
			}
		}
	}
	for (ChildIndex ci = 0; ci < 8; ci++) {
		child = dynamic_cast<HydroFMMGrid*>(get_child(ci));
		if (child != NULL) {
			tag = tag_gen(TAG_MOMENT, get_id(), ci);
			MPI_Irecv(m.ptr(), 1, MPI_comm_child1_t[ci], get_child(ci)->proc(), tag, MPI_COMM_WORLD, recv_request + ci);
		} else {
			recv_request[ci] = MPI_REQUEST_NULL;
			const int xlb = FBW + ci.get_x() * (INX / 2);
			const int ylb = FBW + ci.get_y() * (INX / 2);
			const int zlb = FBW + ci.get_z() * (INX / 2);
			const int xub = xlb + (INX / 2) - 1;
			const int yub = ylb + (INX / 2) - 1;
			const int zub = zlb + (INX / 2) - 1;
			for (int k = zlb; k <= zub; k++) {
				for (int j = ylb; j <= yub; j++) {
					for (int i = xlb; i <= xub; i++) {
						m(i, j, k).X[0] = HydroGrid::xc(i + BW - FBW);
						m(i, j, k).X[1] = HydroGrid::yc(j + BW - FBW);
						m(i, j, k).X[2] = HydroGrid::zc(k + BW - FBW);
						m(i, j, k).M() = (*this)(i + BW - FBW, j + BW - FBW, k + BW - FBW).rho() * dv;
						m(i, j, k).is_leaf = true;
						for (int n = 0; n < 3; n++) {
							for (int l = n; l < 3; l++) {
								m(i, j, k).M2(n, l) = 0.0;
#ifdef USE_HIGH_ORDER_POT
								for (int a = l; a < 3; a++) {
									m(i, j, k).M3(n, l, a) = 0.0;
								}
#endif
							}
						}
					}
				}
			}
		}
	}
	inc_instruction_pointer();
}

void HydroFMMGrid::_4force_recv(int) {
	int tag;
	HydroFMMGrid* child;
	const Real dv = pow(get_dx(), 3);
	for (ChildIndex ci = 0; ci < 8; ci++) {
		child = dynamic_cast<HydroFMMGrid*>(get_child(ci));
		if (child != NULL) {
			tag = tag_gen(TAG_4FORCE, get_id(), ci);
			MPI_Irecv(g.ptr(), 1, MPI_comm_child3_t[ci], get_child(ci)->proc(), tag, MPI_COMM_WORLD, recv_request + ci);
		} else {
			recv_request[ci] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void HydroFMMGrid::moments_recv_wait(int) {
	int flag;
	MPI_Testall(8, recv_request, &flag, MPI_STATUS_IGNORE );
	if (flag) {
		inc_instruction_pointer();
	}
}

void HydroFMMGrid::_4force_recv_wait(int) {
	int flag;
	MPI_Testall(8, recv_request, &flag, MPI_STATUS_IGNORE );
	if (flag) {
		inc_instruction_pointer();
	}
}

void HydroFMMGrid::moments_send(int) {
	if (get_level() != 0) {
		int cnt;
		_3Vec Y;
		Moment tmp;
		Real mass;
		cnt = 0;
		moment_buffer = new Moment[INX * INX * INX / 8];
		for (int k = FBW; k < FNX - FBW; k += 2) {
			for (int j = FBW; j < FNX - FBW; j += 2) {
				for (int i = FBW; i < FNX - FBW; i += 2) {
					tmp.M() = 0.0;
					for (int p = 0; p < 3; p++) {
						for (int q = p; q < 3; q++) {
							tmp.M2(p, q) = 0.0;
#ifdef USE_HIGH_ORDER_POT
							for (int a = q; a < 3; a++) {
								tmp.M3(p, q, a) = 0.0;
							}
#endif
						}
					}
					tmp.X = 0.0;
					tmp.is_leaf = false;
					for (int a = 0; a < 2; a++) {
						for (int b = 0; b < 2; b++) {
							for (int c = 0; c < 2; c++) {
								mass = m(i + a, j + b, k + c).M();
								tmp.M() += mass;
								tmp.X += (m(i + a, j + b, k + c).X * mass);
							}
						}
					}
					if (tmp.M() != 0.0) {
						tmp.X /= tmp.M();
					} else {
						tmp.X = 0.0;
						for (int a = 0; a < 2; a++) {
							for (int b = 0; b < 2; b++) {
								for (int c = 0; c < 2; c++) {
									tmp.X += (m(i + a, j + b, k + c).X * 0.125);
								}
							}
						}
					}
					for (int a = 0; a < 2; a++) {
						for (int b = 0; b < 2; b++) {
							for (int c = 0; c < 2; c++) {
								for (int p = 0; p < 3; p++) {
									for (int q = p; q < 3; q++) {
										Moment& child = m(i + a, j + b, k + c);
										Y = child.X - tmp.X;
										tmp.M2(p, q) += child.M2(p, q);
										tmp.M2(p, q) += child.M() * Y[p] * Y[q];
#ifdef USE_HIGH_ORDER_POT
										for (int r = q; r < 3; r++) {
											tmp.M3(p, q, r) += child.M() * Y[p] * Y[q] * Y[r];
											tmp.M3(p, q, r) += child.M3(p, q, r);
											tmp.M3(p, q, r) += child.M2(p, q) * Y[r];
											tmp.M3(p, q, r) += child.M2(q, r) * Y[p];
											tmp.M3(p, q, r) += child.M2(r, p) * Y[q];
										}
#endif

									}
								}
							}
						}
					}
					moment_buffer[cnt] = tmp;
					cnt++;
				}
			}
		}
//		printf("%i\n", cnt);
		assert(cnt==INX * INX * INX / 8);
		int tag = tag_gen(TAG_MOMENT, get_parent()->get_id(), my_child_index());
		MPI_Isend(moment_buffer, cnt * sizeof(Moment), MPI_BYTE, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
	}
	inc_instruction_pointer();
}

void HydroFMMGrid::_4force_send(int) {
	if (get_level() != 0) {
		int cnt;
		_3Vec Y;
		Vector<Real, 4> tmp;
		Real mass, mtot;
		cnt = 0;
		_4force_buffer = new Vector<Real, 4> [INX * INX * INX / 8];
		for (int k = FBW; k < FNX - FBW; k += 2) {
			for (int j = FBW; j < FNX - FBW; j += 2) {
				for (int i = FBW; i < FNX - FBW; i += 2) {
					tmp = 0.0;
					mtot = 0.0;
					for (int a = 0; a < 2; a++) {
						for (int b = 0; b < 2; b++) {
							for (int c = 0; c < 2; c++) {
								//				mass = m(i + a, j + b, k + c).M();
								//				mtot += mass;
								tmp += g(i + a, j + b, k + c); // * mass;
							}
						}
					}
					tmp /= 8.0; // * mtot;
					_4force_buffer[cnt] = tmp;
					cnt++;
				}
			}
		}
//		printf("%i\n", cnt);
		assert(cnt==INX * INX * INX / 8);
		int tag = tag_gen(TAG_4FORCE, get_parent()->get_id(), my_child_index());
		MPI_Isend(_4force_buffer, cnt * sizeof(Vector<Real, 4> ), MPI_BYTE, get_parent()->proc(), tag, MPI_COMM_WORLD, send_request);
	}
	inc_instruction_pointer();
}

void HydroFMMGrid::moments_send_wait(int) {
	int flag;
	if (get_level() != 0) {
		MPI_Test(send_request, &flag, MPI_STATUS_IGNORE );
		if (flag) {
			delete[] moment_buffer;
			inc_instruction_pointer();
		}
	} else {
		inc_instruction_pointer();
	}
}
void HydroFMMGrid::_4force_send_wait(int) {
	int flag;
	if (get_level() != 0) {
		MPI_Test(send_request, &flag, MPI_STATUS_IGNORE );
		if (flag) {
			delete[] _4force_buffer;
			inc_instruction_pointer();
		}
	} else {
		inc_instruction_pointer();
	}
}

void HydroFMMGrid::moments_communicate_all(int) {
	int dir, tag_send, tag_recv;
	OctNode* neighbor;
	find_neighbors();
	for (dir = 0; dir < 26; dir++) {
		neighbor = neighbors[dir];
		if (neighbor == NULL) {
			recv_request[dir] = MPI_REQUEST_NULL;
			send_request[dir] = MPI_REQUEST_NULL;
		} else {
			assert(dir==face_id[dir]);
			tag_send = tag_gen(TAG_BOUND, neighbor->get_id(), face_opp_id[dir]);
			tag_recv = tag_gen(TAG_BOUND, get_id(), face_id[dir]);
			MPI_Isend(m.ptr(), 1, MPI_send_bnd_t[dir], neighbor->proc(), tag_send, MPI_COMM_WORLD, send_request + dir);
			MPI_Irecv(m.ptr(), 1, MPI_recv_bnd_t[dir], neighbor->proc(), tag_recv, MPI_COMM_WORLD, recv_request + dir);
		}
	}
	inc_instruction_pointer();
}

void HydroFMMGrid::moments_communicate_wait_all(int) {
	int flag_recv, flag_send;
	bool ready;
	MPI_Testall(26, recv_request, &flag_recv, MPI_STATUS_IGNORE );
	MPI_Testall(26, send_request, &flag_send, MPI_STATUS_IGNORE );
	ready = flag_recv && flag_send;
	if (ready) {
		inc_instruction_pointer();
	}
}

void HydroFMMGrid::compute_interactions(int) {
	Moment *n1;
	Moment *n2;
	Taylor *l;
	int xlb, xub, ylb, yub, zlb, zub;
	static const Real delta[3][3] = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
#ifdef USE_FMM_ANGULAR
	Taylor2 *rxl;
//	static const Real eps[3][3][3] = { { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, +1.0 }, { 0.0, -1.0, 0.0 } },
//			{ { 0.0, 0.0, -1.0 }, { 0.0, 0.0, 0.0 }, { +1.0, 0.0, 0.0 } }, { { 0.0, +1.0, 0.0 }, { -1.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } } };
#endif
	const Real dxinv = 1.0 / get_dx();
	for (int k0 = FBW; k0 < FNX - FBW; k0++) {
		for (int j0 = FBW; j0 < FNX - FBW; j0++) {
			for (int i0 = FBW; i0 < FNX - FBW; i0++) {
				if (get_level() != 0) {
					xlb = 2 * ((i0 / 2) - FORDER);
					ylb = 2 * ((j0 / 2) - FORDER);
					zlb = 2 * ((k0 / 2) - FORDER);
					xub = 2 * ((i0 / 2) + FORDER) + 1;
					yub = 2 * ((j0 / 2) + FORDER) + 1;
					zub = 2 * ((k0 / 2) + FORDER) + 1;
				} else {
					xlb = ylb = zlb = FBW;
					xub = yub = zub = FNX - FBW - 1;
				}
				n1 = m.ptr(i0, j0, k0);
				l = L.ptr(i0, j0, k0);
				l->zero();
				Real tmp;
				//Self potential
			//	if (n1->is_leaf) {
			//		const Real factor = (3.0 + 3.0 / sqrt(2.0) + 1.0 / sqrt(3.0)) / 3.0;
			//		(*l)() -= factor * n1->M() / get_dx();
			//	}
				for (int k = zlb; k <= zub; k++) {
					for (int j = ylb; j <= yub; j++) {
						for (int i = xlb; i <= xub; i++) {
							n2 = m.ptr(i, j, k);
							bool interaction_pair;
							if (n1->is_leaf || n2->is_leaf) {
								interaction_pair = !((i == i0) && (j == j0) && (k == k0));
							} else {
								interaction_pair = !((abs(i0 - i) <= FORDER) && (abs(j - j0) <= FORDER) && (abs(k - k0) <= FORDER));
							}
							if (interaction_pair) {
								if (!n1->is_leaf || !n2->is_leaf) {
									const _3Vec Y = n1->X - n2->X;
									const Real r2inv = 1.0 / Y.dot(Y);
									const Real d0 = -sqrt(r2inv);
									const Real d1 = -d0 * r2inv;
									const Real d2 = -3.0 * d1 * r2inv;
									const Real d3 = -5.0 * d2 * r2inv;
									(*l)() += n2->M() * d0;
									for (int a = 0; a < 3; a++) {
										(*l)() += n2->M2(a, a) * d1 * 0.5;
										(*l)(a) += n2->M() * d1 * Y[a];
										for (int b = 0; b < 3; b++) {
											(*l)() += n2->M2(a, b) * d2 * 0.5 * Y[a] * Y[b];
											(*l)(a) += d2 * (n2->M2(b, b) * 0.5 * Y[a] + n2->M2(a, b) * Y[b]);
											for (int c = 0; c < 3; c++) {
#ifdef USE_HIGH_ORDER_POT
												(*l)() += n2->M3(a, b, c) * (1.0 / 6.0)
														* (d3 * Y[a] * Y[b] * Y[c] + d2 * (Y[a] * delta[b][c] + Y[c] * delta[a][b] + Y[b] * delta[c][a]));
#endif
												(*l)(a) += n2->M2(b, c) * d3 * 0.5 * Y[a] * Y[b] * Y[c];
											}
										}
										for (int b = a; b < 3; b++) {
											(*l)(a, b) += n2->M() * (delta[a][b] * d1 + Y[a] * Y[b] * d2);
											for (int c = b; c < 3; c++) {
												(*l)(a, b, c) += n2->M() * (d2 * (delta[a][b] * Y[c] + delta[b][c] * Y[a] + delta[a][c] * Y[b]) + d3 * Y[a] * Y[b] * Y[c]);
											}
										}
									}
								} else {
									Real tmp = n2->M() * dxinv;
									(*l)() += d0_array[abs(k - k0)][abs(j - j0)][abs(i - i0)] * tmp;
									tmp *= dxinv;
									const int k1 = k - k0 + INX;
									const int j1 = j - j0 + INX;
									const int i1 = i - i0 + INX;
									(*l)(0) += d1_array[k1][j1][i1][0] * tmp;
									(*l)(1) += d1_array[k1][j1][i1][1] * tmp;
									(*l)(2) += d1_array[k1][j1][i1][2] * tmp;
								}
							}
						}
					}
				}
#ifdef USE_FMM_ANGULAR
				rxl = RxL.ptr(i0, j0, k0);
				rxl->zero();
				const _3Vec R = n1->X;
				for (int a = 0; a < 2; a++) {
					const int b = 1 - a;
					const Real eps = Real(2 * b - 1);
					(*rxl)() -= eps * R[a] * (*l)(b);
					for (int m = 0; m < 3; m++) {
						(*rxl)(m) -= eps * (+delta[a][m] * (*l)(b) + R[a] * (*l)(b, m));
						for (int n = 0; n <= m; n++) {
							(*rxl)(m, n) -= eps * (delta[a][n] * (*l)(b, m) + delta[a][m] * (*l)(b, n) + R[a] * (*l)(b, n, m));
						}
					}
				}
#endif
			}
		}
	}
	inc_instruction_pointer();
}

void HydroFMMGrid::expansion_recv(int) {
	if (get_level() != 0) {
		const ChildIndex ci = my_child_index();
		const HydroFMMGrid* p = dynamic_cast<const HydroFMMGrid*>(get_parent());
		int i, j, k, sz, tag;
		Real a;
		tag = tag_gen(TAG_EXPANSION, get_id(), ci);
		taylor_buffer = new Taylor[INX * INX * INX / 8];
		MPI_Irecv(taylor_buffer, INX * INX * INX * sizeof(Taylor) / 8, MPI_BYTE, p->proc(), tag, MPI_COMM_WORLD, recv_request + 0);
#ifdef USE_FMM_ANGULAR
		taylor2_buffer = new Taylor2[INX * INX * INX / 8];
		tag = tag_gen(TAG_EXPANSION2, get_id(), ci);
		MPI_Irecv(taylor2_buffer, INX * INX * INX * sizeof(Taylor2) / 8, MPI_BYTE, p->proc(), tag, MPI_COMM_WORLD, recv_request + 1);
#endif
	}
	inc_instruction_pointer();

}

void HydroFMMGrid::expansion_recv_wait(int) {
	bool rc;
	int cnt;
#ifdef USE_FMM_ANGULAR
	int cnt2;
#endif
	if (get_level() != 0) {
		int flag;
#ifdef USE_FMM_ANGULAR
		MPI_Testall(2, recv_request, &flag, MPI_STATUS_IGNORE );
#else
		MPI_Test(recv_request, &flag, MPI_STATUS_IGNORE );
#endif
		rc = flag;
		if (rc) {
			Taylor pL;
			Taylor cL;
			_3Vec Y, Z;
			Real mtmp;
			cnt = 0;
#ifdef USE_FMM_ANGULAR
			cnt2 = 0;
			Taylor2 RxpL;
			Taylor2 RxcL;
#endif
			for (int k0 = FBW; k0 < FNX - FBW; k0 += 2) {
				for (int j0 = FBW; j0 < FNX - FBW; j0 += 2) {
					for (int i0 = FBW; i0 < FNX - FBW; i0 += 2) {
#ifdef USE_FMM_ANGULAR
						RxpL = taylor2_buffer[cnt2];
						cnt2++;
#endif
						pL = taylor_buffer[cnt];
						cnt++;
						Z = 0.0;
						mtmp = 0.0;
						for (int k = k0; k < k0 + 2; k++) {
							for (int j = j0; j < j0 + 2; j++) {
								for (int i = i0; i < i0 + 2; i++) {
									mtmp += m(i, j, k).M();
									Z += (m(i, j, k).X * m(i, j, k).M());
								}
							}
						}
						if (mtmp > 0.0) {
							Z /= mtmp;
						} else {
							Z[0] = m(i0, j0, k0).X[0] + get_dx() / 2.0;
							Z[1] = m(i0, j0, k0).X[1] + get_dx() / 2.0;
							Z[2] = m(i0, j0, k0).X[2] + get_dx() / 2.0;
						}
						for (int k = k0; k < k0 + 2; k++) {
							for (int j = j0; j < j0 + 2; j++) {
								for (int i = i0; i < i0 + 2; i++) {
									Y = m(i, j, k).X - Z;
									cL.zero();
									cL() += pL();
									for (int a = 0; a < 3; a++) {
										cL() += pL(a) * Y[a];
										cL(a) += pL(a);
										for (int b = 0; b < 3; b++) {
											cL() += pL(a, b) * Y[a] * Y[b] * 0.5;
											cL(a) += pL(a, b) * Y[b];
											for (int c = 0; c < 3; c++) {
												cL() += pL(a, b, c) * Y[a] * Y[b] * Y[c] * (1.0 / 6.0);
												cL(a) += pL(a, b, c) * Y[b] * Y[c] * 0.5;
											}
										}
										for (int b = a; b < 3; b++) {
											cL(a, b) += pL(a, b);
											for (int c = 0; c < 3; c++) {
												cL(a, b) += pL(a, b, c) * Y[c];
											}
											for (int c = b; c < 3; c++) {
												cL(a, b, c) += pL(a, b, c);
											}
										}
									}
									L(i, j, k) += cL;
#ifdef USE_FMM_ANGULAR
									RxcL.zero();
									RxcL() += RxpL();
									for (int a = 0; a < 3; a++) {
										RxcL() += RxpL(a) * Y[a];
										RxcL(a) += RxpL(a);
										for (int b = 0; b < 3; b++) {
											RxcL() += RxpL(a, b) * Y[a] * Y[b] * 0.5;
											RxcL(a) += RxpL(a, b) * Y[b];
											for (int c = 0; c < 3; c++) {
											}
										}
										for (int b = a; b < 3; b++) {
											RxcL(a, b) += RxpL(a, b);
										}
									}
									RxL(i, j, k) += RxcL;
#endif
								}
							}
						}
					}
				}
			}
			delete[] taylor_buffer;
#ifdef USE_FMM_ANGULAR
			delete[] taylor2_buffer;
#endif
		}
	} else {
		rc = true;
	}
	if (rc) {
		inc_instruction_pointer();
	}
}

void HydroFMMGrid::expansion_send(int) {
	HydroFMMGrid* child;
	int tag;
	for (int ci = 0; ci < 8; ci++) {
		child = dynamic_cast<HydroFMMGrid*>(get_child(ci));
		if (child == NULL) {
			send_request[ci] = MPI_REQUEST_NULL;
#ifdef USE_FMM_ANGULAR
			send_request[ci + 8] = MPI_REQUEST_NULL;
#endif
		} else {
			tag = tag_gen(TAG_EXPANSION, child->get_id(), ci);
			MPI_Isend(L.ptr(), 1, MPI_comm_taylor_t[ci], child->proc(), tag, MPI_COMM_WORLD, send_request + ci);
#ifdef USE_FMM_ANGULAR
			tag = tag_gen(TAG_EXPANSION2, child->get_id(), ci);
			MPI_Isend(RxL.ptr(), 1, MPI_comm_taylor2_t[ci], child->proc(), tag, MPI_COMM_WORLD, send_request + ci + 8);
#endif
		}
	}
	for (int k = FBW; k < FNX - FBW; k++) {
		for (int j = FBW; j < FNX - FBW; j++) {
			for (int i = FBW; i < FNX - FBW; i++) {
				for (int q = 0; q < 3; q++) {
					g(i, j, k)[q] = -L(i, j, k)(q) * PhysicalConstants::G;
				}
				g(i, j, k)[3] = L(i, j, k)() * PhysicalConstants::G;
			}
		}
	}
	inc_instruction_pointer();

}

void HydroFMMGrid::expansion_send_wait(int) {
	int flag;
#ifdef USE_FMM_ANGULAR
	MPI_Testall(16, send_request, &flag, MPI_STATUS_IGNORE );
#else
	MPI_Testall(8, send_request, &flag, MPI_STATUS_IGNORE );
#endif
	if (flag) {
		inc_instruction_pointer();
	}
}

Real HydroFMMGrid::get_phi(int i, int j, int k) const {
	return g(i + FBW - BW, j + FBW - BW, k + FBW - BW)[3];
}

Real HydroFMMGrid::gx(int i, int j, int k) const {
	return g(i + FBW - BW, j + FBW - BW, k + FBW - BW)[0];
}

Real HydroFMMGrid::gy(int i, int j, int k) const {
	return g(i + FBW - BW, j + FBW - BW, k + FBW - BW)[1];
}

Real HydroFMMGrid::gz(int i, int j, int k) const {
	return g(i + FBW - BW, j + FBW - BW, k + FBW - BW)[2];
}

void HydroFMMGrid::null(int) {
	inc_instruction_pointer();
}

void HydroFMMGrid::step(Real dt) {
	Real start_time;
	Real beta[3] = { 1.0, 0.25, 2.0 / 3.0 };
	HydroGrid::set_dt(dt);
	store();
	for (int i = 0; i < 3; i++) {
		HydroGrid::set_beta(beta[i]);
		start_time = MPI_Wtime();
		substep_driver();
		FMM_solve();
	}
	set_time(get_time() + dt);
}

void HydroFMMGrid::compute_dudt(int dir) {
	HydroGrid::compute_dudt(dir);
	if (dir == 0) {
		_3Vec x;
		Real d;
		Vector<Real, STATE_NF> D, D0;
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					x = HydroGrid::X(i, j, k);
					d = (*this)(i, j, k).rho();
					D = 0.0;
					D0 = get_dudt(i, j, k);
					if (State::cylindrical) {
						Real R = sqrt(x[0] * x[0] + x[1] * x[1]);
#ifdef USE_FMM
						D[State::sy_index] = d * dlz(i, j, k);
						//			printf( "%e %e\n", d * (x[0] * gy(i, j, k) - x[1] * gx(i, j, k)), d*dlz(i, j, k) );
#else
						D[State::sy_index] = d * (x[0] * gy(i, j, k) - x[1] * gx(i, j, k));
#endif
						D[State::sx_index] = d * (x[0] * gx(i, j, k) + x[1] * gy(i, j, k)) / R;
					} else {
						D[State::sx_index] = d * gx(i, j, k);
						D[State::sy_index] = d * gy(i, j, k);
					}
					D[State::sz_index] = d * gz(i, j, k);
					D[State::et_index] = D0[State::pot_index] - (*this)(i, j, k).pot() / d * D0[State::d_index];
					//				D[State::pot_index] = -D0[State::pot_index];
					add_to_dudt(i, j, k, D);
				}
			}
		}
	}
}

bool HydroFMMGrid::check_for_refine() {
	bool rc;
	to_conserved_energy();
	rc = OctNode::check_for_refine();
	if (rc) {
		pot_to_hydro_grid();
		HydroGrid::redistribute_grids();
		int count = OctNode::get_local_node_cnt();
		FMM_solve();
	}
	from_conserved_energy();
	return rc;

}

void HydroFMMGrid::to_conserved_energy() {
	HydroGrid* g0;
	for (int n = 0; n < get_local_node_cnt(); n++) {
		g0 = dynamic_cast<HydroGrid*>(get_local_node(n));
		for (int k = BW - 1; k < GNX - BW + 1; k++) {
			for (int j = BW - 1; j < GNX - BW + 1; j++) {
				for (int i = BW - 1; i < GNX - BW + 1; i++) {
					(*g0)(i, j, k).to_con(g0->X(i, j, k));
				}
			}
		}
	}
}

void HydroFMMGrid::from_conserved_energy() {
	HydroGrid* g0;
	for (int n = 0; n < get_local_node_cnt(); n++) {
		g0 = dynamic_cast<HydroGrid*>(get_local_node(n));
		for (int k = BW - 1; k < GNX - BW + 1; k++) {
			for (int j = BW - 1; j < GNX - BW + 1; j++) {
				for (int i = BW - 1; i < GNX - BW + 1; i++) {
					(*g0)(i, j, k).from_con(g0->X(i, j, k));
				}
			}
		}
	}
}

Vector<Real, 6> HydroFMMGrid::momentum_sum() {
	HydroFMMGrid* g0;
	Real sum[3];
	Real lz, lz2, sx2, sy2, R, x, y;
	lz = 0.0;
	lz2 = 0.0;
	sx2 = sy2 = 0.0;
	sum[0] = sum[1] = sum[2] = 0.0;
	Real norm = 0.0;
	Real norml = 0.0, this_sr, this_lz;
	for (int n = 0; n < get_local_node_cnt(); n++) {
		g0 = dynamic_cast<HydroFMMGrid*>(get_local_node(n));
		for (int k = FBW; k < FNX - FBW; k++) {
			for (int j = FBW; j < FNX - FBW; j++) {
				for (int i = FBW; i < FNX - FBW; i++) {
					if (g0->m(i, j, k).is_leaf) {
#ifdef USE_FMM_ANGULAR
						x = g0->m(i, j, k).X[0];
						y = g0->m(i, j, k).X[1];
						R = sqrt(x * x + y * y);
						this_lz = g0->dlz(i + BW - FBW, j + BW - FBW, k + BW - FBW) * g0->m(i, j, k).M();
						lz2 += this_lz;
						this_sr = (g0->g(i, j, k)[0] * x + g0->g(i, j, k)[1] * y) * g0->m(i, j, k).M();
						this_sr /= R;
						sx2 += x / R * this_sr - y / R / R * this_lz;
						sy2 += y / R * this_sr + x / R / R * this_lz;
#endif
						norml += fabs(g0->dlz(i + BW - FBW, j + BW - FBW, k + BW - FBW) * g0->m(i, j, k).M());
						lz += (g0->g(i, j, k)[1] * g0->m(i, j, k).X[0] - g0->g(i, j, k)[0] * g0->m(i, j, k).X[1]) * g0->m(i, j, k).M();
						sum[0] += g0->g(i, j, k)[0] * g0->m(i, j, k).M();
						sum[1] += g0->g(i, j, k)[1] * g0->m(i, j, k).M();
						sum[2] += g0->g(i, j, k)[2] * g0->m(i, j, k).M();
						norm += sqrt((g0->g(i, j, k)).dot(g0->g(i, j, k))) * g0->m(i, j, k).M();
					}
				}
			}
		}
	}
	Real sum0[3];
	Real norm0 = norm;
	Real lz0 = lz;
	sum0[0] = sum[0];
	sum0[1] = sum[1];
	sum0[2] = sum[2];
	MPI_Allreduce(sum0, sum, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce(&norm0, &norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce(&lz0, &lz, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	lz0 = sx2;
	MPI_Allreduce(&lz0, &sx2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	lz0 = sy2;
	MPI_Allreduce(&lz0, &sy2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	lz0 = lz2;
	MPI_Allreduce(&lz0, &lz2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	norm0 = norml;
	MPI_Allreduce(&norm0, &norml, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	if (MPI_rank() == 0) {
		printf("Momentum Sum = %e %e %e %e %e\n", sum[0] / norm, sum[1] / norm, sum[2] / norm, lz, lz2);
	}
	Vector<Real, 6> ret_vec;
	ret_vec[0] = lz / norml;
	ret_vec[1] = lz2 / norml;
	ret_vec[2] = sum[0] / norm;
	ret_vec[3] = sx2 / norm;
	ret_vec[4] = sum[1] / norm;
	ret_vec[5] = sy2 / norm;
	return ret_vec;
}

Vector<Real, 4> HydroFMMGrid::com_sum() {
	HydroFMMGrid* g0;
	Real sum[4];
	sum[0] = sum[1] = sum[2] = sum[3] = 0.0;
	Real norm = 0.0;
	for (int n = 0; n < get_local_node_cnt(); n++) {
		g0 = dynamic_cast<HydroFMMGrid*>(get_local_node(n));
		const Real dv = pow(g0->get_dx(), 3);
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!g0->zone_is_refined(i, j, k)) {
						sum[0] += (*g0)(i, j, k).rho() * dv * g0->xc(i);
						sum[1] += (*g0)(i, j, k).rho() * dv * g0->yc(j);
						sum[2] += (*g0)(i, j, k).rho() * dv * g0->zc(k);
						sum[3] += (*g0)(i, j, k).rho() * dv;
					}
				}
			}
		}
	}
	Real sum0[4];
	sum0[0] = sum[0];
	sum0[1] = sum[1];
	sum0[2] = sum[2];
	sum0[3] = sum[3];
	MPI_Allreduce(sum0, sum, 4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	sum[0] /= sum[3];
	sum[1] /= sum[3];
	sum[2] /= sum[3];
	if (MPI_rank() == 0) {
		FILE* fp = fopen("com.dat", "at");
		fprintf(fp, "%e %.14e %.14e %.14e  %.14e\n", get_time(), sum[0], sum[1], sum[2], sum[3]);
		fclose(fp);
	}
	Vector<Real, 4> v;
	v[0] = sum[0];
	v[1] = sum[1];
	v[2] = sum[2];
	v[3] = sum[3];
	return v;
}

bool HydroFMMGrid::is_leaf(int i, int j, int k) const {
	const HydroFMMGrid* ptr = this;
	for (int l = 0; l < 3; l++) {
		if (i < FBW) {
			if (ptr->get_sibling(XL) != NULL) {
				ptr = dynamic_cast<const HydroFMMGrid*>(ptr->get_sibling(XL));
				i += INX;
			}
		} else if (i >= FNX - FBW) {
			if (ptr->get_sibling(XU) != NULL) {
				ptr = dynamic_cast<const HydroFMMGrid*>(ptr->get_sibling(XU));
				i -= INX;
			}
		}
		if (j < FBW) {
			if (ptr->get_sibling(YL) != NULL) {
				ptr = dynamic_cast<const HydroFMMGrid*>(ptr->get_sibling(YL));
				j += INX;
			}
		} else if (j >= FNX - FBW) {
			if (ptr->get_sibling(YU) != NULL) {
				ptr = dynamic_cast<const HydroFMMGrid*>(ptr->get_sibling(YU));
				j -= INX;
			}
		}
		if (k < FBW) {
			if (ptr->get_sibling(ZL) != NULL) {
				ptr = dynamic_cast<const HydroFMMGrid*>(ptr->get_sibling(ZL));
				k += INX;
			}
		} else if (k >= FNX - FBW) {
			if (ptr->get_sibling(ZU) != NULL) {
				ptr = dynamic_cast<const HydroFMMGrid*>(ptr->get_sibling(ZU));
				k -= INX;
			}
		}
	}
	return !ptr->zone_is_refined(i + BW - FBW, j + BW - FBW, k + BW - FBW);
}

void HydroFMMGrid::find_neighbors() {
	OctNode* corners[8];
	OctNode* edges[12];
	for (int i = 0; i < 12; i++) {
		edges[i] = NULL;
	}
	for (int i = 0; i < 8; i++) {
		corners[i] = NULL;
	}
	if (get_sibling(XL) != NULL) {
		edges[0] = get_sibling(XL)->get_sibling(YL);
		edges[1] = get_sibling(XL)->get_sibling(YU);
		edges[2] = get_sibling(XL)->get_sibling(ZL);
		edges[3] = get_sibling(XL)->get_sibling(ZU);
	}
	if (get_sibling(XU) != NULL) {
		edges[4] = get_sibling(XU)->get_sibling(YL);
		edges[5] = get_sibling(XU)->get_sibling(YU);
		edges[6] = get_sibling(XU)->get_sibling(ZL);
		edges[7] = get_sibling(XU)->get_sibling(ZU);
	}
	if (get_sibling(YL) != NULL) {
		edges[0] = get_sibling(YL)->get_sibling(XL);
		edges[4] = get_sibling(YL)->get_sibling(XU);
		edges[8] = get_sibling(YL)->get_sibling(ZL);
		edges[9] = get_sibling(YL)->get_sibling(ZU);
	}
	if (get_sibling(YU) != NULL) {
		edges[1] = get_sibling(YU)->get_sibling(XL);
		edges[5] = get_sibling(YU)->get_sibling(XU);
		edges[10] = get_sibling(YU)->get_sibling(ZL);
		edges[11] = get_sibling(YU)->get_sibling(ZU);
	}
	if (get_sibling(ZL) != NULL) {
		edges[2] = get_sibling(ZL)->get_sibling(XL);
		edges[6] = get_sibling(ZL)->get_sibling(XU);
		edges[8] = get_sibling(ZL)->get_sibling(YL);
		edges[10] = get_sibling(ZL)->get_sibling(YU);
	}
	if (get_sibling(ZU) != NULL) {
		edges[3] = get_sibling(ZU)->get_sibling(XL);
		edges[7] = get_sibling(ZU)->get_sibling(XU);
		edges[9] = get_sibling(ZU)->get_sibling(YL);
		edges[11] = get_sibling(ZU)->get_sibling(YU);
	}
	for (int i = 0; i < 8; i++) {
		int l = (i & 1) + ZL;
		int k = ((i & 2) >> 1) + YL;
		int j = ((i & 4) >> 2) + XL;
		corners[i] = NULL;
		if (get_sibling(j) != NULL) {
			if (get_sibling(j)->get_sibling(k) != NULL) {
				corners[i] = get_sibling(j)->get_sibling(k)->get_sibling(l);
			} else if (get_sibling(j)->get_sibling(l) != NULL) {
				corners[i] = get_sibling(j)->get_sibling(l)->get_sibling(k);
			}
		} else if (get_sibling(k) != NULL) {
			if (get_sibling(k)->get_sibling(j) != NULL) {
				corners[i] = get_sibling(k)->get_sibling(j)->get_sibling(l);
			} else if (get_sibling(k)->get_sibling(l) != NULL) {
				corners[i] = get_sibling(k)->get_sibling(l)->get_sibling(j);
			}
		} else if (get_sibling(l) != NULL) {
			if (get_sibling(l)->get_sibling(j) != NULL) {
				corners[i] = get_sibling(l)->get_sibling(j)->get_sibling(k);
			} else if (get_sibling(l)->get_sibling(k) != NULL) {
				corners[i] = get_sibling(l)->get_sibling(k)->get_sibling(j);
			}
		}
	}
	for (int i = 0; i < 6; i++) {
		if (get_sibling(i) != NULL) {
			neighbors[i] = dynamic_cast<HydroFMMGrid*>(get_sibling(i));
		} else {
			neighbors[i] = NULL;
		}
	}
	for (int i = 0; i < 12; i++) {
		if (edges[i] != NULL) {
			neighbors[i + 6] = dynamic_cast<HydroFMMGrid*>(edges[i]);
		} else {
			neighbors[i + 6] = NULL;
		}
	}
	for (int i = 0; i < 8; i++) {
		if (corners[i] != NULL) {
			neighbors[i + 18] = dynamic_cast<HydroFMMGrid*>(corners[i]);
		} else {
			neighbors[i + 18] = NULL;
		}
	}
}

#ifdef USE_FMM_ANGULAR
Real HydroFMMGrid::dlz(int i, int j, int k) const {
	return RxL(i + FBW - BW, j + FBW - BW, k + FBW - BW)() * PhysicalConstants::G;
}
#endif

void HydroFMMGrid::pot_to_hydro_grid() {
	HydroFMMGrid* p;
	HydroGrid* g;
	Real pot;
	for (int n = 0; n < get_local_node_cnt(); n++) {
		p = dynamic_cast<HydroFMMGrid*>(get_local_node(n));
		g = dynamic_cast<HydroGrid*>(get_local_node(n));
		for (int k = BW - 1; k < GNX - BW + 1; k++) {
			for (int j = BW - 1; j < GNX - BW + 1; j++) {
				for (int i = BW - 1; i < GNX - BW + 1; i++) {
					pot = p->get_phi(i, j, k);
					pot *= (*g)(i, j, k).rho();
					pot += (*g)(i, j, k).rot_pot(g->X(i, j, k));
					(*g)(i, j, k).set_pot(pot);
				}
			}
		}
	}
}

void HydroFMMGrid::MPI_datatypes_init() {
	static bool initialized = false;
	if (!initialized) {
		for (int k = 0; k < INX; k++) {
			for (int j = 0; j < INX; j++) {
				for (int i = 0; i < INX; i++) {
					if (!(i == 0 && j == 0 && k == 0)) {
						d0_array[k][j][i] = -1.0 / sqrt(i * i + j * j + k * k);
					}
				}
			}
		}
		for (int k = -INX; k < INX; k++) {
			for (int j = -INX; j < INX; j++) {
				for (int i = -INX; i < INX; i++) {
					if (!(i == 0 && j == 0 && k == 0)) {
						d1_array[k + INX][j + INX][i + INX][0] = -i / pow(i * i + j * j + k * k, 1.5);
						d1_array[k + INX][j + INX][i + INX][1] = -j / pow(i * i + j * j + k * k, 1.5);
						d1_array[k + INX][j + INX][i + INX][2] = -k / pow(i * i + j * j + k * k, 1.5);
					}
				}
			}
		}
		const int lbi0 = FBW;
		const int ubi0 = 2 * FBW - 1;
		const int lbi1 = FNX - 2 * FBW;
		const int ubi1 = FNX - FBW - 1;
		const int lbe0 = 0;
		const int ube0 = FBW - 1;
		const int lbe1 = FNX - FBW;
		const int ube1 = FNX - 1;
		const int hbn0 = FNX - 1 - FBW;
		const int lbn0 = FBW;
		const int hbn1 = FNX - 1 - FBW;
		const int lbn1 = FBW;
		MPI_Type_contiguous(sizeof(Vector<Real, 4> ), MPI_BYTE, &MPI_4force_t);
		MPI_Type_commit(&MPI_4force_t);
		MPI_Type_contiguous(sizeof(Moment), MPI_BYTE, &MPI_moment_t);
		MPI_Type_commit(&MPI_moment_t);
		MPI_Type_contiguous(sizeof(Taylor), MPI_BYTE, &MPI_taylor_t);
		MPI_Type_commit(&MPI_taylor_t);
#ifdef USE_FMM_ANGULAR
		MPI_Type_contiguous(sizeof(Taylor2), MPI_BYTE, &MPI_taylor2_t);
		MPI_Type_commit(&MPI_taylor2_t);
#endif

#define XLYL 0
#define XLYU 1
#define XLZL 2
#define XLZU 3
#define XUYL 4
#define XUYU 5
#define XUZL 6
#define XUZU 7
#define YLZL 8
#define YLZU 9
#define YUZL 10
#define YUZU 11
#define XLYLZL 0
#define XLYLZU 1
#define XLYUZL 2
#define XLYUZU 3
#define XUYLZL 4
#define XUYLZU 5
#define XUYUZL 6
#define XUYUZU 7

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + XL, lbi0, ubi0, lbn0, hbn0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + YL, lbn1, hbn1, lbi0, ubi0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + ZL, lbn1, hbn1, lbn1, hbn1, lbi0, ubi0, MPI_moment_t);

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + XU, lbi1, ubi1, lbn0, hbn0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + YU, lbn1, hbn1, lbi1, ubi1, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + ZU, lbn1, hbn1, lbn1, hbn1, lbi1, ubi1, MPI_moment_t);

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + XL, lbe0, ube0, lbn0, hbn0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + YL, lbn1, hbn1, lbe0, ube0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + ZL, lbn1, hbn1, lbn1, hbn1, lbe0, ube0, MPI_moment_t);

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + XU, lbe1, ube1, lbn0, hbn0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + YU, lbn1, hbn1, lbe1, ube1, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + ZU, lbn1, hbn1, lbn1, hbn1, lbe1, ube1, MPI_moment_t);

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLYL, lbi0, ubi0, lbi0, ubi0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLYU, lbi0, ubi0, lbi1, ubi1, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUYL, lbi1, ubi1, lbi0, ubi0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUYU, lbi1, ubi1, lbi1, ubi1, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLZL, lbi0, ubi0, lbn0, hbn0, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XLZU, lbi0, ubi0, lbn0, hbn0, lbi1, ubi1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUZL, lbi1, ubi1, lbn0, hbn0, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + XUZU, lbi1, ubi1, lbn0, hbn0, lbi1, ubi1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YLZL, lbn0, hbn0, lbi0, ubi0, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YLZU, lbn0, hbn0, lbi0, ubi0, lbi1, ubi1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YUZL, lbn0, hbn0, lbi1, ubi1, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 6 + YUZU, lbn0, hbn0, lbi1, ubi1, lbi1, ubi1, MPI_moment_t);

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLYL, lbe0, ube0, lbe0, ube0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLYU, lbe0, ube0, lbe1, ube1, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUYL, lbe1, ube1, lbe0, ube0, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUYU, lbe1, ube1, lbe1, ube1, lbn0, hbn0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLZL, lbe0, ube0, lbn0, hbn0, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XLZU, lbe0, ube0, lbn0, hbn0, lbe1, ube1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUZL, lbe1, ube1, lbn0, hbn0, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + XUZU, lbe1, ube1, lbn0, hbn0, lbe1, ube1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YLZL, lbn0, hbn0, lbe0, ube0, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YLZU, lbn0, hbn0, lbe0, ube0, lbe1, ube1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YUZL, lbn0, hbn0, lbe1, ube1, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 6 + YUZU, lbn0, hbn0, lbe1, ube1, lbe1, ube1, MPI_moment_t);

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYLZL, lbi0, ubi0, lbi0, ubi0, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYLZU, lbi0, ubi0, lbi0, ubi0, lbi1, ubi1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYUZL, lbi0, ubi0, lbi1, ubi1, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XLYUZU, lbi0, ubi0, lbi1, ubi1, lbi1, ubi1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYLZL, lbi1, ubi1, lbi0, ubi0, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYLZU, lbi1, ubi1, lbi0, ubi0, lbi1, ubi1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYUZL, lbi1, ubi1, lbi1, ubi1, lbi0, ubi0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_send_bnd_t + 18 + XUYUZU, lbi1, ubi1, lbi1, ubi1, lbi1, ubi1, MPI_moment_t);

		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYLZL, lbe0, ube0, lbe0, ube0, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYLZU, lbe0, ube0, lbe0, ube0, lbe1, ube1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYUZL, lbe0, ube0, lbe1, ube1, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XLYUZU, lbe0, ube0, lbe1, ube1, lbe1, ube1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYLZL, lbe1, ube1, lbe0, ube0, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYLZU, lbe1, ube1, lbe0, ube0, lbe1, ube1, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYUZL, lbe1, ube1, lbe1, ube1, lbe0, ube0, MPI_moment_t);
		Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_recv_bnd_t + 18 + XUYUZU, lbe1, ube1, lbe1, ube1, lbe1, ube1, MPI_moment_t);

		const int ds = INX / 2;
		int ci, xlb, xub, ylb, yub, zlb, zub;

		for (int k = 0; k < 2; k++) {
			zlb = FBW + k * (INX / 2);
			zub = zlb + ds - 1;
			for (int j = 0; j < 2; j++) {
				ylb = FBW + j * (INX / 2);
				yub = ylb + ds - 1;
				for (int i = 0; i < 2; i++) {
					xlb = FBW + i * (INX / 2);
					xub = xlb + ds - 1;
					ci = 4 * k + 2 * j + i;
					Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_comm_child1_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_moment_t);
					Array3d<Taylor, FNX, FNX, FNX>::mpi_datatype(MPI_comm_taylor_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_taylor_t);
#ifdef USE_FMM_ANGULAR
					Array3d<Taylor, FNX, FNX, FNX>::mpi_datatype(MPI_comm_taylor2_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_taylor2_t);
#endif
					Array3d<Moment, FNX, FNX, FNX>::mpi_datatype(MPI_comm_child3_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_4force_t);
				}
			}
		}
		face_id[XL] = 0;
		face_id[XU] = 1;
		face_id[YL] = 2;
		face_id[YU] = 3;
		face_id[ZL] = 4;
		face_id[ZU] = 5;
		face_id[XLYL + 6] = 6;
		face_id[XLYU + 6] = 7;
		face_id[XLZL + 6] = 8;
		face_id[XLZU + 6] = 9;
		face_id[XUYL + 6] = 10;
		face_id[XUYU + 6] = 11;
		face_id[XUZL + 6] = 12;
		face_id[XUZU + 6] = 13;
		face_id[YLZL + 6] = 14;
		face_id[YLZU + 6] = 15;
		face_id[YUZL + 6] = 16;
		face_id[YUZU + 6] = 17;
		face_id[XLYLZL + 18] = 18;
		face_id[XLYLZU + 18] = 19;
		face_id[XLYUZL + 18] = 20;
		face_id[XLYUZU + 18] = 21;
		face_id[XUYLZL + 18] = 22;
		face_id[XUYLZU + 18] = 23;
		face_id[XUYUZL + 18] = 24;
		face_id[XUYUZU + 18] = 25;
		face_opp_id[XU] = 0;
		face_opp_id[XL] = 1;
		face_opp_id[YU] = 2;
		face_opp_id[YL] = 3;
		face_opp_id[ZU] = 4;
		face_opp_id[ZL] = 5;
		face_opp_id[XUYU + 6] = 6;
		face_opp_id[XUYL + 6] = 7;
		face_opp_id[XUZU + 6] = 8;
		face_opp_id[XUZL + 6] = 9;
		face_opp_id[XLYU + 6] = 10;
		face_opp_id[XLYL + 6] = 11;
		face_opp_id[XLZU + 6] = 12;
		face_opp_id[XLZL + 6] = 13;
		face_opp_id[YUZU + 6] = 14;
		face_opp_id[YUZL + 6] = 15;
		face_opp_id[YLZU + 6] = 16;
		face_opp_id[YLZL + 6] = 17;
		face_opp_id[XUYUZU + 18] = 18;
		face_opp_id[XUYUZL + 18] = 19;
		face_opp_id[XUYLZU + 18] = 20;
		face_opp_id[XUYLZL + 18] = 21;
		face_opp_id[XLYUZU + 18] = 22;
		face_opp_id[XLYUZL + 18] = 23;
		face_opp_id[XLYLZU + 18] = 24;
		face_opp_id[XLYLZL + 18] = 25;
		initialized = true;
	}
}

#endif
