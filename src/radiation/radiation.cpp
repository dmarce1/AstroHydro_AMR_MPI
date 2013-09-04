/*
 * radiation.cpp
 *
 *  Created on: May 17, 2013
 *      Author: dmarce1
 */

//#define LINEAR_COUPLING
#include "../tag.h"
#include "../indexer3d.h"
#include "radiation.h"
#include "../physical_constants.h"
#ifdef USE_RADIATION

Real Radiation::dt;
MPI_Datatype Radiation::MPI_send_bnd_t[6];
MPI_Datatype Radiation::MPI_recv_bnd_t[6];

Radiation::ifunc_t Radiation::cs[RAD_INS_CNT + 1] = { &Radiation::compute_gradE, &Radiation::gradE_real_boundary_begin_loop,
		&Radiation::gradE_real_boundary_communicate, &Radiation::gradE_real_boundary_wait, &Radiation::gradE_real_boundary_end_loop,
		&Radiation::gradE_amr_boundary_communicate, &Radiation::gradE_amr_boundary_wait, &Radiation::null };

void Radiation::MPI_datatypes_init() {
	static bool initialized = false;
	if (!initialized) {
		const int lbi = 1;
		const int ubi = PNX - 2;
		const int lbe = 0;
		const int ube = PNX - 1;
		MPI_Datatype MPI_3vec_t;
		MPI_Type_contiguous(sizeof(_3Vec), MPI_BYTE, &MPI_3vec_t);
		MPI_Type_commit(&MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd_t + XL, lbi, lbi, lbi, ubi, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd_t + YL, lbe, ube, lbi, lbi, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd_t + ZL, lbe, ube, lbe, ube, lbi, lbi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd_t + XU, ubi, ubi, lbi, ubi, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd_t + YU, lbe, ube, ubi, ubi, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd_t + ZU, lbe, ube, lbe, ube, ubi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd_t + XL, lbe, lbe, lbi, ubi, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd_t + YL, lbe, ube, lbe, lbe, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd_t + ZL, lbe, ube, lbe, ube, lbe, lbe, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd_t + XU, ube, ube, lbi, ubi, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd_t + YU, lbe, ube, ube, ube, lbi, ubi, MPI_3vec_t);
		Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd_t + ZU, lbe, ube, lbe, ube, ube, ube, MPI_3vec_t);
		initialized = true;
	}
}

Real Radiation::change_max() {
	Real derad, q, tmp;
	q = 0.0;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		Radiation* g = dynamic_cast<Radiation*>(get_local_node(l));
		for (int k = 1; k < PNX - 1; k++) {
			for (int j = 1; j < PNX - 1; j++) {
				for (int i = 1; i < PNX - 1; i++) {
					derad = fabs(log(fabs(g->erad(i, j, k) / g->erad1(i, j, k))));
					q = max(derad, q);
				}
			}
		}
	}
	tmp = q;
	MPI_Allreduce(&tmp, &q, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD );
	return q;
}

void Radiation::compute_gradE() {
	Real h = 0.5 / get_dx();
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				gradE(i, j, k)[0] = (erad(i + 1, j, k) - erad(i - 1, j, k)) * h;
				gradE(i, j, k)[1] = (erad(i, j + 1, k) - erad(i, j - 1, k)) * h;
				gradE(i, j, k)[2] = (erad(i, j, k + 1) - erad(i, j, k - 1)) * h;
			}
		}
	}
	ip++;
}

Real Radiation::chi(Real rho) {
	return sigma_T * rho;
}

Real Radiation::kappa_E(Real rho) {
	return sigma_T * rho;
}

Real Radiation::kappa_p(Real rho) {
	return sigma_T * rho;
}

Real Radiation::flux_limiter(Real R) {
	Real Rinv, l;
	if (R < 1.0e-3) {
		l = 1.0 / 3.0;
	} else {
		Rinv = 1.0 / R;
		if (R > 1.0e+3) {
			l = Rinv * (1.0 - Rinv);
		} else {
			l = Rinv * (1.0 / tanh(R) - Rinv);
		}
	}
	return l;
}

void Radiation::compute_coupling() {
	Real T;
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				T = eint(i, j, k) / cv(i, j, k);
				Bp(i, j, k) = PhysicalConstants::sigma * pow(T, 4) / M_PI;
			}
		}
	}
}

void Radiation::compute_difco() {
	const Real dxinv = 1.0 / get_dx();
	Real face_chi, face_rho, gradE_l, gradEt1, gradEt2, abs_gradE, R, lambda, face_E_inv;
	gradEt1 = gradEt2 = 0.0;
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX; i++) {
				face_rho = 0.5 * (rho(i, j, k) + rho(i - 1, j, k));
				face_chi = chi(face_rho);
				face_E_inv = 2.0 / (erad(i, j, k) + erad(i - 1, j, k));
				gradE_l = (erad(i, j, k) - erad(i - 1, j, k)) * dxinv * face_E_inv;
				gradEt1 = 0.5 * (gradE(i, j, k)[1] + gradE(i - 1, j, k)[1]) * face_E_inv;
				gradEt2 = 0.5 * (gradE(i, j, k)[2] + gradE(i - 1, j, k)[2]) * face_E_inv;
				abs_gradE = sqrt(gradE_l * gradE_l + gradEt1 * gradEt1 + gradEt2 * gradEt2);
				R = abs_gradE / face_chi;
				lambda = flux_limiter(R);
				D1(i, j, k) = PhysicalConstants::c * lambda / face_chi;
			}
		}
	}
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				face_rho = 0.5 * (rho(i, j, k) + rho(i, j - 1, k));
				face_chi = chi(face_rho);
				face_E_inv = 2.0 / (erad(i, j, k) + erad(i, j - 1, k));
				gradE_l = (erad(i, j, k) - erad(i, j - 1, k)) * dxinv * face_E_inv;
				gradEt1 = 0.5 * (gradE(i, j, k)[0] + gradE(i, j - 1, k)[0]) * face_E_inv;
				gradEt2 = 0.5 * (gradE(i, j, k)[2] + gradE(i, j - 1, k)[2]) * face_E_inv;
				abs_gradE = sqrt(gradE_l * gradE_l + gradEt1 * gradEt1 + gradEt2 * gradEt2);
				R = abs_gradE / face_chi;
				lambda = flux_limiter(R);
				D2(i, j, k) = PhysicalConstants::c * lambda / face_chi;
			}
		}
	}
	for (int k = 1; k < PNX; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				face_rho = 0.5 * (rho(i, j, k) + rho(i, j, k - 1));
				face_chi = chi(face_rho);
				face_E_inv = 2.0 / (erad(i, j, k) + erad(i, j, k - 1));
				gradE_l = (erad(i, j, k) - erad(i, j, k - 1)) * dxinv * face_E_inv;
				gradEt1 = 0.5 * (gradE(i, j, k)[0] + gradE(i, j, k - 1)[0]) * face_E_inv;
				gradEt2 = 0.5 * (gradE(i, j, k)[1] + gradE(i, j, k - 1)[1]) * face_E_inv;
				abs_gradE = sqrt(gradE_l * gradE_l + gradEt1 * gradEt1 + gradEt2 * gradEt2);
				R = abs_gradE / face_chi;
				lambda = flux_limiter(R);
				D3(i, j, k) = PhysicalConstants::c * lambda / face_chi;
			}
		}
	}
}

void Radiation::compute_RHS() {
	Real eta, kappa, alpha;
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				kappa = kappa_p(rho(i, j, k));
				eta = 16.0 * M_PI * Bp(i, j, k) * kappa * dt / eint(i, j, k);
#ifdef LINEAR_COUPLING
				alpha = 1.0;
#else
				alpha = (4.0 * eint0(i, j, k) / eint(i, j, k) - 3.0);
#endif
				S(i, j, k) = erad0(i, j, k) + alpha * ((4.0 * M_PI * kappa * Bp(i, j, k)) / (1.0 + eta)) * dt;
			}
		}
	}
}

void Radiation::compute_eint() {
	Real eta, a, b;
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				eta = 16.0 * M_PI * Bp(i, j, k) * kappa_p(rho(i, j, k)) * dt / eint(i, j, k);
				a = (PhysicalConstants::c * kappa_E(rho(i, j, k)) * erad(i, j, k)) * dt;
				b = 12.0 * M_PI * Bp(i, j, k) * kappa_p(rho(i, j, k)) * dt;
				eint(i, j, k) = (eint0(i, j, k) + a + b) / (1.0 + eta);
			}
		}
	}
}

Real* Radiation::Jacobian(int i, int j, int k) const {
	static Real a[7];
	Real* aptr = a + 3;
	const Real h = dt / get_dx() / get_dx();
	Real eta;
	eta = 16.0 * M_PI * Bp(i, j, k) * kappa_p(rho(i, j, k)) * dt / eint(i, j, k);
	aptr[-3] = -h * D3(i, j, k);
	aptr[-2] = -h * D2(i, j, k);
	aptr[-1] = -h * D1(i, j, k);
	aptr[+1] = -h * D1(i + 1, j, k);
	aptr[+2] = -h * D2(i, j + 1, k);
	aptr[+3] = -h * D3(i, j, k + 1);
	Real J0 = 1.0;
	J0 += (D3(i, j, k) + D2(i, j, k) + D1(i, j, k) + D1(i + 1, j, k) + D2(i, j + 1, k) + D3(i, j, k + 1)) * h;
	J0 += PhysicalConstants::c * kappa_E(rho(i, j, k)) * dt / (1.0 + eta);
	aptr[0] = J0;
//	printf("%e %e %e %e %e %e %e\n", a[0], a[1], a[2], a[3], a[4], a[5], a[6]);
	return aptr;
}

void Radiation::vdown_init_compute(int) {
	const Real dx2 = get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				if (!poisson_zone_is_refined(i, j, k)) {
					const Real* J = Jacobian(i, j, k);
					Real d = 0.0;
					d -= J[-3] * phi(i, j, k - 1) / J[0];
					d -= J[-2] * phi(i, j - 1, k) / J[0];
					d -= J[-1] * phi(i - 1, j, k) / J[0];
					d -= phi(i, j, k);
					d -= J[+1] * phi(i + 1, j, k) / J[0];
					d -= J[+2] * phi(i, j + 1, k) / J[0];
					d -= J[+3] * phi(i, j, k + 1) / J[0];
					d += S(i, j, k) / J[0];
					dphi(i, j, k) = d;
				}
			}
		}
	}
	inc_instruction_pointer();
}

void Radiation::vdown_init_adjust_send(int) {
	Vector<int, 3> lb, ub;
	Real v;
	int i, cnt, tag, dir;
	for (int face = 0; face < 6; face++) {
		dir = face / 2;
		if (is_amr_bound(face)) {
			mpi_buffer[face] = new Real[PNX * PNX / 4];
			i = 1 + (PNX - 2) * (face & 1);
			cnt = 0;
			for (int k = 1; k < PNX - 1; k += 2) {
				for (int j = 1; j < PNX - 1; j += 2) {
					if (dir == 0) {
						v = +(phi(i, j + 0, k + 0) - phi(i - 1, j + 0, k + 0)) * D1(i, j + 0, k + 0);
						v += (phi(i, j + 1, k + 0) - phi(i - 1, j + 1, k + 0)) * D1(i, j + 1, k + 0);
						v += (phi(i, j + 0, k + 1) - phi(i - 1, j + 0, k + 1)) * D1(i, j + 0, k + 1);
						v += (phi(i, j + 1, k + 1) - phi(i - 1, j + 1, k + 1)) * D1(i, j + 1, k + 1);
					} else if (dir == 1) {
						v = +(phi(j + 0, i, k + 0) - phi(j + 0, i - 1, k + 0)) * D2(j, i, k);
						v += (phi(j + 1, i, k + 0) - phi(j + 1, i - 1, k + 0)) * D2(j + 1, i, k);
						v += (phi(j + 0, i, k + 1) - phi(j + 0, i - 1, k + 1)) * D2(j, i, k + 1);
						v += (phi(j + 1, i, k + 1) - phi(j + 1, i - 1, k + 1)) * D2(j + 1, i, k + 1);
					} else {
						v = +(phi(j + 0, k + 0, i) - phi(j + 0, k + 0, i - 1)) * D3(j, k, i);
						v += (phi(j + 1, k + 0, i) - phi(j + 1, k + 0, i - 1)) * D3(j + 1, k, i);
						v += (phi(j + 0, k + 1, i) - phi(j + 0, k + 1, i - 1)) * D3(j, k + 1, i);
						v += (phi(j + 1, k + 1, i) - phi(j + 1, k + 1, i - 1)) * D3(j + 1, k + 1, i);
					}
					v *= 0.5;
					mpi_buffer[face][cnt] = v;
					cnt++;
				}
			}assert(cnt == (PNX-2)*(PNX-2)/4);
			tag = tag_gen(TAG_DPHI, get_id(), face);
			MPI_Isend(mpi_buffer[face], cnt, MPI_DOUBLE_PRECISION, guard_proc(face), tag, MPI_COMM_WORLD, send_request + face);
		} else {
			send_request[face] = MPI_REQUEST_NULL;
		}
	}
	inc_instruction_pointer();
}

void Radiation::vdown_init_adjust_recv_wait(int) {
	int flag;
	if (amr_cnt > 0) {
		MPI_Testall(amr_cnt, amr_request, &flag, MPI_STATUS_IGNORE );
		if (flag) {
			int cnt, dir;
			Real v, d;
			int i;
			for (int ci = 0; ci < amr_cnt; ci++) {
				dir = amr_face[ci] / 2;
				cnt = 0;
				i = amr_lb[ci][dir] / 2;
				if (i == PNX / 2 - 1) {
					i++;
				} else if (i == PNX - 2) {
					i = PNX - 1;
				}
				if (dir == 0) {
					for (int k = amr_lb[ci][2] / 2; k <= amr_ub[ci][2] / 2; k++) {
						for (int j = amr_lb[ci][1] / 2; j <= amr_ub[ci][1] / 2; j++) {
							v = mpi_buffer[ci][cnt];
							d = -((phi(i, j, k) - phi(i - 1, j, k)) * D1(i, j, k) - v);
							d *= -dt / get_dx() / get_dx();
							if (!poisson_zone_is_refined(i, j, k) && i < PNX - 1) {
								dphi(i, j, k) += d / Jacobian(i, j, k)[0];
							}
							if (!poisson_zone_is_refined(i - 1, j, k) && i > 1) {
								dphi(i - 1, j, k) -= d / Jacobian(i - 1, j, k)[0];
							}
							++cnt;
						}
					}
				} else if (dir == 1) {
					for (int k = amr_lb[ci][2] / 2; k <= amr_ub[ci][2] / 2; k++) {
						for (int j = amr_lb[ci][0] / 2; j <= amr_ub[ci][0] / 2; j++) {
							v = mpi_buffer[ci][cnt];
							d = -((phi(j, i, k) - phi(j, i - 1, k)) * D2(j, i, k) - v);
							d *= -dt / get_dx() / get_dx();
							if (!poisson_zone_is_refined(j, i, k) && i < PNX - 1) {
								dphi(j, i, k) += d / Jacobian(j, i, k)[0];
							}
							if (!poisson_zone_is_refined(j, i - 1, k) && i > 1) {
								dphi(j, i - 1, k) -= d / Jacobian(j, i - 1, k)[0];
							}
							++cnt;
						}
					}
				} else {
					for (int k = amr_lb[ci][1] / 2; k <= amr_ub[ci][1] / 2; k++) {
						for (int j = amr_lb[ci][0] / 2; j <= amr_ub[ci][0] / 2; j++) {
							v = mpi_buffer[ci][cnt];
							d = -((phi(j, k, i) - phi(j, k, i - 1)) * D3(j, k, i) - v);
							d *= -dt / get_dx() / get_dx();
							if (!poisson_zone_is_refined(j, k, i) && i < PNX - 1) {
								dphi(j, k, i) += d / Jacobian(j, k, i)[0];
							}
							if (!poisson_zone_is_refined(j, k, i - 1) && i > 1) {
								dphi(j, k, i - 1) -= d / Jacobian(j, k, i - 1)[0];
							}
							++cnt;
						}
					}
				}
				delete[] mpi_buffer[ci];
			}
			inc_instruction_pointer();
		}
	} else {
		inc_instruction_pointer();
	}
}

void Radiation::relax_compute(int) {
	Real dt_dx2inv = dt / get_dx() / get_dx();
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = (j + k) % 2 + 1; i < PNX - 1; i += 2) {
				const Real* J = Jacobian(i, j, k);
				Real d = 0.0;
				d -= J[-3] * dphi(i, j, k - 1);
				d -= J[-2] * dphi(i, j - 1, k);
				d -= J[-1] * dphi(i - 1, j, k);
				d -= J[0] * dphi(i, j, k);
				d -= J[+1] * dphi(i + 1, j, k);
				d -= J[+2] * dphi(i, j + 1, k);
				d -= J[+3] * dphi(i, j, k + 1);
				d /= J[0];
				d += dphi1(i, j, k);
				dphi(i, j, k) += (3.0/4.0)*d;
			}
		}
	}
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = (j + k + 1) % 2 + 1; i < PNX - 1; i += 2) {
				const Real* J = Jacobian(i, j, k);
				Real d = 0.0;
				d -= J[-3] * dphi(i, j, k - 1);
				d -= J[-2] * dphi(i, j - 1, k);
				d -= J[-1] * dphi(i - 1, j, k);
				d -= J[0] * dphi(i, j, k);
				d -= J[+1] * dphi(i + 1, j, k);
				d -= J[+2] * dphi(i, j + 1, k);
				d -= J[+3] * dphi(i, j, k + 1);
				d /= J[0];
				d += dphi1(i, j, k);
				dphi(i, j, k) += (3.0/4.0)*d;
			}
		}
	}
	inc_instruction_pointer();
}

void Radiation::forces_compute(int) {
	const Real dxinv = 1.0 / get_dx();
	for (int k = 1; k < PNX; k++) {
		for (int j = 1; j < PNX; j++) {
			for (int i = 1; i < PNX; i++) {
				fx(i, j, k) = -D1(i, j, k) * (phi(i, j, k) - phi(i - 1, j, k)) * dxinv;
				fy(i, j, k) = -D2(i, j, k) * (phi(i, j, k) - phi(i, j - 1, k)) * dxinv;
				fz(i, j, k) = -D3(i, j, k) * (phi(i, j, k) - phi(i, j, k - 1)) * dxinv;
			}
		}
	}
	inc_instruction_pointer();
}

void Radiation::residual_error_compute(int) {
	Real sum;
	sum = 0.0;
	const Real h3 = pow(get_dx(), 3);
	const Real h1 = get_dx();
	int i, j, k;
	Real df, eta;
	const Real h1inv = 1.0 / h1;
	for (k = 1; k < PNX - 1; k++) {
		for (j = 1; j < PNX - 1; j++) {
			for (i = 1; i < PNX - 1; i++) {
				if (!poisson_zone_is_refined(i, j, k)) {
					df = (phi(i, j, k) - S(i, j, k));
					df += (fx(i + 1, j, k) - fx(i, j, k)) * h1inv * dt;
					df += (fy(i, j + 1, k) - fy(i, j, k)) * h1inv * dt;
					df += (fz(i, j, k + 1) - fz(i, j, k)) * h1inv * dt;
					eta = 16.0 * M_PI * Bp(i, j, k) * kappa_p(rho(i, j, k)) * dt / eint(i, j, k);
					df += PhysicalConstants::c * kappa_E(rho(i, j, k)) * erad(i, j, k) * dt / (1.0 + eta);
					dphi(i, j, k) = df;
					sum += fabs(dphi(i, j, k)) * h3;
				} else {
					dphi(i, j, k) = 0.0;
				}
			}
		}
	}
	inc_instruction_pointer();
	st0 = sum;
}

void Radiation::boundary_communicate() {
	MPI_datatypes_init();
	Radiation** procs = new Radiation*[get_local_node_cnt()];
	for (int i = 0; i < get_local_node_cnt(); i++) {
		procs[i] = dynamic_cast<Radiation*>(get_local_node(i));
	}
//run_program(list, get_local_node_cnt(), cs, RAD_INS_CNT);
	int last_ip;
	bool done;
	for (int i = 0; i < get_local_node_cnt(); i++) {
		procs[i]->ip = 0;
	}
	do {
		done = true;
		for (int i = 0; i < get_local_node_cnt(); i++) {
			if (procs[i]->ip < RAD_INS_CNT) {
				do {
					last_ip = procs[i]->ip;
					(procs[i]->*cs[procs[i]->ip])();
				} while (last_ip != procs[i]->ip && procs[i]->ip < RAD_INS_CNT);
				if (procs[i]->ip < RAD_INS_CNT) {
					done = false;
				}

			}
		}
	} while (!done);
	delete[] procs;
}

void Radiation::allocate_arrays() {
	Bp.allocate();
	eint.allocate();
	eint0.allocate();
	erad0.allocate();
	eint1.allocate();
	erad1.allocate();
	cv.allocate();
	rho.allocate();
	D1.allocate();
	D2.allocate();
	D3.allocate();
	gradE.allocate();
	MultiGrid::allocate_arrays();
}

void Radiation::deallocate_arrays() {
	Bp.deallocate();
	eint.deallocate();
	eint0.deallocate();
	erad0.deallocate();
	eint1.deallocate();
	erad1.deallocate();
	cv.deallocate();
	rho.deallocate();
	D1.deallocate();
	D2.deallocate();
	D3.deallocate();
	gradE.deallocate();
	MultiGrid::deallocate_arrays();
}

void Radiation::store0() {
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				erad0(i, j, k) = erad(i, j, k);
				eint0(i, j, k) = eint(i, j, k);
			}
		}
	}
}

void Radiation::store1() {
	for (int k = 1; k < PNX - 1; k++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int i = 1; i < PNX - 1; i++) {
				erad1(i, j, k) = erad(i, j, k);
				eint1(i, j, k) = eint(i, j, k);
			}
		}
	}
}

Real Radiation::implicit_solve(Real _dt) {
	const Real chng_max = 0.1;
	MPI_datatypes_init();
	dt = _dt;
	Real err0, err;

	Radiation::boundary_communicate();

	for (int l = 0; l < get_local_node_cnt(); l++) {
		Radiation* g = dynamic_cast<Radiation*>(get_local_node(l));
		g->store0();
		g->store1();
		g->compute_difco();
		g->compute_coupling();
		g->compute_RHS();
	}
	err0 = vcycle();
	if (MPI_rank() == 0) {
		printf("1.0=%e\n", err0);
	}
	do {
#ifndef LINEAR_COUPLING
		for (int l = 0; l < get_local_node_cnt(); l++) {
			Radiation* g = dynamic_cast<Radiation*>(get_local_node(l));
			g->compute_eint();
		}
		//	floors();
		for (int l = 0; l < get_local_node_cnt(); l++) {
			Radiation* g = dynamic_cast<Radiation*>(get_local_node(l));
			g->compute_coupling();
			g->compute_RHS();
			g->store1();
		}
#endif
		err = vcycle();
		if (MPI_rank() == 0) {
			printf("%e\n", err / err0);
		}
	} while (err / err0 > 1.0e-3);
	for (int l = 0; l < get_local_node_cnt(); l++) {
		Radiation* g = dynamic_cast<Radiation*>(get_local_node(l));
		g->compute_eint();
	}
	return dt * min(1.0, chng_max / change_max());
}

Real Radiation::erad(int i, int j, int k) const {
	return get_phi(i, j, k);
}
Real& Radiation::erad(int i, int j, int k) {
	return *(get_phi_ptr(i, j, k));
}

Radiation::Radiation() {
// TODO Auto-generated constructor stub

}

Radiation::~Radiation() {
// TODO Auto-generated destructor stub
}

void Radiation::gradE_amr_boundary_communicate() {
	if (amr_cnt > 0) {
		int cnt, tag_send;
		for (int i = 0; i < amr_cnt; i++) {
			int sz = (PNX - 2) * (PNX - 2);
			mpi_3vec_buffer[i] = new _3Vec[sz];
			cnt = 0;
			for (int l = amr_lb[i][2]; l <= amr_ub[i][2]; l++) {
				for (int k = amr_lb[i][1]; k <= amr_ub[i][1]; k++) {
					for (int j = amr_lb[i][0]; j <= amr_ub[i][0]; j++) {
						mpi_3vec_buffer[i][cnt] = gradE.interp(j, k, l);
						++cnt;
					}
				}
			}
			tag_send = tag_gen(TAG_PHI, amr_id[i], amr_face[i]);
			MPI_Isend(mpi_3vec_buffer[i], cnt, MPI_DOUBLE_PRECISION, amr_child_proc[i], tag_send, MPI_COMM_WORLD, amr_request + i);
		}
	}
	ip++;
}

void Radiation::gradE_amr_boundary_wait() {
	int flag;
	bool rc;
	MPI_Testall(amr_cnt, amr_request, &flag, MPI_STATUS_IGNORE );
	rc = flag;
	if (rc) {
		for (int i = 0; i < amr_cnt; i++) {
			delete[] mpi_3vec_buffer[i];
		}
	}
	if (rc) {
		ip++;
	}
}

void Radiation::gradE_real_boundary_communicate() {
	const int dir = cx;
	int tag_send, tag_recv;
	Radiation* sibling;
	MPI_Datatype send_type, recv_type;
	for (int face = 2 * dir; face < 2 * dir + 2; face++) {
		sibling = dynamic_cast<Radiation*>(get_sibling(face));
		tag_recv = tag_gen(TAG_PHI, get_id(), face);
		if (is_real_bound(face)) {
			tag_send = tag_gen(TAG_PHI, sibling->get_id(), (face ^ 1));
			send_type = MPI_send_bnd_t[face];
			recv_type = MPI_recv_bnd_t[face];
			MPI_Isend(gradE.ptr(), 1, send_type, sibling->proc(), tag_send, MPI_COMM_WORLD, send_request + face);
			MPI_Irecv(gradE.ptr(), 1, recv_type, guard_proc(face), tag_recv, MPI_COMM_WORLD, recv_request + face);
		} else if (!is_phys_bound(face)) {
			MPI_Irecv(gradE.ptr(), 1, MPI_recv_bnd_t[face], guard_proc(face), tag_recv, MPI_COMM_WORLD, recv_request + face);
			send_request[face] = MPI_REQUEST_NULL;
		} else {
			send_request[face] = MPI_REQUEST_NULL;
			recv_request[face] = MPI_REQUEST_NULL;
		}
	}
	ip++;
}

void Radiation::null() {

}

void Radiation::gradE_real_boundary_wait() {
	const int dir = cx;
	bool rc;
	int flag_recv, flag_send;
	MPI_Testall(2, recv_request + 2 * dir, &flag_recv, MPI_STATUS_IGNORE );
	MPI_Testall(2, send_request + 2 * dir, &flag_send, MPI_STATUS_IGNORE );
	rc = flag_recv && flag_send;
	if (rc) {
		ip++;
	}
}

void Radiation::gradE_real_boundary_begin_loop() {
	cx = 0;
	ip++;
	ax = ip;
}

void Radiation::gradE_real_boundary_end_loop() {
	cx++;
	if (cx < 3) {
		ip = ax;
	} else {
		ip++;
	}
}

#endif
