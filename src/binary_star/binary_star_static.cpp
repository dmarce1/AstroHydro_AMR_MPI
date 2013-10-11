#include "binary_star.h"
#include "../hydro_FMM_grid/hydro_FMM_grid.h"

#ifdef HYDRO_GRAV_GRID
#ifdef BINARY_STAR

#define CHECKPT_FREQ 256

//static Real MA = 1.04;
//static Real MD = 0.20;
static Real MA = 1.04;
static Real MD = 0.20;
Real BinaryStar::refine_floor = 1.0e-4;
bool BinaryStar::scf_code = false;
Real BinaryStar::q_ratio = MD / MA;
Real BinaryStar::polyK;
Real BinaryStar::M1 = 5.0e-6, BinaryStar::M2 = 2.5e-6;
Real BinaryStar::Kd, BinaryStar::Ka;
Real BinaryStar::Ax, BinaryStar::Bx, BinaryStar::Cx, BinaryStar::Aphi, BinaryStar::Bphi, BinaryStar::Cphi, BinaryStar::rho_max_a, BinaryStar::rho_max_d;
_3Vec BinaryStar::a0, BinaryStar::d0, BinaryStar::com_vel_correction = 0.0;
Real BinaryStar::lz_t0 = 0.0;
Real BinaryStar::code_to_cm, BinaryStar::code_to_s, BinaryStar::code_to_K, BinaryStar::code_to_g;

void BinaryStar::assign_fracs(Real Hfrac, Real min_phi, Real max_phi) {
	BinaryStar* g;
	Real d, a, minc, maxc, midc, dm, tmp;
	minc = min_phi;
	maxc = max_phi;
	while (log(fabs(minc / maxc)) > 1.0e-6) {
		a = 0.0;
		d = 0.0;
		midc = 0.5 * (maxc + minc);
		for (int l = 0; l < get_local_node_cnt(); l++) {
			g = dynamic_cast<BinaryStar*>(get_local_node(l));
			Real dv = pow(g->get_dx(), 3);
			for (int k = BW; k < GNX - BW; k++) {
				for (int j = BW; j < GNX - BW; j++) {
					for (int i = BW; i < GNX - BW; i++) {
						if (!g->zone_is_refined(i, j, k)) {
							dm = (*g)(i, j, k).rho() * dv;
							if (g->X(i, j, k)[0] < 0.0 || (*g)(i, j, k).phi_eff() > midc) {
								(*g)(i, j, k).set_frac(0, 0.0);
								(*g)(i, j, k).set_frac(1, (*g)(i, j, k).rho());
								d += dm;
							} else {
								(*g)(i, j, k).set_frac(0, (*g)(i, j, k).rho());
								(*g)(i, j, k).set_frac(1, 0.0);
								a += dm;
							}
						}
					}
				}
			}
		}
		tmp = a;
		MPI_Allreduce(&tmp, &a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
		tmp = d;
		MPI_Allreduce(&tmp, &d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
		if (d / (a + d) > Hfrac) {
			minc = midc;
		} else {
			maxc = midc;
		}
		if (MPI_rank() == 0) {
			printf("%e %e %e %e\n", minc, midc, maxc, d / (a + d));
		}
	}
}

Real BinaryStar::virial_error() {
	BinaryStar *g;
	Real tmp, w, t, tot, wtot;
	tot = 0.0;
	wtot = 0.0;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		Real dv = pow(g->get_dx(), 3);
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!g->zone_is_refined(i, j, k)) {
						tot += (*g)(i, j, k).virial(g->X(i, j, k)) * dv;
						wtot += (*g)(i, j, k).phi(g->X(i, j, k)) * dv;
					}
				}
			}
		}
	}
	tmp = tot;
	MPI_Allreduce(&tmp, &tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = wtot;
	MPI_Allreduce(&tmp, &wtot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	return tot / -wtot;
}

void BinaryStar::find_rho_max(Real* rho1, Real* rho2) {
	BinaryStar *g;
	Real tmp;
	*rho1 = 0.0;
	*rho2 = 0.0;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!g->zone_is_refined(i, j, k)) {
						*rho1 = max(*rho1, (*g)(i, j, k).frac(0));
						*rho2 = max(*rho2, (*g)(i, j, k).frac(1));
					}
				}
			}
		}
	}
	tmp = *rho1;
	MPI_Allreduce(&tmp, rho1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD );
	tmp = *rho2;
	MPI_Allreduce(&tmp, rho2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD );
}

#ifdef ZTWD
Real BinaryStar::find_K(int frac, Real phi0, Real center, Real l1_x, Real* span) {
	BinaryStar* g;
	const Real n = 1.5;
	Real K, M, tmp, dv, xmin, xmax, A, df, M0, f;
	_3Vec x0, dx, x;
	bool test;
	x0 = 0.0;
	x0[0] = center;
	if (frac == 0) {
		M0 = M1;
		A = Ka;
	} else {
		M0 = M2;
		A = Kd;
	}
	do {
		xmin = 1.0e+99;
		xmax = -1.0e-99;
		Real dphi;
		M = 0.0;
		df = 0.0;
		for (int l = 0; l < get_local_node_cnt(); l++) {
			g = dynamic_cast<BinaryStar*>(get_local_node(l));
			dv = pow(g->get_dx(), 3);
			for (int k = BW; k < GNX - BW; k++) {
				for (int j = BW; j < GNX - BW; j++) {
					for (int i = BW; i < GNX - BW; i++) {
						if (!g->zone_is_refined(i, j, k)) {
							x = g->X(i, j, k);
							dx = x - x0;
							test = ((frac == 0 && x[0] > l1_x) || (frac == 1 && x[0] < l1_x));
							if (test && (g->geff(i, j, k)).dot(dx) < 0.0) {
								dphi = phi0 - (*g)(i, j, k).phi_eff();
								if (dphi > 0.0) {
									const Real y0 = dphi / (8.0 * A);
									const Real rho0 = pow(y0 * y0 + 2.0 * y0, 3.0 / 2.0);
									df += -3.0 * y0 * (y0 + 1.0) * pow(rho0, 1.0 / 3.0) * dv / A;
									M += rho0 * dv;
									xmax = max(xmax, x[0]);
									xmin = min(xmin, x[0]);
								}
							}
						}
					}
				}
			}
		}
		tmp = df;
		MPI_Allreduce(&tmp, &df, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
		tmp = M;
		MPI_Allreduce(&tmp, &M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
		tmp = xmax;
		MPI_Allreduce(&tmp, &xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD );
		tmp = xmin;
		MPI_Allreduce(&tmp, &xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD );
		*span = xmax - xmin;
		f = M - M0;
		A -= f / df;
	} while (fabs(f / df / A) > 1.0e-6);
	return A;
}

#else

Real BinaryStar::find_K(int frac, Real phi0, Real center, Real l1_x, Real* span) {
	BinaryStar* g;
	const Real n = 1.5;
	Real K, M, tmp, dv, xmin, xmax;
	_3Vec x0, dx, x;
	bool test;
	x0 = 0.0;
	x0[0] = center;
	if (frac == 0) {
		M = M1;
	} else {
		M = M2;
	}
	xmin = 1.0e+99;
	xmax = -1.0e-99;
	Real dphi;
	K = 0.0;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		dv = pow(g->get_dx(), 3);
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!g->zone_is_refined(i, j, k)) {
						x = g->X(i, j, k);
						dx = x - x0;
						test = ((frac == 0 && x[0] > l1_x) || (frac == 1 && x[0] < l1_x));
						if (test && (g->geff(i, j, k)).dot(dx) < 0.0) {
							dphi = phi0 - (*g)(i, j, k).phi_eff();
							if (dphi > 0.0) {
								K += pow(dphi / ((n + 1.0)), n) * dv;
								xmax = max(xmax, x[0]);
								xmin = min(xmin, x[0]);
							}
						}
					}
				}
			}
		}
	}
	tmp = K;
	MPI_Allreduce(&tmp, &K, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = xmax;
	MPI_Allreduce(&tmp, &xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD );
	tmp = xmin;
	MPI_Allreduce(&tmp, &xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD );
	K = pow(K / M, 1.0 / n);
	*span = xmax - xmin;
	return K;
}
#endif

void BinaryStar::find_mass(int frac, Real* mass, Real* com_ptr) {
	BinaryStar* g;
	Real M, tmp, com, dm, dv;
	M = com = 0.0;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		dv = pow(g->get_dx(), 3);
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!g->zone_is_refined(i, j, k)) {
						dm = dv * (*g)(i, j, k).frac(frac);
						M += dm;
						com += dm * g->HydroGrid::xc(i);
					}
				}
			}
		}
	}
	tmp = com;
	MPI_Allreduce(&tmp, &com, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = M;
	MPI_Allreduce(&tmp, &M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	com /= M;
	*com_ptr = com;
	*mass = M;
}

void BinaryStar::next_rho(Real Ka, Real phi0a, Real xa, Real Kd, Real phi0d, Real xd, Real l1_x) {
//	printf( "Next rho...\n");
	_3Vec x, dx, x0a, x0d;
	Real rho1, rho2;
	const Real w = 0.25;
	const Real n = 1.5;
	BinaryStar* g;
	x0a = 0.0;
	x0d = 0.0;
	x0a[0] = xa;
	x0d[0] = xd;
	Real e, K, d, lz, y0, dphi, h;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					x = g->X(i, j, k);
					dx = x - x0d;
					rho1 = rho2 = State::rho_floor / 2.0;
					if (x[0] < l1_x && (g->geff(i, j, k)).dot(dx) < 0.0) {
#ifndef ZTWD
						rho2 = pow(max(phi0d - (*g)(i, j, k).phi_eff(), 0.0) / ((n + 1.0) * Kd), n);
#else
						h = max(phi0d - (*g)(i, j, k).phi_eff(), 0.0);
						y0 = h / (8.0 * Kd);
						rho2 = pow(y0 * y0 + 2.0 * y0, 1.5);
#endif
//						h = max(phi0d - (*g)(i, j, k).phi_eff(), 0.0);
//						y0 = h / (5.0 * Kd);
//						y0 = pow(/*y0 * y0*/ + 2.0 * y0, 1.5);
//						if (rho2 != 0.0) {
//							printf("%e\n", rho2 / y0);
//						}
						rho1 = State::rho_floor;
						K = Kd;
					}
					dx = x - x0a;
					if (x[0] > l1_x && (g->geff(i, j, k)).dot(dx) < 0.0) {
#ifndef ZTWD
						rho1 = pow(max(phi0a - (*g)(i, j, k).phi_eff(), 0.0) / ((n + 1.0) * Ka), n);
#else
						h = max(phi0a - (*g)(i, j, k).phi_eff(), 0.0);
						y0 = h / (8.0 * Ka);
						rho1 = pow(y0 * y0 + 2.0 * y0, 1.5);
#endif
						rho2 = State::rho_floor;
						K = Ka;
					}
					d = rho1 + rho2;
					e = K * pow(d, 1.0 + 1.0 / n) / (State::gamma - 1.0);
					(*g)(i, j, k).set_rho(d * w + (*g)(i, j, k).rho() * (1.0 - w));
					(*g)(i, j, k).set_frac(0, rho1 * w + (*g)(i, j, k).frac(0) * (1.0 - w));
					(*g)(i, j, k).set_frac(1, rho2 * w + (*g)(i, j, k).frac(1) * (1.0 - w));
					d = (*g)(i, j, k).rho();
					e = K * pow(d, 1.0 + 1.0 / n) / (State::gamma - 1.0);
					lz = (x[0] * x[0] + x[1] * x[1]) * d * State::get_omega();
					(*g)(i, j, k).set_et(e);
					if (State::cylindrical) {
						(*g)(i, j, k).set_sx(0.0);
						(*g)(i, j, k).set_sy(lz);
					} else {
						(*g)(i, j, k).set_sx(-d * x[1] * State::get_omega());
						(*g)(i, j, k).set_sy(+d * x[0] * State::get_omega());
					}
					(*g)(i, j, k).set_tau(pow(e, 1.0 / State::gamma));
					(*g)(i, j, k).set_sz(0.0);

				}
			}
		}
	}

}

void BinaryStar::scf_run(int argc, char* argv[]) {
	int maxlev = get_max_level_allowed();
	set_max_level_allowed(min(maxlev, 7));
	int cur_lev = min(maxlev, 7);

	shadow_off();
	setup_grid_structure();

	Real phi_min_a, phi_min_d, l1_phi, l1_x, xd, xa, phi0d, phi0a, l1_x0, com_a, com_d, m_a, m_d, xd0, l2_x, R2, phil1, phil2, l2_phi, Omega, Ra, Rd, Rd0, Ax, Bx,
			Aphi, Bphi, com;
	Real d2 = 0.05;
	Real ff = 2.0;
	Real verr;
	_3Vec O;
	int scf_iter;
	for (scf_iter = 0; scf_iter < 1000; scf_iter++) {
		if (scf_iter % 10 == 0) {
			get_root()->output("S", 2 * scf_iter, GNX, BW);
			if (check_for_refine())
				next_rho(Ka, phi0a, xa, Kd, phi0d, xd, l1_x);
			get_root()->output("S", 2 * scf_iter + 1, GNX, BW);
		}
		if (scf_iter % 20 == 0 && scf_iter != 0) {
			if (cur_lev < maxlev) {
				cur_lev++;
			}
			set_max_level_allowed(cur_lev);
		}
		//	set_poisson_tolerance(1.0e-10);
#ifdef USE_FMM
		FMM_solve();
		max_dt_driver();
#else
		solve_poisson();
#endif
		find_phimins(&phi_min_a, &xa, &phi_min_d, &xd);
		find_mass(0, &m_a, &com_a);
		find_mass(1, &m_d, &com_d);
		find_l(com_a, com_d, &l1_phi, &l1_x, 1);
		phi0d = l1_phi * 0.99 + 0.01 * phi_min_d;
		if (scf_iter == 0) {
			find_l(com_a, com_d, &l2_phi, &l2_x, 2);
			l1_x0 = l1_x;
			xd0 = com_d;
		}
		//	phi0a = phi_min_a * 0.5 + 0.5 * l1_phi;
		phi0a = ff * l1_phi;
		com = (com_a * m_a + com_d * m_d) / (m_a + m_d);
		O = get_origin();
		O[0] += 1.5 * com;
		set_origin(O);
		Ka = find_K(0, phi0a, xa, l1_x, &Ra);
		Kd = find_K(1, phi0d, xd, l1_x, &Rd);
#ifdef ZTWD
		const Real A = 6.0023e+22;
		const Real B = 2.0 * 9.7393e+5;
		const Real G = 6.6745e-8;
		Real g, cm, s;
		g = pow(Kd, -1.5);
		cm = pow(Kd, -0.5);
		s = 1.0;
		g *= pow(A / G, 1.5) / B / B;
		s *= 1.0 / sqrt(B * G);
		cm *= 1.0 / B * sqrt(A / G);
#endif
		ff *= pow(Ka / Kd, 1.0 / 20.0);
		ff = max(1.0, ff);
		//	ff = 1.5;
		//	phi0a = min(phi0d, l1_phi);
		Ra /= 2.0;
		Rd /= 2.0;
		if (scf_iter == 0) {
			Rd0 = Rd;
		}
		Ax = l1_x - 2.0 * Rd0;
		Bx = l1_x;
		Aphi = get_phi_at(Ax, 0.0, 0.0);
		Bphi = get_phi_at(Bx, 0.0, 0.0);
		Omega = sqrt((Aphi - Bphi) / (0.5 * (Ax * Ax - Bx * Bx)));
		State::set_omega(Omega * 0.5 + State::get_omega() * 0.5);
		next_rho(Ka, phi0a, xa, Kd, phi0d, xd, l1_x);
		verr = virial_error();
#ifdef ZTWD
		m_a *= g;
		m_d *= g;
		m_a /= 1.9891e+33;
		m_d /= 1.9891e+33;
		M1 *= pow(MA / m_a, 1.0 / 10.0);
		M2 *= pow(MD / m_d, 1.0 / 10.0);
#endif
		if (MPI_rank() == 0) {
			if (scf_iter % 25 == 0) {
				printf("\n   %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s \n", "g", "cm", "s", "Ma", "Ka", "phi0a",
						"xa", "Md", "Kd", "phi0d", "xd", "Omega", "l1_x", "Rd", "Ra", "com", "virial", "ff", "err");
			}
			printf("%i %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e \n", scf_iter, g, cm, s, m_a, Ka, phi0a, com_a,
					m_d, Kd, phi0d, com_d, Omega, l1_x, Rd, Ra, com, verr, ff, fabs(log(Ka / Kd)));
		}
		if (fabs(log(Ka / Kd)) < 1.0e-4 && scf_iter >= 10) {
			break;
		}
	}
	assign_fracs(0.708333, phi_min_a, phi0a);
	const Real n = 1.5;
	BinaryStar* gg;
	PhysicalConstants::A = Kd;
	const Real A = 6.0023e+22;
	const Real B = 2.0 * 9.7393e+5;
	const Real G = 6.6745e-8;
	Real g = pow(Kd, -1.5);
	Real cm = pow(Kd, -0.5);
	Real s = 1.0;
	g *= pow(A / G, 1.5) / B / B;
	s *= 1.0 / sqrt(B * G);
	cm *= 1.0 / B * sqrt(A / G);
	Real cm3 = cm * cm * cm;
	Real erg = g * cm * cm / s / s;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		gg = dynamic_cast<BinaryStar*>(get_local_node(l));
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					Real e, sr, lz, sz, K, d, O, tau, sx, sy, f1, f2;
					O = State::get_omega();
					_3Vec x = gg->X(i, j, k);
					d = (*gg)(i, j, k).rho();
					lz = (x[0] * x[0] + x[1] * x[1]) * d * O;
					sr = 0.0;
					tau = pow(State::ei_floor, 1.0 / State::gamma);
					sx = -x[1] * d * O;
					sy = +x[0] * d * O;
					sz = 0.0;
					f1 = (*gg)(i, j, k).frac(0);
					f2 = (*gg)(i, j, k).frac(1);
					e = (*gg)(i, j, k).ed() + State::ei_floor;
					d *= g / cm3;
					f1 *= g / cm3;
					f2 *= g / cm3;
					sx *= (g * cm / s) / cm3;
					sy *= (g * cm / s) / cm3;
					sr *= (g * cm / s) / cm3;
					sz *= (g * cm / s) / cm3;
					lz *= (g * cm * cm / s) / cm3;
					e *= erg / cm3;
					(*gg)(i, j, k).set_rho(d);
					(*gg)(i, j, k).set_frac(0, f1);
					(*gg)(i, j, k).set_frac(1, f2);
					(*gg)(i, j, k).set_tau(tau);
					if (State::cylindrical) {
						(*gg)(i, j, k).set_sx(sr);
						(*gg)(i, j, k).set_sy(lz);
					} else {
						(*gg)(i, j, k).set_sx(sx);
						(*gg)(i, j, k).set_sy(sy);
					}
					(*gg)(i, j, k).set_sz(sz);
					(*gg)(i, j, k).set_et(e);
				}
			}
		}
	}
	State::ei_floor *= erg / cm3;
	PhysicalConstants::set_cgs();
	State::rho_floor = 0.0;
	Real tmp = refine_floor;
	refine_floor = 0.0;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		gg = dynamic_cast<BinaryStar*>(get_local_node(l));
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					State::rho_floor = max(State::rho_floor, 1.0e-14 * (*gg)(i, j, k).rho());
					refine_floor = max(refine_floor, tmp * (*gg)(i, j, k).rho());
				}
			}
		}
	}
	tmp = State::rho_floor;
	MPI_Allreduce(&tmp, &State::rho_floor, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD );
	tmp = refine_floor;
	MPI_Allreduce(&tmp, &refine_floor, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD );
	State::set_omega(State::get_omega() / s);
	dynamic_cast<HydroGrid*>(get_root())->HydroGrid::mult_dx(cm);
#ifndef USE_FMM
	dynamic_cast<MultiGrid*>(get_root())->MultiGrid::mult_dx(cm);
#endif
	for (int l = 0; l < get_local_node_cnt(); l++) {
		gg = dynamic_cast<BinaryStar*>(get_local_node(l));
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					(*gg)(i, j, k).floor(gg->HydroGrid::X(i, j, k));
				}
			}
		}
	}
	get_root()->output("S", scf_iter, GNX, BW);
}

void BinaryStar::read_from_file(const char* str, int* i1, int* i2) {
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
	fread(&refine_floor, sizeof(Real), 1, fp);
	fread(&State::ei_floor, sizeof(Real), 1, fp);
	fread(&State::rho_floor, sizeof(Real), 1, fp);
	fread(&com_vel_correction, sizeof(_3Vec), 1, fp);
	fread(&O, sizeof(_3Vec), 1, fp);
	fread(&omega, sizeof(Real), 1, fp);
	set_origin(O);
	fread(&last_dt, sizeof(Real), 1, fp);
	fread(i1, sizeof(int), 1, fp);
	fread(i2, sizeof(int), 1, fp);
	get_root()->read_checkpoint(fp);
	fclose(fp);
	State::set_omega(omega);
#ifndef USE_FMM
	int count = OctNode::get_local_node_cnt();
	for (int i = 0; i < count; i++) {
		dynamic_cast<MultiGrid*>(OctNode::get_local_node(i))->phi_calc_amr_bounds();
	}
#endif
#ifdef USE_FMM
	FMM_solve();
#else
	solve_poisson();
#endif
}

void BinaryStar::write_to_file(int i1, int i2, const char* idname) {
	char* fname;
	Real omega = State::get_omega();
	FILE* fp;
	_3Vec O = get_origin();

//	asprintf(&fname, "mv -f checkpoint.hello.%i.bin checkpoint.goodbye.%i.bin 2> /dev/null\n", MPI_rank(), MPI_rank());
//	system(fname);
//	free(fname);
	char dummy;
	if (MPI_rank() != 0) {
		MPI_Recv(&dummy, 1, MPI_BYTE, MPI_rank() - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}
	asprintf(&fname, "checkpoint.%s.%i.bin", idname, MPI_rank());
	fp = fopen(fname, "wb");
	free(fname);
	fwrite(&refine_floor, sizeof(Real), 1, fp);
	fwrite(&State::ei_floor, sizeof(Real), 1, fp);
	fwrite(&State::rho_floor, sizeof(Real), 1, fp);
	fwrite(&com_vel_correction, sizeof(_3Vec), 1, fp);
	fwrite(&O, sizeof(_3Vec), 1, fp);
	fwrite(&omega, sizeof(Real), 1, fp);
	fwrite(&last_dt, sizeof(Real), 1, fp);
	fwrite(&i1, sizeof(int), 1, fp);
	fwrite(&i2, sizeof(int), 1, fp);
	get_root()->write_checkpoint(fp);
	fclose(fp);
	if (MPI_rank() < MPI_size() - 1) {
		MPI_Send(&dummy, 1, MPI_BYTE, MPI_rank() + 1, 0, MPI_COMM_WORLD );
	}
}

void BinaryStar::run(int argc, char* argv[]) {
//	State::turn_off_total_energy();

	if (scf_code) {
		scf_run(argc, argv);
		return;
	}
#ifndef USE_FMM
	set_poisson_tolerance(1.0e-8);
#endif
	shadow_off();

	Real dt;
	bool do_output, last_step;
	int step_cnt = 0;
	int ostep_cnt = 0;
	int nnodes;

	if (argc == 2) {
		scf_code = true;
		scf_run(argc, argv);
		scf_code = false;
//		setup_grid_structure();
	} else {
		read_from_file(argv[2], &step_cnt, &ostep_cnt);
	}

	PhysicalConstants::set_cgs();

#ifndef USE_FMM
	set_poisson_tolerance(1.0e-10);
#endif
	if (MPI_rank() == 0) {
		printf("Beginning evolution with period = %e, omega = %e\n", 2.0 * M_PI / State::get_omega(), State::get_omega());
	}
#ifndef USE_FMM
	hydro_time = poisson_boundary_time = poisson_interior_time = 0.0;
#endif
	State sum;
	if (get_time() == 0.0) {
		get_root()->output("X", 0.0, GNX, BW);
		integrate_conserved_variables(&sum);
		if (State::cylindrical) {
			lz_t0 = sum[State::sy_index];
		}
	}
	int iter = 0;
	do {
		iter++;
		if (iter == 10) {
			//		break;
		}
		integrate_conserved_variables(&sum);
		if (MPI_rank() == 0) {
			FILE* fp = fopen("sums.dat", "at");
			fprintf(fp, "%e ", get_time());
			for (int l = 0; l < STATE_NF; l++) {
				fprintf(fp, "%.10e ", sum[l]);
			}
			fprintf(fp, "\n");
			fclose(fp);
		}
#ifdef USE_FMM
		momentum_sum();
		com_sum();
#else
		if (MPI_rank() == 0) {
			_3Vec com = get_center_of_mass();
			FILE* fp = fopen("com.dat", "at");
			fprintf(fp, "%e %e %e %e %e %e %e\n", get_time(), com[0], com[1], com[2], com_vel_correction[0], com_vel_correction[1], com_vel_correction[2]);
			fclose(fp);
		}
#endif
		//	State::set_drift_tor(((lz_t0 - sum[State::sy_index]) / sum[State::d_index])/dt);
		Real ofreq = (OUTPUT_TIME_FREQ * (2.0 * M_PI) / State::get_omega());
		if (step_cnt % CHECKPT_FREQ == 0) {
			if (step_cnt % (CHECKPT_FREQ * 2) == 0) {
				write_to_file(step_cnt, ostep_cnt, "hello");
			} else {
				write_to_file(step_cnt, ostep_cnt, "goodbye");
			}
		}
		dt = next_dt(&do_output, &last_step, &ostep_cnt, ofreq);
		step(dt);
#ifndef USE_FMM
		_3Vec com = get_center_of_mass();
		Real com_omega = State::get_omega() * 100.0;
		com_vel_correction -= (com_vel_correction * 2.0 + com * com_omega) * dt * com_omega;
		State::set_drift_vel(-com_vel_correction);
#else
		State::set_drift_vel(_3Vec(0.0));
#endif
		diagnostics(dt);
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
				printf("*");
			}
			get_root()->output("X", nint(HydroGrid::get_time() / ofreq), GNX, BW);
		}
		if (MPI_rank() == 0) {
			printf("\n");
		}
	} while (!last_step);
	/*if (MPI_rank() == 0) {
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
	 }*/
}

void BinaryStar::integrate_conserved_variables(Vector<Real, STATE_NF>* sum_ptr) {
	int i;
	BinaryStar* g;
	Vector<Real, STATE_NF> my_sum = 0.0;
	Real dv;
	for (int i = 0; i < get_local_node_cnt(); i++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(i));
		dv = pow(g->get_dx(), 3);
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!g->zone_is_refined(i, j, k)) {
						for (int l = 0; l < STATE_NF; l++) {
							if (l == State::et_index) {
								my_sum[l] += (*g)(i, j, k).conserved_energy(g->HydroGrid::X(i, j, k)) * dv;
							} else {
								my_sum[l] += (*g)(i, j, k)[l] * dv;
							}
						}
					}
				}
			}
		}
	}
	MPI_Allreduce(&my_sum, sum_ptr, STATE_NF, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
}

void BinaryStar::diagnostics(Real dt) {
	int i;
	BinaryStar* b;
	Real dv;
	State u, dudt;
	_3Vec x;

	Real tmp;
	Real lz = 0.0;
	Real dlz_grav = 0.0;
	Real dlz_hydro = 0.0;
	Real dlz_source = 0.0;
	Real dlz_flow_off = 0.0;
	Real da;
	Real etall = 0.0;
	for (int i = 0; i < get_local_node_cnt(); i++) {
		b = dynamic_cast<BinaryStar*>(get_local_node(i));
		dv = pow(b->get_dx(), 3);
		da = pow(b->get_dx(), 2);
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!b->zone_is_refined(i, j, k)) {
						u = (*b)(i, j, k);
						dudt = b->get_dudt(i, j, k);
						x = b->X(i, j, k);
						etall += u.conserved_energy(x) * dv;
						if (State::cylindrical) {
							lz += u.sy() * dv;
//							dlz_hydro += dudt[State::sy_index] * dv;
							const int l = State::sy_index;
							dlz_hydro += (b->get_flux(0, i + 1, j, k)[l] - b->get_flux(0, i, j, k)[l]) * da;
							dlz_hydro += (b->get_flux(1, i, j + 1, k)[l] - b->get_flux(1, i, j, k)[l]) * da;
							dlz_hydro += (b->get_flux(2, i, j, k + 1)[l] - b->get_flux(2, i, j, k)[l]) * da;
						} else {
							lz += (x[0] * u.sy() - x[1] * u.sx()) * dv;
							int l = State::sy_index;
							dlz_hydro += (b->get_flux(0, i + 1, j, k)[l] - b->get_flux(0, i, j, k)[l]) * da * b->HydroGrid::xc(i);
							dlz_hydro += (b->get_flux(1, i, j + 1, k)[l] - b->get_flux(1, i, j, k)[l]) * da * b->HydroGrid::xc(i);
							dlz_hydro += (b->get_flux(2, i, j, k + 1)[l] - b->get_flux(2, i, j, k)[l]) * da * b->HydroGrid::xc(i);
							dlz_hydro -= (*b)(i, j, k).source(b->HydroGrid::X(i, j, k), get_time())[l] * b->HydroGrid::xc(i) * dv;
							l = State::sx_index;
							dlz_hydro -= (b->get_flux(0, i + 1, j, k)[l] - b->get_flux(0, i, j, k)[l]) * da * b->HydroGrid::yc(j);
							dlz_hydro -= (b->get_flux(1, i, j + 1, k)[l] - b->get_flux(1, i, j, k)[l]) * da * b->HydroGrid::yc(j);
							dlz_hydro -= (b->get_flux(2, i, j, k + 1)[l] - b->get_flux(2, i, j, k)[l]) * da * b->HydroGrid::yc(j);
							dlz_hydro += (*b)(i, j, k).source(b->HydroGrid::X(i, j, k), get_time())[l] * b->HydroGrid::yc(j) * dv;
						}
#ifdef USE_FMM
						dlz_grav += b->dlz(i, j, k) * u.rho() * dv;
#else
						dlz_grav += b->glz(i, j, k) * u.rho() * dv;
#endif
						if (State::cylindrical) {
							if (b->is_phys_bound(XL) && i == BW) {
								dlz_flow_off -= (b->get_flux(0, i, j, k))[State::sy_index] * da;
							}
							if (b->is_phys_bound(XU) && i == GNX - BW - 1) {
								dlz_flow_off += (b->get_flux(0, i + 1, j, k))[State::sy_index] * da;
							}
							if (b->is_phys_bound(YL) && j == BW) {
								dlz_flow_off -= (b->get_flux(1, i, j, k))[State::sy_index] * da;
							}
							if (b->is_phys_bound(YU) && j == GNX - BW - 1) {
								dlz_flow_off += (b->get_flux(1, i, j + 1, k))[State::sy_index] * da;
							}
							if (b->is_phys_bound(ZL) && k == BW) {
								dlz_flow_off -= (b->get_flux(2, i, j, k))[State::sy_index] * da;
							}
							if (b->is_phys_bound(ZU) && k == GNX - BW - 1) {
								dlz_flow_off += (b->get_flux(2, i, j, k + 1))[State::sy_index] * da;
							}
						} else {
							if (b->is_phys_bound(XL) && i == BW) {
								dlz_flow_off -= (b->get_flux(0, i, j, k))[State::sy_index] * da * b->HydroGrid::xf(i);
								dlz_flow_off += (b->get_flux(0, i, j, k))[State::sx_index] * da * b->HydroGrid::yc(j);
							}
							if (b->is_phys_bound(XU) && i == GNX - BW - 1) {
								dlz_flow_off += (b->get_flux(0, i + 1, j, k))[State::sy_index] * da * b->HydroGrid::xf(i + 1);
								dlz_flow_off -= (b->get_flux(0, i + 1, j, k))[State::sx_index] * da * b->HydroGrid::yc(j);
							}
							if (b->is_phys_bound(YL) && j == BW) {
								dlz_flow_off -= (b->get_flux(1, i, j, k))[State::sy_index] * da * b->HydroGrid::xc(i);
								dlz_flow_off += (b->get_flux(1, i, j, k))[State::sx_index] * da * b->HydroGrid::yf(j);
							}
							if (b->is_phys_bound(YU) && j == GNX - BW - 1) {
								dlz_flow_off += (b->get_flux(1, i, j + 1, k))[State::sy_index] * da * b->HydroGrid::xc(i);
								dlz_flow_off -= (b->get_flux(1, i, j + 1, k))[State::sx_index] * da * b->HydroGrid::yf(j + 1);
							}
							if (b->is_phys_bound(ZL) && k == BW) {
								dlz_flow_off -= (b->get_flux(2, i, j, k))[State::sy_index] * da * b->HydroGrid::xc(i);
								dlz_flow_off += (b->get_flux(2, i, j, k))[State::sx_index] * da * b->HydroGrid::yc(j);
							}
							if (b->is_phys_bound(ZU) && k == GNX - BW - 1) {
								dlz_flow_off += (b->get_flux(2, i, j, k + 1))[State::sy_index] * da * b->HydroGrid::xc(i);
								dlz_flow_off -= (b->get_flux(2, i, j, k + 1))[State::sx_index] * da * b->HydroGrid::yc(j);
							}
						}
					}
				}
			}
		}
	}
	tmp = etall;
	MPI_Allreduce(&tmp, &etall, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = lz;
	MPI_Allreduce(&tmp, &lz, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = dlz_flow_off;
	MPI_Allreduce(&tmp, &dlz_flow_off, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = dlz_hydro;
	MPI_Allreduce(&tmp, &dlz_hydro, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = dlz_grav;
	MPI_Allreduce(&tmp, &dlz_grav, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	dlz_hydro -= dlz_flow_off;
//	dlz_hydro -= dlz_grav;
	if (MPI_rank() == get_root()->proc()) {
		etall += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::et_index];
		etall += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::pot_index];
		//	lz += dynamic_cast<const BinaryStar*>(get_root())->get_flow_off()[State::sy_index];
		FILE* fp = fopen("diag.dat", "at");
		fprintf(fp, "%.12e %.18e %.18e  %.12e\n", get_time(), lz + get_flow_off()[State::sy_index], etall, get_flow_off()[State::d_index]);
		fclose(fp);
	}
}

void BinaryStar::find_phimins(Real* phi_min_a, Real* a_x, Real* phi_min_d, Real* d_x) {
	BinaryStar* g;
	Vector<int, 3> loc;
	Real pma = +1.0e+99;
	Real pmd = +1.0e+99;
	int i, j, k;
	int plane_index;
	Real xa, xd;
	Real r[2], r0[2];
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		loc = g->get_location();
		plane_index = (1 << g->get_level()) / 2;
		if (loc[1] == plane_index && loc[2] == plane_index) {
			j = k = BW;
			for (i = BW; i < GNX - BW; i++) {
				if (!g->zone_is_refined(i, j, k) && (*g)(i, j, k).rho() > 10.0 * State::rho_floor) {
					if ((*g)(i, j, k).frac(0) < (*g)(i, j, k).frac(1)) {
						if ((*g)(i, j, k).phi_eff() < pmd) {
							pmd = (*g)(i, j, k).phi_eff();
							xd = g->HydroGrid::xc(i);
						}
					} else {
						if ((*g)(i, j, k).phi_eff() < pma) {
							pma = (*g)(i, j, k).phi_eff();
							xa = g->HydroGrid::xc(i);
						}
					}
				}
			}
		}
	}
	r[0] = pma;
	r[1] = xa;
	MPI_Allreduce(r, r0, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD );
	*phi_min_a = r0[0];
	*a_x = r0[1];
	r[0] = pmd;
	r[1] = xd;
	MPI_Allreduce(r, r0, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD );
	*phi_min_d = r0[0];
	*d_x = r0[1];
	//	printf("%e %e \n", xd, pmd);
	//	printf("phi_min_a = %e\n", phi_min_a);
}

Real BinaryStar::find_acc_com() {
	BinaryStar* g;
	Real com = 0.0;
	Real tmp, dv, mtot = 0.0, dm;
	int i, j, k;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		dv = pow(g->get_dx(), 3);
		for (k = BW; k < GNX - BW; k++) {
			for (j = BW; j < GNX - BW; j++) {
				for (i = BW; i < GNX - BW; i++) {
					dm = (*g)(i, j, k).frac(0) * dv;
					com += g->HydroGrid::xc(i) * dm;
					mtot += dm;
				}
			}
		}
	}
	tmp = com;
	MPI_Allreduce(&tmp, &com, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = mtot;
	MPI_Allreduce(&tmp, &mtot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	return com / mtot;
}

void BinaryStar::next_omega(Real Kd) {
	BinaryStar* g;
	int i, j, k;
	Real frot, ftot, x, hp, hm, phip, phim, tmp;
	Real n = 1.5, d, dp, dm;
#ifdef USE_FMM
	const int o = 0;
#else
	const int o = BW - 1;
#endif
	frot = ftot = 0.0;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		Real dv = g->get_dx() * g->get_dx() * g->get_dx();
		for (int k = BW; k < GNX - BW; k++) {
			for (int j = BW; j < GNX - BW; j++) {
				for (int i = BW; i < GNX - BW; i++) {
					if (!g->zone_is_refined(i, j, k) && g->HydroGrid::xc(i) < 0.0) {
						d = (*g)(i, j, k).rho();
						dp = (*g)(i + 1, j, k).rho();
						dm = (*g)(i - 1, j, k).rho();
						hp = Kd * (n + 1.0) * pow(max(dp, 0.0), 1.0 / n);
						hm = Kd * (n + 1.0) * pow(max(dm, 0.0), 1.0 / n);
						phip = g->get_phi(i - o + 1, j - o, k - o);
						phim = g->get_phi(i - o - 1, j - o, k - o);
						tmp = d * g->HydroGrid::xc(i) * dv;
						ftot += d * ((hp - hm) + (phip - phim)) / (2.0 * g->get_dx()) * dv;
						frot += tmp;
					}
				}
			}
		}
	}
	tmp = ftot;
	MPI_Allreduce(&tmp, &ftot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = frot;
	MPI_Allreduce(&tmp, &frot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	//printf("%e %e\n", frot, ftot);
	Real oldO = State::get_omega();
	Real newO = sqrt(ftot / frot);
	State::set_omega(oldO * pow(newO / oldO, 0.05));
}

void BinaryStar::find_l(Real m1_x, Real m2_x, Real* l1_phi, Real* l1_x, int lnum) {
	BinaryStar* g;
	Vector<int, 3> loc;
	Real phi_max = -1.0e+99;
	int i, j, k;
	int plane_index;
	Real x, x_max;
	Real r[2], r0[2];
	bool test;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<BinaryStar*>(get_local_node(l));
		loc = g->get_location();
		plane_index = (1 << g->get_level()) / 2;
		if (loc[1] == plane_index && loc[2] == plane_index) {
			j = k = BW;
			for (i = BW; i < GNX - BW; i++) {
				if (!g->zone_is_refined(i, j, k)) {
					x = g->HydroGrid::xc(i);
					if (lnum == 1) {
						test = x > m2_x && x < m1_x;
					} else {
						test = x < m2_x;
					}
					if (test) {
						if ((*g)(i, j, k).phi_eff() > phi_max) {
							phi_max = (*g)(i, j, k).phi_eff();
							x_max = x;
						}
					}
				}
			}
		}
	}
	r[0] = phi_max;
	r[1] = x_max;
	MPI_Allreduce(r, r0, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD );
	*l1_phi = r0[0];
	*l1_x = r0[1];
}
#endif
#endif
