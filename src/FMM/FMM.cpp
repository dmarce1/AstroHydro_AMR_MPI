#include <stdlib.h>
#include <stdio.h>
#include "FMM.h"

#define NCHILD 8
#define NFACE 6
#define NEDGE 12
#define NVERTEX 8
#define NNEIGHBORS (NFACE+NVERTEX+NEDGE)
#define NINTERACTIONS 189

#define X0 0.2
#define Y0 0.00
#define Z0 0.00

const Real kD[3][3] = { { 1., 0, 0 }, { 0, 1., 0 }, { 0, 0, 1. } };

void FMM::compute_moments() {
	_3Vec dX;
	for (int n = 0; n < 3; n++) {
		for (int m = 0; m < 3; m++) {
			M2[n][m] = 0.0;
		}
	}
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->compute_moments();
		}
		Xcom = 0.0;
		M0 = 0.0;
		for (int i = 0; i < NCHILD; i++) {
			Xcom += children[i]->Xcom * children[i]->M0;
			M0 += children[i]->M0;
		}
		if (M0 == 0.0) {
			Xcom = X;
		} else {
			Xcom /= M0;
			for (int i = 0; i < NCHILD; i++) {
				dX = Xcom - children[i]->Xcom;
				for (int n = 0; n < 3; n++) {
					for (int m = 0; m < 3; m++) {
						M2[n][m] += children[i]->M2[n][m];
						M2[n][m] += dX[n] * dX[m] * children[i]->M0;
					}
				}
			}
		}
	} else {
		M0 = source;
		Xcom = X;
	}
}
void FMM::compute_interactions() {
	FMM* I;
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->compute_interactions();
		}
	}
	C0 = 0.0;
	for (int n = 0; n < 3; n++) {
		C1[n] = 0.0;
		for (int m = 0; m < 3; m++) {
			C2[n][m] = 0.0;
			for (int k = 0; k < 3; k++) {
				C3[n][m][k] = 0.0;
			}
		}
	}
	for (int i = 0; i < interaction_count; i++) {
		I = interactions[i];
		_3Vec R = Xcom - I->Xcom;
		const Real rinv = 1.0 / sqrt(R.dot(R));
		const Real rinv2 = rinv * rinv;
		const Real rinv3 = rinv2 * rinv;
		const Real rinv5 = rinv2 * rinv3;
		const Real rinv7 = rinv2 * rinv5;
		const Real d0 = -rinv;
		const Real d1 = +rinv3;
		const Real d2 = -3.0 * rinv5;
		const Real d3 = +15.0 * rinv7;
		C0 += I->M0 * d0;
		for (int n = 0; n < 3; n++) {
			C0 += I->M2[n][n] * d1 * 0.5;
			C1[n] += I->M0 * d1 * R[n];
			for (int m = 0; m < 3; m++) {
				C0 += I->M2[n][m] * d2 * 0.5 * R[n] * R[m];
				C1[n] += I->M2[m][m] * d2 * 0.5 * R[n];
				C1[n] += I->M2[n][m] * d2 * 1.0 * R[m];
				C2[n][m] += I->M0 * kD[n][m] * d1;
				C2[n][m] += I->M0 * R[n] * R[m] * d2;
				for (int q = 0; q < 3; q++) {
					C1[n] += I->M2[m][q] * d3 * 0.5 * R[n] * R[m] * R[q];
					C3[n][m][q] += I->M0 * d2 * kD[n][m] * R[q];
					C3[n][m][q] += I->M0 * d2 * kD[m][q] * R[n];
					C3[n][m][q] += I->M0 * d2 * kD[q][n] * R[m];
					C3[n][m][q] += I->M0 * d3 * R[n] * R[m] * R[q];
				}
			}
		}
	}
}

void FMM::local_expansion() {
	if (parent != NULL) {
		FMM* p = parent;
		_3Vec Y = Xcom - p->Xcom;
		C0 += p->C0;
		for (int n = 0; n < 3; n++) {
			C0 += Y[n] * p->C1[n];
			C1[n] += p->C1[n];
			for (int m = 0; m < 3; m++) {
				C0 += Y[n] * Y[m] * p->C2[n][m] * 0.5;
				C1[n] += Y[m] * p->C2[n][m];
				C2[n][m] += p->C2[n][m];
				for (int q = 0; q < 3; q++) {
					C0 += Y[n] * Y[m] * Y[q] * p->C3[n][m][q] * (1.0 / 6.0);
					C1[n] += Y[m] * Y[q] * p->C3[n][m][q] * 0.5;
					C2[n][m] += Y[q] * p->C3[n][m][q];
					C3[n][m][q] += p->C3[n][m][q];
				}
			}
		}
	}
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->local_expansion();
		}
	}
}

void FMM::compute_gravity() {
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->compute_gravity();
		}
	}
	phi = C0;
	for (int i = 0; i < 3; i++) {
		g[i] = -C1[i];
	}
}

template<int N>
void FMM::form_tree(FMM* tree, Real rho[N][N][N]) {
	tree->refine();
	bool hasrho;
	static Real rho1[N / 2][N / 2][N / 2];
	for (int ci = 0; ci < NCHILD; ci++) {
		int oi = ((ci & 1) >> 0) * (N / 2);
		int oj = ((ci & 2) >> 1) * (N / 2);
		int ok = ((ci & 4) >> 2) * (N / 2);
		hasrho = false;
		for (int i = 0; i < N / 2; i++) {
			for (int j = 0; j < N / 2; j++) {
				for (int k = 0; k < N / 2; k++) {
					rho1[i][j][k] = rho[i + oi][j + oj][k + ok];
					if (rho1[i][j][k] > 0.0) {
						hasrho = true;
					}
				}
			}
		}
		//	hasrho = true;
		if (hasrho) {
			form_tree<N / 2>(tree->children[ci], rho1);
		} else {
			tree->children[ci]->source = 0.0;
		}
	}
}

template<>
void FMM::form_tree<1>(FMM* tree, Real rho[1][1][1]) {
	tree->source = rho[0][0][0] * pow(tree->dx, 3);
}

double ___mtot;

void FMM::run(int, char*[]) {
#define NX 128
	FMM* root;
	static Real rho[NX][NX][NX];
	Real x, y, z, r, m;
	m = 0.0;
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NX; j++) {
			for (int k = 0; k < NX; k++) {
				x = (Real(i - NX / 2) / Real(NX / 2) + 1.0 / NX) * 0.5;
				y = (Real(j - NX / 2) / Real(NX / 2) + 1.0 / NX) * 0.5;
				z = (Real(k - NX / 2) / Real(NX / 2) + 1.0 / NX) * 0.5;
				x -= X0;
				y -= Y0;
				z -= Z0;
				r = sqrt(x * x + y * y + z * z);
				if (r < TEST_RADIUS) {
					//		const double dx = 0.5 / NX;
					const double v0 = (4.0 / 3.0) * M_PI * pow(TEST_RADIUS, 3);
					rho[i][j][k] = 1.0 / v0;
				} else {
					rho[i][j][k] = 0.0;
				}
				m += rho[i][j][k] * pow(1.0 / NX, 3);
			}
		}
	}
	___mtot = m;
	printf("mass = %e\n", m);
	root = new FMM();
	printf("Form tree\n");
	form_tree<NX>(root, rho);
	printf("Find neighbors\n");
	root->find_neighbors();
	printf("Compute interaction list\n");
	root->compute_interaction_list();
	printf("Compute moments\n");
	root->compute_moments();
	printf("Compute interaction\n");
	root->compute_interactions();
	printf("Local expansion\n");
	root->local_expansion();
	printf("Compute gravity\n");
	root->compute_gravity();
	Vector<Real, 4> e = root->error();
	printf("Error is %e, %e, %e, %e \n", e[0], e[1], e[2], e[3]);
	printf("dMomentum is %e %e %e\n", root->dS[0], root->dS[1], root->dS[2]);
}

Vector<Real, 4> FMM::error() {
	Real a[4], n, m;
	Vector<Real, 4> e;
	if (refined) {
		e = 0.0;
		dS = 0.0;
		for (int i = 0; i < NCHILD; i++) {
			e += children[i]->error();
			dS += children[i]->dS;
		}
	} else {
		m = 1.0;
		X[0] -= X0;
		X[1] -= Y0;
		X[2] -= Z0;
		Real r = sqrt(X.dot(X));
		if (r < TEST_RADIUS) {
			a[0] = -1.0 / (2.0 * TEST_RADIUS) * (3.0 - r * r / TEST_RADIUS / TEST_RADIUS);
			a[1] = -X[0] / (TEST_RADIUS * TEST_RADIUS * TEST_RADIUS);
			a[2] = -X[1] / (TEST_RADIUS * TEST_RADIUS * TEST_RADIUS);
			a[3] = -X[2] / (TEST_RADIUS * TEST_RADIUS * TEST_RADIUS);
		} else {
			a[0] = -1.0 / r;
			a[1] = -X[0] / (r * r * r);
			a[2] = -X[1] / (r * r * r);
			a[3] = -X[2] / (r * r * r);
		}
		a[0] *= ___mtot;
		a[1] *= ___mtot;
		a[2] *= ___mtot;
		a[3] *= ___mtot;
		n = phi;
		e[0] = fabs(log(n / a[0])) * dx * dx * dx;
		e[1] = fabs(g[0] - a[1]) * dx * dx * dx;
		e[2] = fabs(g[1] - a[2]) * dx * dx * dx;
		e[3] = fabs(g[2] - a[3]) * dx * dx * dx;
		dS = g * source;
		if (X[2] < dx && X[1] < dx && X[2] > 0.0 * dx && X[1] > 0.0 * dx) {
			printf("%e %e %e %e %e %e %e %e %e %e %e %e\n", X[0], X[1], X[2], a[0], n, source, g[0], a[1], g[1], a[2], g[2], a[3]);
		}
	}
	return e;
}

FMM::FMM(Real _dx, int lev) {
	init();
	level = lev;
	dx = _dx;
	refined = false;
	parent = NULL;
	X[0] = X[1] = X[2] = 0.0;
	for (int i = 0; i < NCHILD; i++) {
		children[i] = NULL;
	}
	for (int i = 0; i < NFACE; i++) {
		siblings[i] = NULL;
	}
	for (int i = 0; i < NINTERACTIONS; i++) {
		interactions[i] = NULL;
	}
	for (int i = 0; i < NNEIGHBORS; i++) {
		neighbors[i] = NULL;
	}
	int n = 0;
}

void FMM::refine() {
	refined = true;
//	printf( "%e\n", dx);
	for (int i = 0; i < NCHILD; i++) {
		children[i] = new FMM(dx * 0.5, level + 1);
		children[i]->parent = this;
	}
	for (int i = 0; i < NCHILD; i++) {
		_3Vec xc = X;
		xc[0] += 0.25 * dx * (((i & 1) != 0) ? 1.0 : -1.0);
		xc[1] += 0.25 * dx * (((i & 2) != 0) ? 1.0 : -1.0);
		xc[2] += 0.25 * dx * (((i & 4) != 0) ? 1.0 : -1.0);
		children[i]->X = xc;
	}
	for (int i = 0; i < NCHILD; i++) {
		for (int mask = 1, dir = 0; dir < 3; mask <<= 1, dir++) {
			const int oface = (dir << 1) + ((i & mask) ? 1 : 0);
			const int nface = (dir << 1) + ((i & mask) ? 0 : 1);
			if (siblings[oface] != NULL) {
				children[i]->siblings[oface] = siblings[oface]->children[i ^ mask];
				if (siblings[oface]->refined) {
					siblings[oface]->children[i ^ mask]->siblings[nface] = children[i];
				}
			} else {
				children[i]->siblings[oface] = NULL;
			}
			children[i]->siblings[nface] = children[i ^ mask];
		}
	}
}

#define XL 0
#define XU 1
#define YL 2
#define YU 3
#define ZL 4
#define ZU 5

void FMM::find_neighbors() {
	FMM* ptr;
	int i = 0;
	neighbors[i++] = siblings[XL];
	neighbors[i++] = siblings[XU];
	neighbors[i++] = siblings[YL];
	neighbors[i++] = siblings[YU];
	neighbors[i++] = siblings[ZL];
	neighbors[i++] = siblings[ZU];
	if (siblings[XL] != NULL) {
		neighbors[i++] = siblings[XL]->siblings[YL];
	} else if (siblings[YL] != NULL) {
		neighbors[i++] = siblings[YL]->siblings[XL];
	}
	if (siblings[XU] != NULL) {
		neighbors[i++] = siblings[XU]->siblings[YU];
	} else if (siblings[YU] != NULL) {
		neighbors[i++] = siblings[YU]->siblings[XU];
	}
	if (siblings[XU] != NULL) {
		neighbors[i++] = siblings[XU]->siblings[YL];
	} else if (siblings[YL] != NULL) {
		neighbors[i++] = siblings[YL]->siblings[XU];
	}
	if (siblings[XL] != NULL) {
		neighbors[i++] = siblings[XL]->siblings[YU];
	} else if (siblings[YU] != NULL) {
		neighbors[i++] = siblings[YU]->siblings[XL];
	}

	if (siblings[XL] != NULL) {
		neighbors[i++] = siblings[XL]->siblings[ZL];
	} else if (siblings[ZL] != NULL) {
		neighbors[i++] = siblings[ZL]->siblings[XL];
	}
	if (siblings[XU] != NULL) {
		neighbors[i++] = siblings[XU]->siblings[ZU];
	} else if (siblings[ZU] != NULL) {
		neighbors[i++] = siblings[ZU]->siblings[XU];
	}
	if (siblings[XU] != NULL) {
		neighbors[i++] = siblings[XU]->siblings[ZL];
	} else if (siblings[ZL] != NULL) {
		neighbors[i++] = siblings[ZL]->siblings[XU];
	}
	if (siblings[XL] != NULL) {
		neighbors[i++] = siblings[XL]->siblings[ZU];
	} else if (siblings[ZU] != NULL) {
		neighbors[i++] = siblings[ZU]->siblings[XL];
	}

	if (siblings[ZL] != NULL) {
		neighbors[i++] = siblings[ZL]->siblings[YL];
	} else if (siblings[YL] != NULL) {
		neighbors[i++] = siblings[YL]->siblings[ZL];
	}
	if (siblings[ZU] != NULL) {
		neighbors[i++] = siblings[ZU]->siblings[YU];
	} else if (siblings[YU] != NULL) {
		neighbors[i++] = siblings[YU]->siblings[ZU];
	}
	if (siblings[ZU] != NULL) {
		neighbors[i++] = siblings[ZU]->siblings[YL];
	} else if (siblings[YL] != NULL) {
		neighbors[i++] = siblings[YL]->siblings[ZU];
	}
	if (siblings[ZL] != NULL) {
		neighbors[i++] = siblings[ZL]->siblings[YU];
	} else if (siblings[YU] != NULL) {
		neighbors[i++] = siblings[YU]->siblings[ZL];
	}

	if (siblings[XL] != NULL) {
		if (siblings[XL]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[YL]->siblings[ZL];
		} else if (siblings[XL]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[ZL]->siblings[YL];
		}
	} else if (siblings[YL] != NULL) {
		if (siblings[YL]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[XL]->siblings[ZL];
		} else if (siblings[YL]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[ZL]->siblings[XL];
		}
	} else if (siblings[ZL] != NULL) {
		if (siblings[ZL]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[YL]->siblings[XL];
		} else if (siblings[ZL]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[XL]->siblings[YL];
		}
	}

	if (siblings[XU] != NULL) {
		if (siblings[XU]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[YL]->siblings[ZL];
		} else if (siblings[XU]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[ZL]->siblings[YL];
		}
	} else if (siblings[YL] != NULL) {
		if (siblings[YL]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[XU]->siblings[ZL];
		} else if (siblings[YL]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[ZL]->siblings[XU];
		}
	} else if (siblings[ZL] != NULL) {
		if (siblings[ZL]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[YL]->siblings[XU];
		} else if (siblings[ZL]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[XU]->siblings[YL];
		}
	}

	if (siblings[XL] != NULL) {
		if (siblings[XL]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[YU]->siblings[ZL];
		} else if (siblings[XL]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[ZL]->siblings[YU];
		}
	} else if (siblings[YU] != NULL) {
		if (siblings[YU]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[XL]->siblings[ZL];
		} else if (siblings[YU]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[ZL]->siblings[XL];
		}
	} else if (siblings[ZL] != NULL) {
		if (siblings[ZL]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[YU]->siblings[XL];
		} else if (siblings[ZL]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[XL]->siblings[YU];
		}
	}

	if (siblings[XU] != NULL) {
		if (siblings[XU]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[YU]->siblings[ZL];
		} else if (siblings[XU]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[ZL]->siblings[YU];
		}
	} else if (siblings[YU] != NULL) {
		if (siblings[YU]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[XU]->siblings[ZL];
		} else if (siblings[YU]->siblings[ZL] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[ZL]->siblings[XU];
		}
	} else if (siblings[ZL] != NULL) {
		if (siblings[ZL]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[YU]->siblings[XU];
		} else if (siblings[ZL]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[ZL]->siblings[XU]->siblings[YU];
		}
	}

	if (siblings[XL] != NULL) {
		if (siblings[XL]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[YL]->siblings[ZU];
		} else if (siblings[XL]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[ZU]->siblings[YL];
		}
	} else if (siblings[YL] != NULL) {
		if (siblings[YL]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[XL]->siblings[ZU];
		} else if (siblings[YL]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[ZU]->siblings[XL];
		}
	} else if (siblings[ZU] != NULL) {
		if (siblings[ZU]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[YL]->siblings[XL];
		} else if (siblings[ZU]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[XL]->siblings[YL];
		}
	}

	if (siblings[XU] != NULL) {
		if (siblings[XU]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[YL]->siblings[ZU];
		} else if (siblings[XU]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[ZU]->siblings[YL];
		}
	} else if (siblings[YL] != NULL) {
		if (siblings[YL]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[XU]->siblings[ZU];
		} else if (siblings[YL]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[YL]->siblings[ZU]->siblings[XU];
		}
	} else if (siblings[ZU] != NULL) {
		if (siblings[ZU]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[YL]->siblings[XU];
		} else if (siblings[ZU]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[XU]->siblings[YL];
		}
	}

	if (siblings[XL] != NULL) {
		if (siblings[XL]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[YU]->siblings[ZU];
		} else if (siblings[XL]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[ZU]->siblings[YU];
		}
	} else if (siblings[YU] != NULL) {
		if (siblings[YU]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[XL]->siblings[ZU];
		} else if (siblings[YU]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[ZU]->siblings[XL];
		}
	} else if (siblings[ZU] != NULL) {
		if (siblings[ZU]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[YU]->siblings[XL];
		} else if (siblings[ZU]->siblings[XL] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[XL]->siblings[YU];
		}
	}

	if (siblings[XU] != NULL) {
		if (siblings[XU]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[YU]->siblings[ZU];
		} else if (siblings[XU]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[ZU]->siblings[YU];
		}
	} else if (siblings[YU] != NULL) {
		if (siblings[YU]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[XU]->siblings[ZU];
		} else if (siblings[YU]->siblings[ZU] != NULL) {
			neighbors[i++] = siblings[YU]->siblings[ZU]->siblings[XU];
		}
	} else if (siblings[ZU] != NULL) {
		if (siblings[ZU]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[YU]->siblings[XU];
		} else if (siblings[ZU]->siblings[XU] != NULL) {
			neighbors[i++] = siblings[ZU]->siblings[XU]->siblings[YU];
		}
	}

	for (int j = i; j < NNEIGHBORS; j++) {
		neighbors[j] = NULL;
	}

	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->find_neighbors();
		}
	}
}

void FMM::compute_interaction_list() {
	FMM* ptr;
	if (parent != NULL) {
		interaction_count = 0;
		for (int i = 0; i < NNEIGHBORS; i++) {
			if (parent->neighbors[i] != NULL) {
				if (parent->neighbors[i]->refined) {
					for (int j = 0; j < NCHILD; j++) {
						ptr = parent->neighbors[i]->children[j];
						if (!is_neighbor(ptr) || !refined || !ptr->refined) {
							interactions[interaction_count++] = ptr;
						}
					}
				}
			}
		}
		for (int i = 0; i < NCHILD; i++) {
			ptr = parent->children[i];
			if (ptr != this && (!refined || !ptr->refined)) {
				interactions[interaction_count++] = ptr;
			}
		}
	}
	/********************************************************************************************/
	if (refined) {
		for (int ci = 0; ci < 8; ci++) {
			children[ci]->compute_interaction_list();
		}
	}
	/********************************************************************************************/
}

bool FMM::is_neighbor(const FMM* f) const {
	_3Vec dist;
	bool near;
	assert( f->level == level);
	dist = f->X - X;
	near = false;
	near = ((fabs(dist[0]) < 1.5 * dx) && (fabs(dist[1]) < 1.5 * dx) && (fabs(dist[2]) < 1.5 * dx));
	return near;
}

FMM::~FMM() {
}

void FMM::init() {
	static bool initialized = false;
	if (!initialized) {
	}
}

void FMM::derefine() {
	refined = false;
	for (int i = 0; i < NCHILD; i++) {
		for (int mask = 1, dir = 0; dir < 3; mask <<= 1, dir++) {
			const int oface = (dir << 1) + ((i & mask) ? 1 : 0);
			const int nface = (dir << 1) + ((i & mask) ? 0 : 1);
			if (siblings[oface] != NULL) {
				if (siblings[oface]->refined) {
					siblings[oface]->children[i ^ mask]->siblings[nface] = NULL;
				}
			}
		}
	}
	for (int i = 0; i < NCHILD; i++) {
		delete children[i];
	}
}
