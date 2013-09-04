#include <stdlib.h>
#include <stdio.h>
#include "FMM.h"

#define NCHILD 8
#define NFACE 6
#define NEDGE 12
#define NVERTEX 8
#define NNEIGHBORS (NFACE+NVERTEX+NEDGE)
#define NINTERACTIONS 189

Real FMM::factorial[2 * (FMM_POLE_MAX + 1) + 1];
SphericalHarmonic FMM::Y(FMM_POLE_MAX + 1);

Real neg_one_pow(int i) {
	return Real(1 - 2 * (i & 1));
}

Real FMM::A(int n, int m) {
	return neg_one_pow(n) / sqrt(factorial[n - m] * factorial[n + m]);
}

void FMM::compute_moments() {
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->compute_moments();
		}
		Complex tmp;
		Real r;
		_3Vec dX;
		Xcom = 0.0;
		source = 0.0;
		for (int i = 0; i < NCHILD; i++) {
			Xcom += children[i]->Xcom * children[i]->source;
			source += children[i]->source;
		}
		if (source == 0.0) {
			Xcom = X;
		} else {
			Xcom /= source;
		}
		for (int j = 0; j <= FMM_POLE_MAX; j++) {
			for (int k = -j; k <= j; k++) {
				M[j][k] = 0.0;
				for (int i = 0; i < NCHILD; i++) {
					dX = children[i]->Xcom - Xcom;
					r = sqrt(dX.dot(dX));
					if (r != 0.0) {
						Y.generate(dX);
						for (int n = 0; n <= j; n++) {
							for (int m = -n; m <= +n; m++) {
								if (abs(k - m) <= j - n) {
									tmp = children[i]->M[j - n][k - m] * A(n, m) * A(j - n, k - m) * pow(r, n) * Y(n, -m) / A(j, k);
									M[j][k] += Complex::ipow(abs(k) - abs(m) - abs(k - m)) * tmp;
								}
							}
						}
					} else {
						M[j][k] += children[i]->M[j][k];
					}
				}
			}
		}
	} else {
		Xcom = X;
		M[0][0].set_real(-source);
		M[0][0].set_imag(0.0);
		for (int l = 1; l <= FMM_POLE_MAX; l++) {
			for (int m = -l; m <= l; m++) {
				M[l][m] = 0.0;
			}
		}
	}
}

void FMM::compute_interactions() {
	Complex tmp;
	Real r;
	_3Vec dX;
	for (int j = 0; j <= FMM_POLE_MAX; j++) {
		for (int k = -j; k <= j; k++) {
			L[j][k] = 0.0;
			for (int i = 0; i < interaction_count; i++) {
				dX = interactions[i]->Xcom - Xcom;
				Y.generate(dX);
				r = sqrt(dX.dot(dX));
				for (int n = 0; n <= FMM_POLE_MAX - j; n++) {
					for (int m = -n; m <= +n; m++) {
						tmp = interactions[i]->M[n][m] * A(n, m) * A(j, k) * pow(r, -(j + n + 1)) * Y(j + n, m - k) / A(j + n, m - k);
						L[j][k] += Complex::ipow(abs(k - m) - abs(k) - abs(m)) * neg_one_pow(n) * tmp;
					}
				}

			}
		}
	}
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->compute_interactions();
		}
	}

}

void FMM::local_expansion() {
	Complex tmp;
	Real r;
	_3Vec dX;
	if (parent != NULL) {
		for (int j = 0; j <= FMM_POLE_MAX; j++) {
			for (int k = -j; k <= j; k++) {
				dX = parent->Xcom - Xcom;
				r = sqrt(dX.dot(dX));
				if (r != 0.0) {
					Y.generate(dX);
					for (int n = j; n <= FMM_POLE_MAX; n++) {
						for (int m = -n; m <= +n; m++) {
							if (abs(m - k) <= n - j) {
								tmp = parent->L[n][m] * A(j, k) * A(n - j, m - k) * pow(r, n - j) * Y(n - j, m - k) / A(n, m);
								L[j][k] += Complex::ipow(abs(m) - abs(k - m) - abs(k)) * neg_one_pow(n + j) * tmp;
							}
						}
					}
				} else {
					L[j][k] += parent->L[j][k];
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
	Real direct, r;
	_3Vec dist;
	Complex direct_complex;
	if (!refined) {
		direct = 0.0;
		for (int i = 0; i < NNEIGHBORS; i++) {
			if (neighbors[i] != NULL) {
				dist = X - neighbors[i]->X;
				r = sqrt(dist.dot(dist));
				direct += -neighbors[i]->source / r;
			}
		}
		direct_complex.set_imag(0.0);
		direct_complex.set_real(direct);
		L[0][0] += direct_complex;
	} else {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->compute_gravity();
		}
	}
	phi = L[0][0].real();
}

template<int N>
void FMM::form_tree(FMM* tree, Real rho[N][N][N]) {
	tree->refine();
	static Real rho1[N / 2][N / 2][N / 2];
	for (int ci = 0; ci < NCHILD; ci++) {
		int oi = ((ci & 1) >> 0) * (N / 2);
		int oj = ((ci & 2) >> 1) * (N / 2);
		int ok = ((ci & 4) >> 2) * (N / 2);
		for (int i = 0; i < N / 2; i++) {
			for (int j = 0; j < N / 2; j++) {
				for (int k = 0; k < N / 2; k++) {
					rho1[i][j][k] = rho[i + oi][j + oj][k + ok];
				}
			}
		}
		form_tree<N / 2>(tree->children[ci], rho1);
	}
}

template<>
void FMM::form_tree<1>(FMM* tree, Real rho[1][1][1]) {
	tree->source = rho[0][0][0] * pow(tree->dx, 3);
}

void FMM::run(int, char*[]) {
#define NX 32
	FMM* root;
	static Real rho[NX][NX][NX];
	Real x, y, z, r, m;
	m = 0.0;
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NX; j++) {
			for (int k = 0; k < NX; k++) {
				x = Real(i - NX / 2) / Real(NX / 2);
				y = Real(j - NX / 2) / Real(NX / 2);
				z = Real(k - NX / 2) / Real(NX / 2);
				r = sqrt(x * x + y * y + z * z);
				if (r < 1.0 / NX) {
					rho[i][j][k] = pow(1.0 / NX, -3);
				} else {
					rho[i][j][k] = 0.0;
				}
				m += rho[i][j][k] * pow(1.0 / NX, 3);
			}
		}
	}
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
	printf("Error is %e\n", root->error());
}

Real FMM::error() const {
	Real e, a, n, m;
	if (refined) {
		e = 0.0;
		for (int i = 0; i < NCHILD; i++) {
			e += children[i]->error();
		}
	} else {
		m = 1.0;
		Real r = sqrt(X.dot(X));
		if (r < dx) {
			a = -3.0 / 2.0 / dx;
		} else {
			a = -1.0 / r;
		}
		n = phi;
		e = pow(1.0 - n / a, 2) * pow(dx, 3);
		if (fabs(X[2]) < dx) {
			printf("%e %e %e %e %e %e\n", Xcom[0], Xcom[1], Xcom[2], a, n, source);
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
	for (int l = 0; l <= FMM_POLE_MAX; l++) {
		M[l] = mspace + n + l;
		L[l] = lspace + n + l;
		n += 2 * l + 1;
	}
}

void FMM::refine() {
	refined = true;
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
		neighbors[i++] = siblings[XL]->siblings[YU];
		neighbors[i++] = siblings[XL]->siblings[ZL];
		neighbors[i++] = siblings[XL]->siblings[ZU];
		if (siblings[XL]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[YL]->siblings[ZL];
			neighbors[i++] = siblings[XL]->siblings[YL]->siblings[ZU];
		}
		if (siblings[XL]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[XL]->siblings[YU]->siblings[ZL];
			neighbors[i++] = siblings[XL]->siblings[YU]->siblings[ZU];
		}
	}
	if (siblings[XU] != NULL) {
		neighbors[i++] = siblings[XU]->siblings[YL];
		neighbors[i++] = siblings[XU]->siblings[YU];
		neighbors[i++] = siblings[XU]->siblings[ZL];
		neighbors[i++] = siblings[XU]->siblings[ZU];
		if (siblings[XU]->siblings[YL] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[YL]->siblings[ZL];
			neighbors[i++] = siblings[XU]->siblings[YL]->siblings[ZU];
		}
		if (siblings[XU]->siblings[YU] != NULL) {
			neighbors[i++] = siblings[XU]->siblings[YU]->siblings[ZL];
			neighbors[i++] = siblings[XU]->siblings[YU]->siblings[ZU];
		}
	}
	if (siblings[YL] != NULL) {
		neighbors[i++] = siblings[YL]->siblings[ZL];
		neighbors[i++] = siblings[YL]->siblings[ZU];
	}
	if (siblings[YU] != NULL) {
		neighbors[i++] = siblings[YU]->siblings[ZL];
		neighbors[i++] = siblings[YU]->siblings[ZU];
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
						if (!is_neighbor(ptr)) {
							interactions[interaction_count++] = ptr;
						}
					}
				}
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
		factorial[0] = 1;
		for (int l = 1; l <= 2 * FMM_POLE_MAX; l++) {
			factorial[l] = l * factorial[l - 1];
		}
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
