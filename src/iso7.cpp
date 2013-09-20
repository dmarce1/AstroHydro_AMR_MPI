#include "iso7.h"
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

#define nrowmax 1000

__attribute__((constructor))
void initialize_iso7() {
	init_iso7_();
	read_helm_table_();
	net_initialize_();
}

void nuclear_burn(double dt, double density, double ein, double T, double xx[ISO7_N], double* Tout, double *dout, double* eout, double yy[ISO7_N]) {
	int nok, nbad;
	double conserv;
	double abar = 0.0, zbar = 0.0;
	double cs, p, tmp;
	ein /= density;
	burner_(&dt, &T, &density, &ein, xx, Tout, dout, eout, yy, &conserv, &nok, &nbad);
	*eout *= density;
}

void invert_eos_ed(double* T, double* cs, double* p, double d, double e, double abar, double zbar) {
	double dpdd, dpdt, dedt, dpde, cs2;
	int fail;
#ifdef STELLAR_STATE
	d = max(d, State::rho_floor);
	e = max(e, eos_energy(d, 2.0e+3, abar, zbar));
	invert_helm_ed_c_interface_(&d, T, &e, p, &abar, &zbar, &dpdd, &dpdt, &dedt, &fail);
	dpde = dpdt / dedt;
	cs2 = (*p / (d * d)) * dpde + dpdd;
	*cs = sqrt(cs2);
	if (fail == 1) {
		printf( "%e %e %e %e\n", d, *T, abar, zbar);
		abort();
	}
#endif
}

double eos_energy(double d, double T, double abar, double zbar) {
	double e, p;
	int fail;
	double dpdd, dpdt, dedt, dpde, cs2;
	helmeos_c_interface_(&d, &T, &e, &p, &abar, &zbar, &dpdd, &dpdt, &dedt, &fail);
	if (fail == 1) {
		printf( "%e %e %e %e\n", d, T, abar, zbar);
		abort();
	}
	return e;
}

double eos_pressure(double d, double T, double abar, double zbar) {
	double e, p;
	int fail;
	double dpdd, dpdt, dedt, dpde, cs2;
	helmeos_c_interface_(&d, &T, &e, &p, &abar, &zbar, &dpdd, &dpdt, &dedt, &fail);
	if (fail == 1) {
		printf( "%e %e %e %e\n", d, T, abar, zbar);
		abort();
	}
	return p;
}

double eos_soundspeed(double d, double T, double abar, double zbar) {
	double e, p, cs;
	double dpdd, dpdt, dedt, dpde, cs2;
	int fail;
	helmeos_c_interface_(&d, &T, &e, &p, &abar, &zbar, &dpdd, &dpdt, &dedt, &fail);
	dpde = dpdt / dedt;
	cs2 = (p / (d * d)) * dpde + dpdd;
	cs = sqrt(cs2);
	if (fail == 1) {
		printf( "%e %e %e %e\n", d, T, abar, zbar);
		abort();
	}
	return cs;
}

void eos(double* e, double* p, double* cs, double d, double T, double abar, double zbar) {
	double dpdd, dpdt, dedt, dpde, cs2;
	int fail;
	helmeos_c_interface_(&d, &T, e, p, &abar, &zbar, &dpdd, &dpdt, &dedt, &fail);
	dpde = dpdt / dedt;
	cs2 = (*p / (d * d)) * dpde + dpdd;
	*cs = sqrt(cs2);
	if (fail == 1) {
		printf( "%e %e %e %e\n", d, T, abar, zbar);
		abort();
	}
}

double nuclear_burn_dtmax(double d, double T, double abar, double zbar, double xx[ISO7_N]) {
	double e, p, dt;
	double dpdd, dpdt, dedt, dpde;
	int fail;
	helmeos_c_interface_(&d, &T, &e, &p, &abar, &zbar, &dpdd, &dpdt, &dedt, &fail);
	burn_dt_(&d, &e, &T, xx, &dt);
	return dt;
}
