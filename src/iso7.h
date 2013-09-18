#ifndef ISOISO7_N_______GH
#define ISOISO7_N_______GH

#define ISO7_N 7

#include "physical_constants.h"

typedef struct {
	const char* name;
	double A;
	int Z;

} iso7_t;

//static iso7_t Iso7topes[ISO7_N] = { { "He", 4.002602, 2 }, { "C", 12.0107, 6 }, { "O", 15.9994, 8 }, { "Ne", 20.1797, 10 }, { "Mg", 24.3050, 12 },
//		{ "Si", 28.0855, 14 }, { "Ni", 58.6934, 28 } };

static iso7_t Iso7topes[ISO7_N] = { { "He", 4.0, 2 }, { "C", 12.0, 6 }, { "O", 16.0, 8 }, { "Ne", 20.0, 10 }, { "Mg", 24.0, 12 }, { "Si", 28.0, 14 }, { "Ni",
		56.0, 28 } };

extern "C" {
void read_helm_table_(void);
void burner_(double* tstep, double* tin, double* din, double* ein, double xin[ISO7_N], double* tout, double* dout, double* eout, double xout[ISO7_N],
		double* conserv, int* nok, int* nbad);
void net_initialize_(void);
void init_iso7_(void);
void invert_helm_ed_c_interface_(const double* d, double* T, const double* e, double* p, const double* abar, const double* zbar, double* dpdd, double* dpdt,
		double* dedt, int* failure);
void helmeos_c_interface_(const double* d, const double* T, double* e, double* p, const double* abar, const double* zbar, double* dp_dd, double* dp_dt,
		double* de_dt, int* fail);
void burn_dt_(const double* din, const double* ein, const double* tin, double xin[ISO7_N], double* dt);
}

void nuclear_burn(double dt, double density, double ein, double T, double xx[ISO7_N], double* Tout, double *dout, double* eout, double yy[ISO7_N]);
void invert_eos_ed(double* T, double* cs, double* p, double d, double e, double abar, double zbar);
double eos_energy(double d, double T, double abar, double zbar);
double eos_pressure(double d, double T, double abar, double zbar);
double eos_soundspeed(double d, double T, double abar, double zbar);
void eos(double* e, double* p, double* cs, double d, double T, double abar, double zbar);
double nuclear_burn_dtmax(double d, double T, double abar, double zbar, double xx[ISO7_N]);

#endif
