#ifndef FMM_H_
#define FMM_H_

#include "../vector.h"
#include "../real.h"
#include "../complex.h"
#include "../legendre.h"

#define FMM_POLE_MAX 2

class FMM {
private:
	static Real factorial[2 * (FMM_POLE_MAX + 1) + 1];
	static SphericalHarmonic Y;
	bool refined;
	int level;
	Real source;
	FMM* children[8];
	FMM* siblings[6];
	FMM* neighbors[6 + 12 + 8];
	FMM* interactions[216];
	int interaction_count;
	FMM* parent;
	_3Vec X, Xcom;
	Real dx;
	Complex* M[FMM_POLE_MAX + 1];
	Complex* L[FMM_POLE_MAX + 1];
	Complex mspace[(FMM_POLE_MAX + 1) * (FMM_POLE_MAX + 1)];
	Complex lspace[(FMM_POLE_MAX + 1) * (FMM_POLE_MAX + 1)];
	Real phi;
public:
	static Real A(int, int);
	static void init();
	void set_source(Real);
	FMM(Real dx = 1.0, int level = 0);
	virtual ~FMM();
	void refine();
	void derefine();
	void compute_interaction_list();
	void local_expansion();
	bool is_neighbor(const FMM*) const;
	void compute_moments();
	void compute_interactions();
	void compute_gravity();
	void find_neighbors();
	template<int N>
	static void form_tree(FMM* tree, Real rho[N][N][N]);
	static void run(int, char*[]);
	Real error() const;
};

#endif
