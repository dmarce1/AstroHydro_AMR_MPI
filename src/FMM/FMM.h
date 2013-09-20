#ifndef FMM_H_
#define FMM_H_

#include "../vector.h"
#include "../real.h"

#define FMM_POLE_MAX 1
#define TEST_RADIUS 0.25

class FMM {
private:
bool refined;
	int level;
	Real source;
	FMM* children[8];
	FMM* siblings[6];
	FMM* neighbors[6 + 12 + 8];
	FMM* interactions[189];
	int interaction_count;
	FMM* parent;
	_3Vec X, Xcom;
	Real dx;
	double M0;
	double M2[3][3];
	double C0;
	double C1[3];
	double C2[3][3];
	double C3[3][3][3];
	/*Complex* M[FMM_POLE_MAX + 1];
	Complex* L[FMM_POLE_MAX + 1];
	Complex mspace[(FMM_POLE_MAX + 1) * (FMM_POLE_MAX + 1)];
	Complex lspace[(FMM_POLE_MAX + 1) * (FMM_POLE_MAX + 1)];*/
	Real phi;
	_3Vec g;
public:
	Vector<Real,3> dS;
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
	Real error() ;
};

#endif
