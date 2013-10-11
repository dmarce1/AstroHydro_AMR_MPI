/*
 * Moment.h
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#ifndef MOMENT_H_
#define MOMENT_H_
#define USE_HIGH_ORDER_POT

#include "../defs.h"
#include "../vector.h"
#include "../real.h"

class Moment {
private:
	Real m;
	Real m2[6];
#ifdef USE_HIGH_ORDER_POT
	Real m3[10];
#endif
public:
	bool is_leaf;
	Moment& operator=(const Moment& m);
	_3Vec X;
	Moment();
	Real M() const;
	Real M2(int i, int j) const;
	Real& M();
	Real& M2(int i, int j);
#ifdef USE_HIGH_ORDER_POT
	Real M3(int i, int j, int k) const;
	Real& M3(int i, int j, int k);
#endif
	virtual ~Moment();
};

#endif /* MOMENT_H_ */
