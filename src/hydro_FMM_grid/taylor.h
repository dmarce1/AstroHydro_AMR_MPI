/*
 * taylor.h
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#ifndef TAYLOR_H_
#define TAYLOR_H_

#include "../defs.h"
#include "../real.h"

class Taylor {
private:
	Real l0;
	Real l1[3];
	Real l2[6];
	Real l3[10];
	Real l4[15];
public:
	Taylor& operator=(const Taylor& t);
	Taylor& operator+=(const Taylor& t);
	void zero();
	Taylor();
	Real operator()() const;
	Real operator()(int i) const;
	Real operator()(int i, int j) const;
	Real& operator()();
	Real& operator()(int i);
	Real& operator()(int i, int j);
	Real operator()(int i, int j, int k) const;
	Real& operator()(int i, int j, int k);
	virtual ~Taylor();
};

#endif /* TAYLOR_H_ */
