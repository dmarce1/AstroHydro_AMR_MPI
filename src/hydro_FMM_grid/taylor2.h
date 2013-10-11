/*
 * taylor.h
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#ifndef TAYLOR2_H_
#define TAYLOR2_H_

#include "../defs.h"
#include "../real.h"

class Taylor2 {
private:
	Real l0;
	Real l1[3];
	Real l2[6];
public:
	Taylor2& operator=(const Taylor2& t);
	Taylor2& operator+=(const Taylor2& t);
	void zero();
	Taylor2();
	Real operator()() const;
	Real operator()(int i) const;
	Real operator()(int i, int j) const;
	Real& operator()();
	Real& operator()(int i);
	Real& operator()(int i, int j);
	Real operator()(int i, int j, int k) const;
	Real& operator()(int i, int j, int k);
	virtual ~Taylor2();
};
#endif /* TAYLOR_H_ */

