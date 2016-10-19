/*
 * MathFunctions.h
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#ifndef FUNCTIONS_MATHFUNCTIONS_H_
#define FUNCTIONS_MATHFUNCTIONS_H_
#include <math.h>
#include "../auxiliaries/Configuration.h"

class MathFunctions {
private:
	MathFunctions();
public:
	virtual ~MathFunctions();

	/**
	* Ackley's function. The minimum is f(0, 0) = 0.
	* Recommendation: vars should be [-4, 4].
	*/
	inline static double ackley(double *vars) {
		Ackley ack;
		double s = 0;

		for (int i = 0; i < ack.NGENES; i++)
			s += pow(vars[i], 2);
		double p1 = -20 * exp(-0.2 * sqrt(1.0 / ack.NGENES * s));

		s = 0;
		for (int i = 0; i < ack.NGENES; i++)
			s += cos(2 * M_PI * vars[i]);
		double p2 = -exp(1.0 / ack.NGENES * s) + 20 + M_E;
		return p1 + p2;
	}

};


#endif /* FUNCTIONS_MATHFUNCTIONS_H_ */
