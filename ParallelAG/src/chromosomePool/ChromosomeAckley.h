/*
 * ChromosomeAckley.h
 *
 *  Created on: Oct 8, 2016
 *      Author: reginaldo
 */

#ifndef CHROMOSOMEPOOL_CHROMOSOMEACKLEY_H_
#define CHROMOSOMEPOOL_CHROMOSOMEACKLEY_H_

#include "../abstract/Chromosome.h"

class ChromosomeAckley: public Chromosome {
private:


public:
	ChromosomeAckley(double *);
	virtual ~ChromosomeAckley();

	static ChromosomeAckley *createChromosome();

	void evaluate();
};

#endif /* CHROMOSOMEPOOL_CHROMOSOMEACKLEY_H_ */
