/*
 * Configuration.h
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#ifndef AUXILIARIES_CONFIGURATION_H_
#define AUXILIARIES_CONFIGURATION_H_

/*Genetic algorithm parameters*/
const int NCHROMOSOMES = 500; // must be even ('cause of crossover)
const int NGENERATIONS = 1000;
const double CROSSOVERRATE = 0.8;
const double MUTATIONRATE = 0.01;
const double SD = 0.1;
const int RANK = 3;

/*Output format parameters*/
const int PRECISION = 3;
const int WIDTH = PRECISION + 4;

/*Thread parameters*/
const int NTHREADS = 2;

/*Auxiliary classes*/
class Ackley {
public:
	static const int ID = 0;
	static const int NGENES = 500;
	static const int INFIMUM = -40;
	static const int MAXIMUM = 40;
};


#endif /* AUXILIARIES_CONFIGURATION_H_ */
