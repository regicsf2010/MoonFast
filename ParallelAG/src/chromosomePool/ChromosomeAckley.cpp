/*
 * ChromosomeAckley.cpp
 *
 *  Created on: Oct 8, 2016
 *      Author: reginaldo
 */

#include "ChromosomeAckley.h"
#include "../functions/MathFunctions.h"
#include "../auxiliaries/Configuration.h"

ChromosomeAckley::ChromosomeAckley(double *genes) : Chromosome(genes) {

}

ChromosomeAckley::~ChromosomeAckley() {

}


void ChromosomeAckley::evaluate() {
	double fitness = MathFunctions::ackley(Chromosome::getGenes());
	Chromosome::setFitness(fitness);
}

ChromosomeAckley *ChromosomeAckley::createChromosome() {
	return new ChromosomeAckley(Chromosome::initializeGenesAtRandom(Ackley::NGENES, Ackley::INFIMUM, Ackley::MAXIMUM));
}
