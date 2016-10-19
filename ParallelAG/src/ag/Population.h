/*
 * Population.h
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#ifndef AG_POPULATION_H_
#define AG_POPULATION_H_

#include "../abstract/Chromosome.h"

class Population {
private:
	Chromosome **chromos;

	int functionID;
	double fitnessMean;
	double fitnessStd;

	// Muste be static because of the call
	static Chromosome *chromosomeFactory(const int);
	static Chromosome *chromosomeEmptyFactory(const int);
	int nGenes;

public:
	Population(int functionID);
	virtual ~Population();

	static Population *createPopulation(const int, bool);

	Chromosome *getChromosome(const int &) const;

	void calculateFitnessMean();
	void calculateFitnessStd();
	double getFitnessMean();
	double getFitnessStd();

	void setChromosome(const int &, Chromosome *);
	int getFittest();
	int getFunctionID();
	int getNGenes();

	Chromosome **getChromosomes();

	void print();
};

#endif /* AG_POPULATION_H_ */
