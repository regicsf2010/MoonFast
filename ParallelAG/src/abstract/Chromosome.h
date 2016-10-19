/*
 * Chromosome.h
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

class Chromosome {
private:
	double fitness;
	double *genes;
	bool busy;

public:
	Chromosome(double *);
	virtual ~Chromosome();

	void setFitness(double);
	double getFitness() const;

	void setGene(const int &, const double &);
	double getGene(const int &) const;
	double *getGenes() const;

	void setBusy(const bool &);
	bool isBusy() const;

	static double *initializeGenesAtRandom(const int, const int, const int);
	static int getNGenes(const int &);

	void print(const int &);

	virtual void evaluate() = 0;
};

#endif /* CHROMOSOME_H_ */
