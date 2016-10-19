/*
 * MoonSlow.h
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#ifndef AG_MOONFAST_H_
#define AG_MOONFAST_H_

class Population;
class Chromosome;

struct DataF {
	Population *p;
	int t;
};

typedef struct DataF DataF;

void *innerCalculateFitness(void *);

class MoonFast {
private:
	Population *p;
	Population *selected;
	int functionID;

	static bool compare(Chromosome *, Chromosome *);

public:
	MoonFast(const int);
	virtual ~MoonFast();

	void initializePopulation();
	void calculateFitness(Population *);
	Population *parentSelection(Population *);
	Population *crossover(Population *);
	void mutation(Population *);
	Population *survivorSelection(Population *, Population *);
	Population *run();
};

#endif /* AG_MOONFAST_H_ */
