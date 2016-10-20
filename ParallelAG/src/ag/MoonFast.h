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

#include "../auxiliaries/Configuration.h"
#include <pthread.h>

struct DataF {
	Population *p;
	int t;
	pthread_mutex_t *mutex;
	pthread_cond_t *cond;
};

typedef struct DataF DataF;

void *innerCalculateFitness(void *);
void *innerCalculateFitnessMeanAndStd(void *);

class MoonFast {
private:
	Population *p;
	Population *selected;
	int functionID;

	static bool compare(Chromosome *, Chromosome *);

	pthread_t threads[NTHREADS];
	pthread_t meanStd_t[2];
	pthread_attr_t attr;
	pthread_mutex_t mutex;
	pthread_cond_t cond;

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
	void calculateFitnessMeanAndStd(Population *);
};

#endif /* AG_MOONFAST_H_ */
