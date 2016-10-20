/*
 * MoonSlow.cpp
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#include <iostream>
using namespace std;

#include "MoonFast.h"
#include "Population.h"
#include "../auxiliaries/Rand.h"

#include <algorithm>

MoonFast::MoonFast(const int functionID) {
	this->functionID = functionID;
	p = 0;
	selected = 0;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&cond, NULL);
}

MoonFast::~MoonFast() {
	if(p)
		delete p;
	if(selected)
		delete selected;
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&cond);
}

void MoonFast::initializePopulation() {
	this->p = Population::createPopulation(this->functionID, false);
}

void *innerCalculateFitness(void *data) {
	DataF *d = (DataF *) data;
	int partition = NCHROMOSOMES / NTHREADS;
	for (int i = partition * d->t; i < partition * d->t + partition; ++i)
		d->p->getChromosome(i)->evaluate();
	pthread_exit(NULL);
}

void MoonFast::calculateFitness(Population *pop) {
	DataF data;
	data.p = pop;
	for (int t = 0; t < NTHREADS; t++) {
		data.t = t;
		pthread_create(&threads[t], &attr, innerCalculateFitness, (void *)&data);
	}

	for(int t = 0; t < NTHREADS; t++)
		pthread_join(threads[t], NULL);
}

void *innerCalculateFitnessMeanAndStd(void *data) {
	DataF *d = (DataF *) data;
	if(d->t == 0) {
		pthread_cond_wait(d->cond, d->mutex);
		d->p->calculateFitnessStd();
	} else {
		d->p->calculateFitnessMean();
		pthread_cond_signal(d->cond);
	}
	pthread_exit(NULL);
}

void MoonFast::calculateFitnessMeanAndStd(Population *pop) {
	DataF data;
	data.p = pop;
	data.mutex = &mutex;
	data.cond = &cond;

	for (int t = 0; t < 2; ++t) {
		data.t = t;
		pthread_create(&meanStd_t[t], &attr, innerCalculateFitnessMeanAndStd, (void *)&data);
	}

	for(int t = 0; t < 2; t++)
		pthread_join(meanStd_t[t], NULL);
}


Population *MoonFast::parentSelection(Population *pop) {
	Population *selected = Population::createPopulation(pop->getFunctionID(), true);

	int idx1 = -1, idx2 = -1;
	Rand *r = new Rand();
	r->SetSeed();

	for (int i = 0; i < NCHROMOSOMES; i++) {
		idx1 = r->RandInt(NCHROMOSOMES);

		for (int j = 1; j < RANK; j++) { // j starts at 1
			idx2 = r->RandInt(NCHROMOSOMES);
			if(pop->getChromosome(idx1)->getFitness() >  pop->getChromosome(idx2)->getFitness())
				idx1 = idx2;
		}

		selected->setChromosome(i, pop->getChromosome(idx1));
	}

	delete r;
	return selected;
}

Population *MoonFast::crossover(Population *parents) {
	Population *offspring = Population::createPopulation(parents->getFunctionID(), true);
	int nGenes = parents->getNGenes();
	Rand *r = new Rand();
	r->SetSeed();

	for (int i = 0; i < NCHROMOSOMES; i += 2) {
		for (int j = 0; j < nGenes; j++) {
			double a1 = parents->getChromosome(i)->getGene(j);
			double a2 = parents->getChromosome(i + 1)->getGene(j);

			if(r->Uniform() < CROSSOVERRATE) {
				offspring->getChromosome(i)->setGene(j, a1 + (a2 - a1) * r->Uniform());
				offspring->getChromosome(i + 1)->setGene(j, a2 + (a1 - a2) * r->Uniform());
			} else {
				offspring->getChromosome(i)->setGene(j, a1);
				offspring->getChromosome(i + 1)->setGene(j, a2);
			}

		}
	}

	delete r;
	delete parents;
	return offspring;
}

void MoonFast::mutation(Population *offspring) {
	Rand *r = new Rand();
	r->SetSeed();
	int nGenes = offspring->getNGenes();
	double val = 0, newVal = 0;

	for (int i = 0; i < NCHROMOSOMES; i++) {
		for (int j = 0; j < nGenes; j++) {

			if(r->Uniform() < MUTATIONRATE) {
				val = offspring->getChromosome(i)->getGene(j);
				newVal = val + SD * r->Normal();
				offspring->getChromosome(i)->setGene(j, newVal);
			}

		}
	}
	delete r;
}

Population *MoonFast::survivorSelection(Population *p1, Population *p2) {
	Population *survivors = Population::createPopulation(p1->getFunctionID(), true);

	sort(p1->getChromosomes(), p1->getChromosomes() + NCHROMOSOMES, compare);
	sort(p2->getChromosomes(), p2->getChromosomes() + NCHROMOSOMES, compare);

	for (int i = 0; i < NCHROMOSOMES; i++) {
		if(p1->getChromosome(i)->getFitness() < p2->getChromosome(i)->getFitness())
			survivors->setChromosome(i, p1->getChromosome(i));
		else
			survivors->setChromosome(i, p2->getChromosome(i));
	}
	delete p1;
	delete p2;
	return survivors;
}

bool MoonFast::compare(Chromosome *x, Chromosome *y){
	return x->getFitness() < y->getFitness();
}

Population *MoonFast::run() {
	this->initializePopulation();
	this->calculateFitness(p);

	for (int i = 0; i < NGENERATIONS; ++i) {
		selected = this->parentSelection(p);
		selected = this->crossover(selected);
		this->mutation(selected);
		this->calculateFitness(selected);
		p = this->survivorSelection(p, selected);
	}
	return p;
}


