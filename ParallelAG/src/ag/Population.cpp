/*
 * Population.cpp
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#include "Population.h"
#include "../auxiliaries/Configuration.h"
#include "../chromosomePool/ChromosomeAckley.h"
#include <math.h>
#include <iostream>
#include <limits>
using namespace std;

Population::Population(int functionID) {
	this->functionID = functionID;
	this->fitnessMean = 0;
	this->fitnessStd = 0;

	this->chromos = new Chromosome*[NCHROMOSOMES];

	this->nGenes = Chromosome::getNGenes(functionID);
}

Population::~Population() {
	for (int i = 0; i < NCHROMOSOMES; ++i)
		if(chromos[i])
			delete chromos[i];
	if(chromos)
		delete [] chromos;
}

Chromosome *Population::chromosomeFactory(const int functionID) {
	switch(functionID) {
	case Ackley::ID:
		return ChromosomeAckley::createChromosome();
	default:
		return 0;
	}
}

/**
 * This function creates an empty chromosome.
 * 'Empty' here means genes equal to zero in all dimensions.
 */
Chromosome *Population::chromosomeEmptyFactory(const int functionID) {
	switch(functionID) {
	case Ackley::ID:
		return new ChromosomeAckley(new double[Ackley::NGENES]{0, 0});
	default:
		return 0;
	}
}

Population *Population::createPopulation(const int functionID, bool empty) {
	Population *p = new Population(functionID);
	if(!empty){
		for (int i = 0; i < NCHROMOSOMES; ++i)
			p->chromos[i] = chromosomeFactory(functionID);
	} else { // create empty chromosomes of the type functionID
		for (int i = 0; i < NCHROMOSOMES; ++i)
			p->chromos[i] = chromosomeEmptyFactory(functionID);
	}

	return p;
}

Chromosome *Population::getChromosome(const int &pos) const {
	return this->chromos[pos];
}

void Population::calculateFitnessMean() {
	double mean = 0;
	for (int i = 0; i < NCHROMOSOMES; i++)
		mean += this->chromos[i]->getFitness();
	this->fitnessMean = mean / NCHROMOSOMES;
}

void Population::calculateFitnessStd() {
	double std = 0;
	for (int i = 0; i < NCHROMOSOMES; i++)
		std += pow(this->chromos[i]->getFitness() - this->fitnessMean, 2);
	this->fitnessStd = sqrt(std / (NCHROMOSOMES - 1));
}

double Population::getFitnessMean() {
	return this->fitnessMean;
}

double Population::getFitnessStd() {
	return this->fitnessStd;
}

void Population::setChromosome(const int &pos, Chromosome *c) {
	for (int i = 0; i < this->nGenes; ++i)
		this->chromos[pos]->setGene(i, c->getGene(i));
	this->chromos[pos]->setFitness(c->getFitness());
	this->chromos[pos]->setBusy(c->isBusy());
}

int Population::getFittest() {
	int idx = 0;
	double bestFitness = numeric_limits<double>::max();

	for (int i = 0; i < NCHROMOSOMES; i++) {
		if(bestFitness > this->chromos[i]->getFitness()) {
			bestFitness = this->chromos[i]->getFitness();
			idx = i;
		}
	}
	return idx;
}

int Population::getFunctionID() {
	return this->functionID;
}

int Population::getNGenes() {
	return this->nGenes;
}

Chromosome **Population::getChromosomes() {
	return this->chromos;
}

void Population::print() {
	for (int i = 0; i < NCHROMOSOMES; ++i)
		this->chromos[i]->print(this->functionID);
}


