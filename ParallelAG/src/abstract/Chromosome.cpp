/*
 * Chromosome.cpp
 *
 *  Created on: Oct 7, 2016
 *      Author: reginaldo
 */

#include "Chromosome.h"
#include "../auxiliaries/Rand.h"
#include "../auxiliaries/Configuration.h"
#include <iostream>
#include <iomanip>
using namespace std;

Chromosome::Chromosome(double *genes) {
	this->genes = genes;
	this->busy = false;
	this->fitness = 0;
}

Chromosome::~Chromosome() {
	if(genes)
		delete [] genes;
}

void Chromosome::setFitness(double fitness) {
	this->fitness = fitness;
}

double Chromosome::getFitness() const {
	return this->fitness;
}

void Chromosome::setGene(const int &pos, const double &value) {
	this->genes[pos] = value;
}

double Chromosome::getGene(const int &pos) const {
	return this->genes[pos];
}

double *Chromosome::getGenes() const {
	return this->genes;
}

void Chromosome::setBusy(const bool &busy) {
	this->busy = busy;
}

bool Chromosome::isBusy() const {
	return this->busy;
}

double *Chromosome::initializeGenesAtRandom(const int nGenes, const int infimum, const int maximum) {
	double *genes = new double[nGenes];
	double val = 0;
	Rand *r = new Rand();
	r->SetSeed();
	for (int i = 0; i < nGenes; i++) {
		val = r->Uniform();
		genes[i] = (1 - val) * infimum + val * maximum;
	}
	delete r;
	return genes;
}

int Chromosome::getNGenes(const int &functionID) {
	switch(functionID) {
	case Ackley::ID:
		return Ackley::NGENES;
	default:
		return 0;
	}
}

void Chromosome::print(const int &functionID) {
	int nGenes = getNGenes(functionID);
	cout << "g = {";
	for (int g = 0; g < nGenes; ++g) {
		cout << setw(WIDTH) << setprecision(PRECISION) << fixed << this->getGene(g);
		if(g != nGenes - 1)
			cout << ", ";
	}
	cout << "} | f(x) = " << this->getFitness() << endl;
}
