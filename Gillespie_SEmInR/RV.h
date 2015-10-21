/*
 *  RV.h
 *  Gillespie
 *
 *  Created by David Champredon  on 13-04-19.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

// THINK OF USING "Boost Random Number Library"

#ifndef RV_H
#define RV_H

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <random>

//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

//#include "dcTools.h"

using namespace std;

#define PI 3.141592654

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)



// === Random number generator seed ====

void force_seed_reset();
void force_seed_reset(unsigned int);


// === Continuous support ====

double	uniform(int seed);				// Uniform in [0;1]
double	uniform01();

double	expo(double lambda);
double	normal(double mean, double stddev);

double	gamma(double shape, double scale);
double	beta(double a, double b);


// === Discrete support ====

int		uniformInt(int min, int max); // Uniform INTEGER in [min;max]
int		geometric(double p);
int		binom(double p, int N);	// Binomial(p,N) N trials with probability of success p
int		binom_old(int seed, double p, int N);	// Binomial(p,N) N trials with probability of success p
unsigned long	poisson(double expectedValue);
int		poisson_old(int seed, double expectedValue);

// Explicitly define values andassociated probabilities
// (histogram-like)
int		probaHistogramInt(vector<int> values, vector<double> probas);


vector<unsigned int> multinomial(unsigned int N, 
								 vector<double> proba); // Multinomial (N,p) ; returns vector of size proba.size()


// === GSL wrappers ===


//gsl_rng * GSL_generator(unsigned int seed);
//
//// Multinomial (N,p) ; returns vector of size proba.size()
//vector<unsigned int> multinomial_gsl(gsl_rng * r,unsigned int N, vector<double> proba);



#endif