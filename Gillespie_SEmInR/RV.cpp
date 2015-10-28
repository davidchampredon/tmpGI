/*
 *  RV.cpp
 *
 *
 *  Created by David Champredon  on 13-04-19.
 *  Copyright 2013-14. All rights reserved.
 *
 */

#include "RV.h"
#include "dcTools.h"
#include "globalVar.h"

using namespace std;


void force_seed_reset()
{
	_RANDOM_GENERATOR.seed(_RANDOM_SEED);
}

void force_seed_reset(unsigned int manual_seed)
{
	_RANDOM_GENERATOR.seed(manual_seed);
}



/* ****************************************************
 CONTINUOUS SUPPORT
 ****************************************************/


double uniform01()
{
	uniform_real_distribution<double> dist(0,1);
	return dist(_RANDOM_GENERATOR);
}


double expo(double lambda)
{
	exponential_distribution<double> dist(lambda);
	return dist(_RANDOM_GENERATOR);
	
}


double normal(double mean, double sd)
{
	normal_distribution<double> dist(mean, sd);
	return dist(_RANDOM_GENERATOR);
}

double gamma(double shape, double scale)
{
	gamma_distribution<double> dist(shape,scale);
	return dist(_RANDOM_GENERATOR);
}

double beta(double a, double b)
{
	gamma_distribution<double> X(a,1.0);
	gamma_distribution<double> Y(b,1.0);
	
	double xx = X(_RANDOM_GENERATOR);
	double yy = Y(_RANDOM_GENERATOR);
	
	return xx/(xx+yy);
}




/* ****************************************************
			DISCRETE SUPPORT
 ****************************************************/


int	uniformInt(int nmin, int nmax)
{
	uniform_int_distribution<> dist(nmin,nmax);
	return dist(_RANDOM_GENERATOR);
}

int geometric(double p)
{
	geometric_distribution<int> dist(p);
	return dist(_RANDOM_GENERATOR);
}

int	binom(double p, int N)
{
	stopif( (p<0 || p>1),
		   "probability ("+to_string(p)+") not in [0;1]");

	binomial_distribution<int> dist(N,p);
	return dist(_RANDOM_GENERATOR);
}


unsigned long poisson(double expectedValue)
{
	poisson_distribution<unsigned long> dist(expectedValue);
	return dist(_RANDOM_GENERATOR);
}



vector<unsigned int> multinomial(unsigned int N, vector<double> proba)
{
	// Retrieve the dimension that defines the random variables vector
	unsigned long dim = proba.size();
	
	// Healthy checks
	if (dim==0 || N==0)
	{
		cout << endl << "ERROR [multinomial]: can't generate with dimension=0"<<endl;
		if (dim==0) cout << "size of proba vector = 0 !"<<endl;
		if (N==0) cout << "N = 0 !"<<endl;
		exit(1);
	}
	double s=0;
	for (int i=0; i<proba.size(); i++) {s+=proba[i];}
	if (fabs(s-1)>0.0001)
	{
		cout << endl << "ERROR [multinomial]: probabilities do not add-up to 1"<<endl;
		displayVector(proba);
		exit(1);
	}
	
	// Partition interval [0;1] with respect to probabilities
	vector<double> a(dim+1);
	a[0]=0;
	for (int i=1; i<dim+1; i++) 
	{
		a[i] = a[i-1]+proba[i-1];
	}
	
	// initiate the result vector 

	vector<unsigned int> res(dim,0); // all values at 0
	
	// Draw the N trials
	
	for (int k=0; k<N; k++) 
	{
		double u = uniform01();
		
		// Find the interval where "u" is.
		int m;
		for (m=1; a[m]<u; m++) {}
		
		// "u" is between a[m-1] and a[m]
		// so the (m-1)th value of the stochastic vector is increased
		
		res[m-1]++;
	}
	
	return res;
}



int	probaHistogramInt(vector<int> values, vector<double> probas)
{
	
	// --- Healthy checks ---
	if (values.size() != probas.size())
	{
		cout<<endl<<"ERROR [probaHistogramInt]: size of 'values' ("<< values.size()
		<<") and 'probas' ("<< values.size()<<") not the same!"<<endl;
		exit(1);
	}
	
	if (sumElements(probas)!=1)
	{
		cout<<endl<<"ERROR [probaHistogramInt]: Probabilities do not add up to 1!"<<endl;
		exit(1);
	}
	
	// ----------------------
	
	double s=probas[0];
	
	int k = 0;
	int kk = 0;
	
	bool found = false;
	
	double u = uniform01();
	
	while(k<probas.size() && !found)
	{
		if (u<s) {
			kk=k;
			found = true;
		}
		
		k++;
		s += probas[k];
	}
	
	return values[kk];
}
	



double uniform(int seed)
{
	
	// Generator: L'Ecuyer (good quality)
	// 'seed' initializes the generator
	// source: Numerical Recipies in C (2nd edition)
	
	seed = -seed; // initialization must be negative
	
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	
	if (seed <= 0)
	{
		if (-(seed)<1) seed=1;
		else seed=-(seed);
		idum2=(seed);
		for (j=NTAB+7;j>=0;j--)
		{
			k=(seed)/IQ1;
			seed=IA1*(seed-k*IQ1)-k*IR1;
			if ( seed < 0 ) seed+=IM1;
			if (j<NTAB) iv[j]= seed;
		}
		iy=iv[0];
	}
	k=(seed)/IQ1;
	seed=IA1*(seed-k*IQ1)-k*IR1;
	if (seed<0) seed += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2<0) idum2+=IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j]= seed;
	if (iy<1) iy+=IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double expo_old(double lambda,int seed)
{
	double val;
	
	stopif(!(lambda>0),"Exp rand. var. parameter not positive!");
	
	val=(-log(uniform(seed))/lambda);
	
	return val;
}


double normal_old(double mean, double stddev)
{
	// BOX-MULLER algorithm
	
	double u1 = uniform(rand());
	double u2 = uniform(rand());
	
	return mean + stddev*(sqrt(-2*log(u1))*cos(2*PI*u2) );
}
int	binom_old(int seed, double p, int N)
{
	int b=0;
	
	for (int i=0; i<N; i++)
	{
		double u = uniform(seed+i);
		if (u<p) b++;
	}
	
	return b;
}

int poisson_old(int seed, double expectedValue)
{
	int n = 0; //counter of iteration
	double limit;
	double x;  //pseudo random number
	limit = exp(-expectedValue);
	x = uniform(seed);
	while (x > limit)
	{
		n++;
		x *= uniform(seed+rand());
	}
	return n;
}
