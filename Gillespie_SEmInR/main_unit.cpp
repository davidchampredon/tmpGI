//
//  main_unit.cpp
//  Gillespie_SEmInR
//
//  Created by David CHAMPREDON on 2015-03-04.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//


#include <stdlib.h>
#include <iostream>
#include "simulator.h"
#include "individual.h"
#include "mc.h"

int main(int argc, const char * argv[]) {
	
	system("pwd");
	system("date");
	
	// For performance monitoring
	// - do not delete -
	timeval tim;
	gettimeofday(&tim, NULL);
	double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	// ------------------------------------------
	
	
	// Read the job number
	int jobnum = atoi(argv[1]);
	string fileparam = argv[2];
	
	// Read main simulation parameters from file
	

	
	double horizon			= getParameterFromFile("horizon", fileparam);
	unsigned long popSize	= getParameterFromFile("popSize", fileparam);
	double R0				= getParameterFromFile("R0", fileparam);;
	double latent_mean		= getParameterFromFile("latent_mean", fileparam);;
	double infectious_mean	= getParameterFromFile("infectious_mean", fileparam);;
	int nE					= getParameterFromFile("nE", fileparam);
	int nI					= getParameterFromFile("nI", fileparam);
	unsigned long mc_iter	= getParameterFromFile("mc_iter", fileparam);
	int njobs				= getParameterFromFile("njobs", fileparam);
	
	unsigned long initInfectious	= getParameterFromFile("init_I1", fileparam);
	bool calc_WIW_Re				= getParameterFromFile("calc_WIW_Re", fileparam);
	bool doExact					= getParameterFromFile("doExact", fileparam);
	double timeStepTauLeap			= getParameterFromFile("timeStepTauLeap", fileparam);
	
	int mc_job = int(mc_iter/njobs);
	
	
	// Derive other variables
	double sigma0	= 1/latent_mean;
	double gamma0	= 1/infectious_mean;
	double beta		= R0*gamma0;
	
	
	vector<double> sigma(nE);
	vector<double> gamma(nI);
	for (int i=0; i<nE; i++) sigma[i]=sigma0*nE;
	for (int i=0; i<nI; i++) gamma[i]=gamma0*nI;
	
	// Simulation
	
	simulator SIM(beta, sigma, gamma, popSize, nE, nI);
	
	MC_run(SIM, mc_job, horizon,initInfectious,
		   jobnum,fileparam,calc_WIW_Re,
		   doExact,timeStepTauLeap);
	
	cout<<endl<<"--- Job #"<<jobnum<<" finished!"<<endl;
	
	
	// --------------------------------------------------------------
	// COMPUTER TIME MONITORING - do not delete!
	
	gettimeofday(&tim, NULL);
	double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	int minutes = (int)((t2-t1)/60.0);
	double sec = (t2-t1)-minutes*60.0;
	cout << endl << " - - - Computational time : ";
	cout << minutes<<" min "<<sec<<" sec" << endl;
	
	// --------------------------------------------------------------
	
	
	
	return 0;
}
