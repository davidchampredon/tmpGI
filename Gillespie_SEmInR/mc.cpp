//
//  mc.cpp
//  Gillepsie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "mc.h"
#include "globalVar.h"


void MC_run(simulator S,
			unsigned long iter_mc,
			double horizon,
			unsigned long initInfectious,
			int jobnum,
			string param_set,
			bool calc_WIW_Re,
			bool doExact,
			double timeStepTauLeap)
{
	/// Run several iterations of the simulations (Monte Carlo)
	/// Can use either the exact Gillespie algo, or tau leap approximation
	
	if (doExact) MC_run_exact(S, iter_mc, horizon, initInfectious, jobnum, param_set, calc_WIW_Re);
	if (!doExact)MC_run_tauLeap(S, iter_mc, horizon, timeStepTauLeap, initInfectious, jobnum, param_set, calc_WIW_Re);
}





void MC_run_exact(simulator S, unsigned long iter_mc,
				  double horizon, unsigned long initInfectious,
				  int jobnum, string param_set,
				  bool calc_WIW_Re)
{
	// Forces the seed to be different for each job
	// (each job is executed independently)
	force_seed_reset(jobnum*7);
	
	unsigned long cnt=0;
	
	// Monte-Carlo loop
	for(unsigned long i=0; i<iter_mc; i++)
	{
		cout<<endl<<"MC "<<i+1<<"/"<<iter_mc<< " (exact)"<<endl;
		
		S.run_exact(horizon,initInfectious,calc_WIW_Re);
		save_all_outputs(S, param_set, jobnum, iter_mc, i);
		S.displayInfo();
		// counts the number of fizzles
		if (S.get_cumIncidence()[S.get_cumIncidence().size()-1] < S.get_popSize()/20) cnt++;
	}
	cout<<endl<<"Fizzle proportion = "<<(double)(cnt)/iter_mc<<endl;
}



void MC_run_tauLeap(simulator S, unsigned long iter_mc,
					double horizon, double timeStep,
					unsigned long initInfectious,
					int jobnum, string param_set,
					bool calc_WIW_Re)
{
	// Forces the seed to be different for each job
	// (each job is executed independently)
	force_seed_reset(jobnum*7);
	
	unsigned long cnt=0;
	
	// Monte-Carlo loop
	for(unsigned long i=0; i<iter_mc; i++)
	{
		cout<<endl<<"MC "<<i+1<<"/"<<iter_mc<<" (tau leap "<<timeStep<<")"<<endl;
		
		S.run_tauLeap(horizon, timeStep, initInfectious, calc_WIW_Re);
		save_all_outputs(S, param_set, jobnum, iter_mc, i);
		S.displayInfo();
		// counts the number of fizzles
		if (S.get_cumIncidence()[S.get_cumIncidence().size()-1] < S.get_popSize()/20) cnt++;
	}
	cout<<endl<<"Fizzle proportion = "<<(double)(cnt)/iter_mc<<endl;
}





void save_all_outputs(simulator S,
					  string param_set, int jobnum,
					  unsigned long iter_mc,
					  unsigned long current_mc)
{
	/// Save all relevant outputs to files
	
	string f_prev	= param_set + "__prev";
	string f_cumInc	= param_set + "__cumInc";
	string f_nS		= param_set + "__nS";
	string f_nR		= param_set + "__nR";
	string f_GIbck	= param_set + "__GIbck";
	string f_GIfwd	= param_set + "__GIfwd";
	string f_Reff	= param_set + "__Reff";
	
	string post =  to_string((jobnum-1)*iter_mc+current_mc+1) + ".out";
	
	// Prevalence
	string tmp_prev = _DIR_OUT +f_prev + post;
	S.save_prevalence(tmp_prev);
	
	// Cumulative Incidence
	string tmp_cumInc = _DIR_OUT +f_cumInc + post;
	S.save_cumIncidence(tmp_cumInc);
	
	// Susceptible
	string tmp_nS = _DIR_OUT +f_nS + post;
	S.save_nS(tmp_nS);
	
	// Recovered
	string tmp_nR = _DIR_OUT +f_nR + post;
	S.save_nR(tmp_nR);
	
	// Generation interval
	string tmp_gibck = _DIR_OUT +f_GIbck + post;
	S.save_GIbck(tmp_gibck);
	
	string tmp_gifwd = _DIR_OUT +f_GIfwd + post;
	S.save_GIfwd(tmp_gifwd);
	
	// Realized Reff
	string tmp_Reff = _DIR_OUT +f_Reff + post;
	S.save_Reff(tmp_Reff);
	
}




void test(simulator S)
{
	for (int i=0;i<10;i++) {
		cout<<i<<"->"<<uniform01()<<endl;
		
	}
}


void MC_run__OLD(simulator S, unsigned long iter_mc,
				 double horizon, unsigned long initInfectious,
				 int jobnum, string param_set,
				 bool calc_WIW_Re)
{
	// Forces the seed to be different for each job
	// (each job is executed independently)
	force_seed_reset(jobnum*7);
	
	unsigned long cnt=0;
	
	// Monte-Carlo loop
	for(unsigned long i=0; i<iter_mc; i++)
	{
		cout<<endl<<"MC "<<i+1<<"/"<<iter_mc<<endl;
		
		S.run_exact(horizon,initInfectious,calc_WIW_Re);
		save_all_outputs(S, param_set, jobnum, iter_mc, i);
		S.displayInfo();
		// counts the number of fizzles
		if (S.get_cumIncidence()[S.get_cumIncidence().size()-1] < S.get_popSize()/20) cnt++;
	}
	cout<<endl<<"Fizzle proportion = "<<(double)(cnt)/iter_mc<<endl;
}
