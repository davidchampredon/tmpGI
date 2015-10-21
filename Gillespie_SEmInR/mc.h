//
//  mc.h
//  Gillepsie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#ifndef __Gillepsie_SEmInR__mc__
#define __Gillepsie_SEmInR__mc__

#include <stdio.h>
#include "simulator.h"


void MC_run(simulator S, unsigned long iter_mc,
			double horizon, unsigned long initInfectious,
			int jobnum, string param_set,
			bool calc_WIW_Re,
			bool doExact, double timeStepTauLeap);

void MC_run__OLD(simulator S, unsigned long iter_mc,
			double horizon, unsigned long initInfectious,
			int jobnum, string param_set,
			bool calc_WIW_Re);


void MC_run_exact(simulator S, unsigned long iter_mc,
				  double horizon, unsigned long initInfectious,
				  int jobnum, string param_set,
				  bool calc_WIW_Re);


void MC_run_tauLeap(simulator S,
					unsigned long iter_mc,
					double horizon,
					double timeStep,
					unsigned long initInfectious,
					int jobnum,
					string param_set,
					bool calc_WIW_Re);



void save_all_outputs(simulator S,
					  string param_set, int jobnum,
					  unsigned long iter_mc,
					  unsigned long current_mc);

void test(simulator S);

#endif /* defined(__Gillepsie_SEmInR__mc__) */
