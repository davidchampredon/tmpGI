//
//  individual.cpp
//  Gillespie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "individual.h"



void individual::create()
{
	_maxInfectiousStatus	= 99999;
	_infectiousStatus		= 0;	// susceptible by default
	_timeDiseaseAcquisition	= 9E9;	// very large time means not infected yet
	
	_infectorID				= 0;
	
	_GIbck					= -999.999;
	_GIfwd.resize(0);
}



void individual::incrementStatus()
{
	if(_infectiousStatus<_maxInfectiousStatus) _infectiousStatus++;
}


void individual::displayInfo()
{
	cout<<endl;
	cout<<"Infectious status = " << _infectiousStatus << endl;
	cout<<"Infection time = " << _timeDiseaseAcquisition << endl;
	cout<<"Infector ID = " << _infectorID << endl;
	cout<<"GIbck = " << _GIbck<< endl;
	cout<<"GIfwd:";
	if(_GIfwd.size()>0) displayVector(_GIfwd);
	cout<<endl;
}