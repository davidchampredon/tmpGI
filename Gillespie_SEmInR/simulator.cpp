//
//  simulator.cpp
//  Gillespie_SEmInR
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "simulator.h"


simulator::simulator(double beta,
					 vector<double>sigma,
					 vector<double>gamma,
					 unsigned long popSize,
					 int nE, int nI)
{
	/// Constructs a simulator object
	
	_beta = beta;
	_sigma = sigma;
	_gamma = gamma;
	
	_popSize = popSize;
	_indiv.resize(popSize);
	
	_nE = nE;
	_nI = nI;
	
	for(unsigned long i=0; i<popSize; i++)	{
		_indiv[i].create();
		_indiv[i].set_ID(i);
		_indiv[i].set_maxInfectiousStatus(nE+nI+1);
	}
}



void simulator::displayInfo()
{
	coutline(40);
	cout<<" === SIMULATOR INFO ==="<<endl;
	
	cout << "Population size: "<<_indiv.size()<<endl;
	
	cout << "# E compartments:" << _nE <<endl;
	cout << "# I compartments:" << _nI <<endl;
	
	cout << "sigma rates:"; displayVector(_sigma);
	cout << "gamma rates:"; displayVector(_gamma);
	double R0 = _beta/_gamma[0]*_nI;
	cout << "R0 = "<<R0<<endl;
	
	cout << "Number of events: "<<_time.size()<<endl;
	cout << "Last date: "<<_time[_time.size()-1]<<endl;
	cout << "Cumul incidence: "<<_cumIncidence[_cumIncidence.size()-1]<<endl;
	coutline(40);
}




double simulator::eventRate_infection(){
	return at_least_one_S_and_I()?_beta*_count_I*_count_S/_popSize:0.0;
}


double simulator::eventRate_latencyProgress(unsigned int i){
	
	/// Returns the rate to progress from the ith latent stage
	/// (stages from '0' to '_nE-1')
	
	stopif(i>=_nE, "index out of bounds");
	return _sigma[i]*census_status(i+1);
}


double simulator::eventRate_infectiousProgress(unsigned int i){
	
	/// Returns the rate to progress from the ith infectious stage
	/// (stages from '0' to '_nI-1')
	
	stopif(i>=_nI, "index out of bounds");
	return _gamma[i]*census_status(_nE+1+i);
}



double simulator::eventRate_all_latencyProgress(){
	double A_E = 0.0;
	for(int i=1; i<=_nE; i++) A_E += _sigma[i-1]*census_status(i);
	return A_E;
}


double simulator::eventRate_all_infectiousProgress(){
	double A_I = 0.0;
	for(int i=_nE+1; i<=_nE+_nI; i++)
		A_I += _gamma[i-1-_nE]*census_status(i);
	return A_I;
}



vector<unsigned long> simulator::drawNumberEvents_tauLeap(double timestep)
{
	/// Draws the number of events during a tau-leap
	/// for every event
	
	vector<unsigned long> res;
	
	// New infection event
	res.push_back(poisson(eventRate_infection()*timestep));
	
	// Migrations through the E[k]
	for(int i=0; i<_nE; i++) res.push_back(poisson(eventRate_latencyProgress(i)*timestep));
	
	// Migrations through the I[k]
	for(int i=0; i<_nI; i++) res.push_back(poisson(eventRate_infectiousProgress(i)*timestep));

	return res;
}






double simulator::drawNextEventInterval()
{
	/// Draw the next event time
	
	// Infection
	double A_S = _beta*_count_I*_count_S/_popSize;
	
	//DEBUG
	if(_count_I !=_prevalence[_prevalence.size()-1]){
	cout<< "_count_I:"<<_count_I << "; prevalence[t]:"<<_prevalence[_prevalence.size()-1]<<endl;
	}
	
	// Latency progressions
	double A_E = 0.0;
	for(int i=1; i<=_nE; i++)
		A_E += _sigma[i-1]*census_status(i);
	
	// Infectiousness progressions
	// and ultimately recovery
	double A_I = 0.0;
	for(int i=_nE+1; i<=_nE+_nI; i++)
		A_I += _gamma[i-1-_nE]*census_status(i);
	
	// Sum of all rates
	double A = A_S+A_E+A_I;
	
	double u = uniform01();
	
	return -log(u)/A;
}



unsigned int simulator::drawEventType()
{
	/// Determine which event type occurs
	/// (after an event time has been drawn)
	
	
	vector<double> x;
	
	// sufficiently small such that probability to fall
	// b/w x and x+tiny is pratically 0 but large enough
	// that a 'double comparison' recognizes it
	// (used when no indiv in a given status)
	double tiny = 1E-10;
	
	// Anchor at 0
	x.push_back(0.0);
	
	// New infection event
	double tmp = at_least_one_S_and_I()?_beta*_count_I*_count_S/_popSize:tiny;
	x.push_back(tmp);
	
	// Migrations through the E[k]
	for(int i=0; i<_nE; i++)
	{
		double tmp = (census_status(i+1)>0)?_sigma[i]*census_status(i+1):tiny;
		x.push_back(tmp);
	}
	
	// Migrations through the I[k]
	for(int i=0; i<_nI; i++)
	{
		double tmp = (census_status(_nE+1+i)>0)?_gamma[i]*census_status(_nE+1+i):tiny;
		x.push_back(tmp);
	}
	
	// draw a uniform on sum of all rates and
	// see where it lands -> that's the drawn event type
	double u = uniform01()*sumElements(x);
	
	//	cout<<"drawEventType: "<<u<<" :: " <<sumElements(x)<<endl;
	
	
	// vector of cumulative sum of x elements
	vector<double> cumx;
	for(int i=0; i<x.size(); i++)
	{
		// initialize sum at 0
		cumx.push_back(0.0);
		// Calculate cumulative value
		for(int j=0; j<=i; j++) cumx[i]+=x[j];
	}
	
	unsigned int idx = 0;
	bool found = false;
	for(int i=0; i<cumx.size(); i++)
	{
		if(cumx[i]<u)
		{
			idx = i;
			found = true;
		}
	}
	stopif(!found, "can't find index");
	
	// Make sure there is an individual
	// in the drawn event (they may be very close, see 'tiny')
	if (census_status(idx)==0)
	{
		unsigned int idx2 = idx;
		while (census_status(idx2)==0) idx2=idx2-1;
		idx=idx2;
	}
	
	return idx;
}



void simulator::clean_start()
{
	/// Makes sure the simulation starts clean
	
	_prevalence.clear();
	_time.clear();
	_cumIncidence.clear();
	_nS.clear();
	_nR.clear();
	_WIW_times.clear();
	_WIW.clear();
	_Reff.clear();

	for(int i=0;i<_popSize; i++) _indiv[i].create();
}


void simulator::initialize(unsigned long initInfectious)
{
	/// Initialize before running the simulation
	
	clean_start();
	
	// Set initial infectious individuals
	
	for(int i=0; i<initInfectious; i++)
	{
		// Warning: initial infectious individuals
		// are put in I[1] (_not_ exposed/latent stage)
		_indiv[i].set_infectiousStatus(_nE+1);
		_indiv[i].set_timeDiseaseAcquisition(0.0);
	}
	
	// intialize counts
	_count_I = initInfectious;
	_count_S = _popSize - initInfectious;
	_count_R = 0;
	
	// Initialization
	
	double t=0.0;
	_time.push_back(t);
	_cumIncidence.push_back(initInfectious);
	_prevalence.push_back(initInfectious);
	_nS.push_back(_popSize - initInfectious);
	_nR.push_back(0);

}


void simulator::actionOnEvent(unsigned int eventType, double time_event){
	
	/// Performs updates on individuals given an event type
	
	// Event is NOT an new infection
	if(eventType>0){
		// retrieve all the relevant IDs
		vector<unsigned long> x = census_ID(eventType);
		
		if(x.size()>0){
			// Pick one randomly
			unsigned long ID_selected = extractElementRandom(x);
			// Move this individual to the next stage
			_indiv[ID_selected].incrementStatus();
		}
	}
	
	// Event is an new infection
	if(eventType==0 && _count_S>0){
		// retrieve all the susceptible IDs
		vector<unsigned long> x = census_ID(0);
		
		// Pick new infectee randomly
		unsigned long ID_infectee = extractElementRandom(x);
		
		// retrieve all the infectious IDs
		vector<unsigned long> x_I = census_ID_I();
		
		// Pick new infector randomly
		unsigned long ID_infector = extractElementRandom(x_I);
		
		// Now that we know who infected whom,
		// update all relevant quantities
		
		// > infectious status and infector ID
		set_infectiousStatus(ID_infectee, 1);
		set_infectorID(ID_infectee, ID_infector);
		
		// > timing of infections from both view points
		set_timeDiseaseAcquisition(ID_infectee, time_event);
		set_timeDiseaseTransmit(ID_infector, time_event);
		
		// > generation intervals
		double gi = time_event - get_timeDiseaseAcquisition(ID_infector);
		set_GIbck(ID_infectee, gi);
		set_GIfwd(ID_infector, gi);
	}
	
	
	// -- Keep track of S and I counts --
	// transmission:
	if(eventType==0 && _count_S>0) _count_S--;
	
	// from last E to I[1]:
	if(eventType==_nE) _count_I++;
	
	// from last I to R
	if(eventType==(_nE+_nI) && _count_I>0){
		_count_I--;
		_count_R++;
	}
}



void simulator::run(double horizon,
					unsigned long initInfectious,
					bool calc_WIW_Re,
					bool doExact,
					double timeStep)
{
	/// Run a simulation
	/// either 'exact' Gillespie, or approximation with 'tau leap'
	
	if(doExact) run_exact(horizon, initInfectious, calc_WIW_Re);
	if(!doExact) run_tauLeap(horizon, timeStep, initInfectious, calc_WIW_Re);
}


void simulator::run_exact(double horizon,
						  unsigned long initInfectious,
						  bool calc_WIW_Re)
{
	/// Run the epidemic simulation with the exact Gillespie algorithm
	
	initialize(initInfectious);
	double t = 0.0;
	
	// Simulation
	
	while ( t<horizon && at_least_one_S_and_E_or_I() )
	{
		// Draw the next event time based on the sum of specified rates
		double dt = drawNextEventInterval();
		double time_event = t+dt;
		
		// Draw the type of event based on each specified rates
		unsigned int et	= drawEventType();
		
		// Update individuals status based on event
		actionOnEvent(et, time_event);
		
		// Update WIW matrix
		// (only at rounded times)
		if (calc_WIW_Re && t-dt<int(t) && int(t)<t)	calc_WIW(t);
		
		// update times
		_time.push_back(time_event);
		t += dt;

		// calculate cumulative incidence
		double inc = (et==0)?1:0;
		inc += _cumIncidence[_cumIncidence.size()-1];
		_cumIncidence.push_back(inc);
		
		
		// update prevalence time series
		_prevalence.push_back(_count_I);
		
		// update susceptible counts time series
		_nS.push_back(_count_S);
		_nR.push_back(_count_R);
		
	} // end-while-loop
	
	// Calculate time series of case effective reproductive number _Reff
	// (_WIW must be calculated until horizon)
	if(calc_WIW_Re) {
		calc_Reff();
		// Integrity check
		if (_WIW.size()>0){
			unsigned long n_wiw = _WIW[_WIW.size()-1].countNonZeroElements();
			unsigned long n_cuminc =_cumIncidence[_cumIncidence.size()-1];
			stopif(n_wiw>n_cuminc, "cumulative incidences (_WIW vs _cumIncidence) not consistent!");
		}
	}
}



void simulator::run_tauLeap(double horizon,
							double timestepSize,
							unsigned long initInfectious,
							bool calc_WIW_Re)
{
	/// Run the epidemic simulation with the
	/// tau-leap Poisson approximation of Gillespie algorithm
	
	initialize(initInfectious);
	double t = 0.0;
	
	while ( t<horizon && at_least_one_S_and_E_or_I() )
	{
		double					time_event	= t+timestepSize;
		vector<unsigned long>	nEvents		= drawNumberEvents_tauLeap(timestepSize);
		
		for(int event=0; event<nEvents.size(); event++){	// <-- loop on all event types
			for (int k=0; k<nEvents[event]; k++) {			// <-- do action for as many number of events was drawn for this event type
				actionOnEvent(event, time_event);
			}
		}
		
		// update times
		_time.push_back(time_event);
		t += timestepSize;
		
		// update cumulative incidence
		double previous_cumInc = _cumIncidence[_cumIncidence.size()-1];
		_cumIncidence.push_back(previous_cumInc + nEvents[0]);

		// update prevalence time series
		_prevalence.push_back(_count_I);
	
		// update susceptible counts time series
		_nS.push_back(_count_S);
		_nR.push_back(_count_R);
		
	} // end while
	
	// Calculate time series of case effective reproductive number _Reff
	// (_WIW must be calculated until horizon)
	if(calc_WIW_Re) {
		calc_Reff();
		// Integrity check
		if (_WIW.size()>0){
			unsigned long n_wiw = _WIW[_WIW.size()-1].countNonZeroElements();
			unsigned long n_cuminc =_cumIncidence[_cumIncidence.size()-1];
			stopif(n_wiw>n_cuminc, "cumulative incidences (_WIW vs _cumIncidence) not consistent!");
		}
	}
}



void simulator::run__OLD(double horizon,
					unsigned long initInfectious,
					bool calc_WIW_Re)
{
	/// Run the epidemic simulation with the exact Gillespie algorithm
	
	clean_start();
	
	// Set initial infectious individuals
	
	for(int i=0; i<initInfectious; i++)
	{
		// Warning: initial infectious individuals
		// are put in I[1] (_not_ exposed/latent stage)
		_indiv[i].set_infectiousStatus(_nE+1);
		_indiv[i].set_timeDiseaseAcquisition(0.0);
	}
	
	// intialize counts
	_count_I = initInfectious;
	_count_S = _popSize - initInfectious;
	_count_R = 0;
	
	// Initialization
	
	double t=0.0;
	_time.push_back(t);
	_cumIncidence.push_back(initInfectious);
	_prevalence.push_back(initInfectious);
	_nS.push_back(_popSize - initInfectious);
	_nR.push_back(0);
	
	// Simulation
	
	while ( t<horizon && at_least_one_S_and_E_or_I() )
	{
		// Draw the next event time based on the sum of specified rates
		double dt = drawNextEventInterval();
		double time_event = t+dt;
		
		// Draw the type of event based on each specified rates
		unsigned int et	= drawEventType();
		
		// --- Event is NOT a new infection
		if(et>0)
		{
			// retrieve all the relevant IDs
			vector<unsigned long> x = census_ID(et);
			
			// Pick one randomly
			unsigned long ID_selected = extractElementRandom(x);
			
			// Move this individual to the next stage
			_indiv[ID_selected].incrementStatus();
		}
		
		
		// --- Event is an new infection
		if(et==0){
			// retrieve all the susceptible IDs
			vector<unsigned long> x = census_ID(0);
			
			// Pick new infectee randomly
			unsigned long ID_infectee = extractElementRandom(x);
			
			// retrieve all the infectious IDs
			vector<unsigned long> x_I = census_ID_I();
			
			// Pick new infector randomly
			unsigned long ID_infector = extractElementRandom(x_I);
			
			
			// Now that we know who infected whom,
			// update all relevant quantities
			
			// > infectious status and infector ID
			set_infectiousStatus(ID_infectee, 1);
			set_infectorID(ID_infectee, ID_infector);
			
			// > timing of infections from both view points
			set_timeDiseaseAcquisition(ID_infectee, time_event);
			set_timeDiseaseTransmit(ID_infector, time_event);
			
			// > generation intervals
			double gi = time_event - get_timeDiseaseAcquisition(ID_infector);
			set_GIbck(ID_infectee, gi);
			set_GIfwd(ID_infector, gi);
		}
		
		// Update WIW matrix
		// (only at rounded times)
		if (calc_WIW_Re && t-dt<int(t) && int(t)<t)	calc_WIW(t);
		
		
		// update times
		_time.push_back(time_event);
		t += dt;
		
		// -- Keep track of S and I counts --
		// transmission:
		if(et==0) _count_S--;
		
		// from last E to I[1]:
		if(et==_nE) _count_I++;
		
		// from last I to R
		if(et==(_nE+_nI)){
			_count_I--;
			_count_R++;
		}
		
		// update prevalence time series
		_prevalence.push_back(_count_I);
		
		// update susceptible counts time series
		_nS.push_back(_count_S);
		_nR.push_back(_count_R);
		
		// calculate cumulative incidence
		double inc = (et==0)?1:0;
		inc += _cumIncidence[_cumIncidence.size()-1];
		_cumIncidence.push_back(inc);
		
	} // end-while-loop
	
	// Calculate time series of case effective reproductive number _Reff
	// (_WIW must be calculated until horizon)
	if(calc_WIW_Re) calc_Reff();
	
	// Integrity check
	if (calc_WIW_Re && _WIW.size()>0){
		unsigned long n_wiw = _WIW[_WIW.size()-1].countNonZeroElements();
		unsigned long n_cuminc =_cumIncidence[_cumIncidence.size()-1];
		stopif(n_wiw>n_cuminc, "cumulative incidences (_WIW vs _cumIncidence) not consistent!");
	}
}



void simulator::calc_WIW(double t)
{
	/// Calculate matrix of 'Who Infected Who' at different time points
	
	Matrix M(_popSize,_popSize);
	M.setAllValues(0.0);
	
	for(int j=0; j<_popSize; j++){
		// only if the individual has been infected
		if (_indiv[j].get_infectiousStatus()>0)
		{
			unsigned long ii = _indiv[j].get_infectorID();
			M(ii,j) = _indiv[j].get_GIbck();
		}
	}
	_WIW.push_back(M);
	_WIW_times.push_back(t);
}



void simulator::calc_Reff()
{
	/// Calculate the realized Effective Reproductive number
	/// It can be calculated only when WIW matrix is updated
	/// and at horizon of the simulation
	
	if (_WIW.size()>0)
	{
	 vector<double> n_2nd_cases;
	 vector<double> acq_time;
	 
	 for (unsigned long i=0; i<_popSize; i++)
	 {
		 if (_indiv[i].get_infectiousStatus()>0)
		 {
			 // counts the number of infectees for everyone
			 unsigned long c_i = _WIW[_WIW.size()-1].countNonZeroElements_line(i);
			 n_2nd_cases.push_back((double)c_i);
			 
			 // retrieve disease acquisition time:
			 acq_time.push_back(_indiv[i].get_timeDiseaseAcquisition());
		 }
	 }
		_Reff.addColVector(acq_time);
		_Reff.addColVector(n_2nd_cases);
	}
}


bool simulator::at_least_one_S_and_I()
{
	/// Tests if there is at least one susceptible AND one infectious
	/// individual in the whole population
	
	bool res_S = false;
	bool res_I = false;
	
	for(int i=0;i<_popSize; i++)
	{
		unsigned int istatus = _indiv[i].get_infectiousStatus();
		
		if (istatus==0) res_S = true;
		
		if (istatus>=_nE+1 &&
			istatus<=_nE+_nI) res_I = true;
		
		if(res_S && res_I) break;
	}
	
	return res_S*res_I;
}


bool simulator::at_least_one_S_and_E_or_I()
{
	/// Tests if there is at least one susceptible AND
	/// one infectious or latent individual in the whole population
	
	bool res_S = false;
	bool res_EI = false;
	
	for(int i=0;i<_popSize; i++)
	{
		unsigned int istatus = _indiv[i].get_infectiousStatus();
		
		if (istatus==0) res_S = true;
		
		if (istatus>= 1 &&
			istatus<=_nE+_nI) res_EI = true;
		
		if(res_S && res_EI) break;
	}
	
	return res_S*res_EI;
}


unsigned long simulator::census_status(unsigned int a, unsigned int b)
{
	/// counts individuals b/w infectious status a and b (both included)
	
	unsigned long cnt = 0;
	
	for(int i=0;i<_popSize; i++)
	{
		unsigned int istatus = _indiv[i].get_infectiousStatus();
		
		if (istatus>=a &&
			istatus<=b)
		{
			cnt++;
		}
	}
	return cnt;
}



unsigned long simulator::census_status(unsigned int a)
{
	/// counts individuals of infectious status a
	
	unsigned long cnt = 0;
	
	for(int i=0;i<_popSize; i++){
		if (_indiv[i].get_infectiousStatus()==a) cnt++;
	}
	return cnt;
}



unsigned long simulator::census_I()
{
	/// counts individuals in I[k] for all k
	
	return census_status(_nE+1,_nE+_nI);
}


unsigned long simulator::census_S()
{
	/// counts susceptible individuals
	
	return census_status(0);
}



vector<unsigned long> simulator::census_ID(unsigned int a)
{
	/// Retrieve IDs of all individuals of infectious status 'a'
	
	vector<unsigned long> res;
	
	for(int i=0;i<_popSize; i++)
	{
		if (_indiv[i].get_infectiousStatus()==a)
		{
			res.push_back(_indiv[i].get_ID());
		}
	}
	return res;
}



vector<unsigned long> simulator::census_ID_I()
{
	/// Retrieve IDs of all infectious individuals
	
	vector<unsigned long> res;
	
	for(int i=0;i<_popSize; i++)
	{
		if (_indiv[i].get_infectiousStatus()>=_nE+1 &&
			_indiv[i].get_infectiousStatus()<=_nE+_nI )
		{
			res.push_back(_indiv[i].get_ID());
		}
	}
	return res;
}





void simulator::set_GIbck(unsigned long IDindiv, double gi)
{
	/// Set backward generation interval to value 'gi' for individual ID 'IDindiv'
	
	for(int i=0; i<_popSize; i++)
	{
		if(_indiv[i].get_ID()==IDindiv) {
			_indiv[i].set_GIbck(gi);
			break;
		}
	}
	
}

void simulator::set_GIfwd(unsigned long IDindiv, double gi)
{
	/// Set backward generation interval to value 'gi' for individual ID 'IDindiv'
	
	for(int i=0; i<_popSize; i++)
	{
		if(_indiv[i].get_ID()==IDindiv){
			_indiv[i].set_GIfwd_incr(gi);
			break;
		}
	}
	
}



double simulator::get_timeDiseaseAcquisition(unsigned long ID)
{
	double res = -9.99;
	
	for(int i=0; i<_popSize; i++){
		if(_indiv[i].get_ID()==ID){
			res=_indiv[i].get_timeDiseaseAcquisition();
			break;
		}
	}
	return res;
}




void simulator::set_infectiousStatus(unsigned long IDindiv, unsigned int s)
{
	for(int i=0; i<_popSize; i++){
		if(_indiv[i].get_ID()==IDindiv){
			_indiv[i].set_infectiousStatus(s);
			break;
		}
	}
}



void simulator::set_infectorID(unsigned long IDinfectee, unsigned long IDinfector)
{
	for(int i=0; i<_popSize; i++){
		if(_indiv[i].get_ID()==IDinfectee){
			_indiv[i].set_infectorID(IDinfector);
			break;
		}
	}
}


void simulator::set_timeDiseaseAcquisition(unsigned long IDinfectee, double t)
{
	for(int i=0; i<_popSize; i++)
	{
		if(_indiv[i].get_ID()==IDinfectee){
			_indiv[i].set_timeDiseaseAcquisition(t);
			break;
		}
	}
}

void simulator::set_timeDiseaseTransmit(unsigned long IDinfector, double t)
{
	for(int i=0; i<_popSize; i++){
		if(_indiv[i].get_ID()==IDinfector){
			_indiv[i].set_timeDiseaseTransmit(t);
			break;
		}
	}
}



vector<unsigned long> simulator::census()
{
	/// counts everyone in each status
	
	vector<unsigned long> res;
	
	for(int s=0; s<= _nE+_nI+1; s++)
	{
		res.push_back(census_status(s));
	}
	
	return res;
}



vector<unsigned long> simulator::census_not_R()
{
	/// counts everyone in each status except 'removed' individuals
	
	vector<unsigned long> res;
	
	for(int s=0; s< _nE+_nI+1; s++)
	{
		res.push_back(census_status(s));
	}
	
	return res;
}





void simulator::save_GIbck(string filename){
	/// Save the temporal evolution of
	/// backward GI to a file

	ofstream f(filename);
	
	for (int i=0; i<_popSize;  i++){
		double tmp = _indiv[i].get_timeDiseaseAcquisition();
		
		if(tmp>0){
			f << tmp;
			f << ",";
			f << _indiv[i].get_GIbck() << endl;
		}
	}
}




void simulator::save_GIfwd(string filename){
	/// Save the temporal evolution of
	/// forward GI to a file
	
	ofstream f(filename);
	
	double R0 = _beta/_gamma[0]*_nI;
	// make sure the max length
	// of GIfwd is sufficient (for all MC iters!)
	unsigned long gimaxsize = 20*(1+R0);
	
	for (int i=0; i<_popSize;  i++)
		gimaxsize = max(gimaxsize,_indiv[i].get_GIfwd().size());
	
	// Now save all GI (even if there's not as much as 'gimaxsize')
	
	for (int i=0; i<_popSize;  i++)
	{
		double tmp = _indiv[i].get_timeDiseaseAcquisition();
		unsigned long ntransm = _indiv[i].get_GIfwd().size();
		
		if(ntransm>0){
			f << tmp <<",";
			
			for(int k=0;k<gimaxsize;k++)
			{
				if(k<ntransm) f<<_indiv[i].get_GIfwd()[k];
				if(k>=ntransm) f<<"";
				
				if(k<gimaxsize-1) f<<",";
			}
			f<<endl;
		}
	}
}



void simulator::save_prevalence(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _prevalence.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_prevalence[t]<<endl;
	}
}


void simulator::save_cumIncidence(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _cumIncidence.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_cumIncidence[t]<<endl;
	}
}


void simulator::save_nS(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _nS.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_nS[t]<<endl;
	}
}


void simulator::save_nR(string filename)
{
	ofstream f(filename);
	
	stopif(_time.size() != _nR.size(),"simulation inconsistent");
	
	for (int t=0; t<_time.size(); t++) {
		f<<_time[t]<<","<<_nR[t]<<endl;
	}
}



void simulator::save_Reff(string filename)
{
	stopif( (_Reff.getNbRows() ==0) && (_WIW.size()>0) ,"Reff not calculated --> cannot save associated data!");
	
	if(_Reff.getNbRows()>0)	_Reff.WriteToFileCSV(filename);
	
	if(_Reff.getNbRows()==0)
	{
		// Make sure the file is not completely empty...
		vector<double> tmp(2,0.0);
		_Reff.addRowVector(tmp);
		_Reff.WriteToFileCSV(filename);
	}
	
}