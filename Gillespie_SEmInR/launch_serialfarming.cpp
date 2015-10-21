//
//  launch_serialfarming.cpp
//
//
//  Created by David CHAMPREDON on 2015-03-04.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "dcTools.h"
#include <string.h>

int main(int argc, const char * argv[]) {
	
	system("pwd");
	system("date");
	
	
	// Read main simulation parameters from file
	
	vector<string> param_list;
	vectorFromCSVfile_string(param_list, "param_all_list.csv", 1);
	
	// Loop on all parameter sets requested
	
	for(int s=0; s<param_list.size(); s++){
		
		string fileparam = param_list[s];
		cout<<"Launching parameter set #"<<s+1<<"/";
		cout<< param_list.size() <<" ("<<fileparam<<")"<<endl;
		
		int njobs = getParameterFromFile("njobs", fileparam);
		
		
		for (int j=1; j<=njobs; j++)
		{
			cout << endl << " > Launching job #"<<j<<"/"<<njobs<<endl;
			string cmd = "./gilSEIR_U " + int2string(j) + " " + fileparam + " &";
			cout<<cmd<<endl;
			system(cmd.c_str());
		}
	}
	
	return 0;
}
