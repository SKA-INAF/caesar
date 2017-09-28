// **************************************************************************
// * License and Disclaimer                                                 *
// *                                                                        *
// * Copyright 2016 Simone Riggi																			      *
// *																																	      *
// * This file is part of SKA DSH.LMC 																		  *
// * SKA DSH.LMC is free software: you can redistribute it and/or modify it *
// * under the terms of the GNU General Public License as published by      *
// * the Free Software Foundation, either * version 3 of the License,       * 
// * or (at your option) any later version.                                 * 
// * Caesar is distributed in the hope that it will be useful, but 			    *
// * WITHOUT ANY WARRANTY; without even the implied warranty of             * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
// * See the GNU General Public License for more details. You should        * 
// * have received a copy of the GNU General Public License along with      * 
// * DSH.LMC. If not, see http://www.gnu.org/licenses/.                      *
// **************************************************************************
#include "gtest/gtest.h"

#include <TestEnvironment.h>


/*
#include <iostream>
#include <getopt.h>

using namespace std;

static const struct option options_tab[] = {
  { "datadir", required_argument, 0, 'd' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<"-d, --datadir \t Test data directory/path"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()
*/

int main(int argc, char **argv) {

	/*
	//Get options
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}
	int c = 0;
  int option_index = 0;

	std::string datadir;

	while((c = getopt_long(argc, argv, "d:",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}		
			case 'd':	
			{
				datadir= std::string(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while

	cout<<"*** ARGS ***"<<endl;
	cout<<"test data dir: "<<datadir<<endl;
	cout<<"*************"<<endl;
	*/

	//## Initialize GoogleTest
  ::testing::InitGoogleTest(&argc, argv);

	//Init test environment
	::testing::AddGlobalTestEnvironment(new TestEnvironment());

  return RUN_ALL_TESTS();
}
