// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************

#include <Image.h>
#include <SFinder.h>

#include <ConfigParser.h>
#include <Logger.h>
#include <CodeUtils.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace Caesar;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-c, --config \t Config file containing option settings"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string configFileName= "";

//Variables
//...

//Functions
int ParseOptions(int argc, char *argv[]);


int main(int argc, char *argv[])
{
	//================================
	//== Parse command line options
	//================================
	if(ParseOptions(argc,argv)<0){
		ERROR_LOG("Failed to parse command line options!");
		return -1;
	}
	
	//=======================
	//== Run SourceFinder
	//=======================
	INFO_LOG("Starting source finding");
	SFinder finder;
	if(finder.Run()<0){
		ERROR_LOG("Source finding failed!");
		return -1;
	}

	INFO_LOG("End source finding");

	return 0;

}//close main



int ParseOptions(int argc, char *argv[])
{
	
	//## Check args
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//## Parse options
	std::string configFileName= "";
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hc:",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}
			case 'h':
			{
      	Usage(argv[0]);	
				exit(0);
			}
    	case 'c':	
			{
				configFileName= std::string(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
 
	//## Read config options 
	if(ConfigParser::Instance().Parse(configFileName)<0){
		cerr<<"ERROR: Failed to parse config options!"<<endl;
		return -1;
	}
	PRINT_OPTIONS();

	//=======================
	//== Init Logger 
	//=======================
	//Get main logger options
	int loggerTarget= 0;
	if(GET_OPTION_VALUE(loggerTarget,loggerTarget)<0){
		cerr<<"ERROR: Failed to get loggerTarget option!"<<endl;
		return -1;
	}
	std::string loggerTag= "";
	std::string logLevel= "";
	if(GET_OPTION_VALUE(loggerTag,loggerTag)<0){
		cerr<<"ERROR: Failed to get loggerTag option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(logLevel,logLevel)<0){
		cerr<<"ERROR: Failed to get logLevel option!"<<endl;
		return -1;
	}

	//Init logger
	if(loggerTarget==eCONSOLE_TARGET){
		std::string consoleTarget= "";
		GET_OPTION_VALUE(consoleTarget,consoleTarget);
		LoggerManager::Instance().CreateConsoleLogger(logLevel,loggerTag,consoleTarget);
	}
	else if(loggerTarget==eFILE_TARGET){
		std::string logFile= "";
		std::string maxLogFileSize= "";
		bool appendToLogFile= false;
		int maxBackupLogFiles= 1;
		GET_OPTION_VALUE(logFile,logFile);
		GET_OPTION_VALUE(appendToLogFile,appendToLogFile);
		GET_OPTION_VALUE(maxLogFileSize,maxLogFileSize);
		GET_OPTION_VALUE(maxBackupLogFiles,maxBackupLogFiles);
		LoggerManager::Instance().CreateFileLogger(logLevel,loggerTag,logFile,appendToLogFile,maxLogFileSize,maxBackupLogFiles);
	}
	else if(loggerTarget==eSYSLOG_TARGET){
		std::string syslogFacility= "";
		GET_OPTION_VALUE(syslogFacility,syslogFacility);
		LoggerManager::Instance().CreateSysLogger(logLevel,loggerTag,syslogFacility);
	}
	else{
		cerr<<"ERROR: Failed to initialize logger!"<<endl;
		return -1;
	}

	//=======================
	//== Init thread numbers 
	//=======================
	int nThreads;
	if(GET_OPTION_VALUE(nThreads,nThreads)<0){
		ERROR_LOG("Failed to get nThreads option!");
		return -1;
	}
	if(nThreads>0) SysUtils::SetOMPThreads(nThreads);

	return 0;

}//close ParseOptions()

