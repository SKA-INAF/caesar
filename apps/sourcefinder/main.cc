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
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
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
	cout<<"-o, --options \t Show all defined configuration options"<<endl;
	cout<<"-v, --version \t Print software version information"<<endl;
	cout<<"-a, --authors \t Print authors information"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
	{ "options", no_argument, 0, 'o' },
	{ "version", no_argument, 0, 'v' },
	{ "authors", no_argument, 0, 'a' },
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
	//== Print logo
	//================================
	SysUtils::PrintAsciiLogo();

	//================================
	//== Parse command line options
	//================================
	if(ParseOptions(argc,argv)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to parse command line options!");
		#endif
		return -1;
	}
	
	//=======================
	//== Run SourceFinder
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Starting source finding");
	#endif
	SFinder* finder= new SFinder;

	int status= 0;

	if(finder->Run()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Source finding failed!");
		#endif
		status= 1;
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("End source finding");
		#endif
	}

	//=======================
	//== Clear
	//=======================
	if(finder){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Clear finder");
		#endif
		delete finder;
		finder= 0;
	}

	return status;

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

	while((c = getopt_long(argc, argv, "hc:ova",options_tab, &option_index)) != -1) {
    
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
			case 'o':
			{
				cout<<"== CONFIG OPTIONS =="<<endl;
				ConfigParser::Instance().PrintOptions();
				cout<<"========================"<<endl;
				exit(0);
			}
			case 'v':
			{
				std::string software_version= CAESAR_VERSION + std::string("_") + CAESAR_SHA_VERSION;
				cout<<"caesar "<<software_version<<endl;
				exit(0);
			}
			case 'a':
			{
				std::string authors_str= CAESAR_AUTHORS;
				std::string contacts_str= CAESAR_AUTHOR_CONTACTS;
		
				std::vector<std::string> authors= CodeUtils::SplitStringOnPattern(authors_str, ';');
				std::vector<std::string> contacts= CodeUtils::SplitStringOnPattern(contacts_str, ';');

				for(size_t i=0;i<authors.size();i++)
				{
					std::string author_raw= authors[i];
					std::vector<std::string> author_fields= CodeUtils::SplitStringOnPattern(author_raw, '-');
					std::string author= CodeUtils::JoinStringVec(author_fields, " ");
					std::string contact= contacts[i];
					cout<<author<<" ("<<contact<<")"<<endl;
				}

				exit(0);
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
	cout<<endl;
	cout<<"== RUN CONFIG OPTIONS =="<<endl;
	PRINT_OPTIONS();
	cout<<"========================"<<endl;
	cout<<endl;


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
	#ifdef LOGGING_ENABLED
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
	#endif

	//=======================
	//== Init thread numbers 
	//=======================
	int nThreads;
	if(GET_OPTION_VALUE(nThreads,nThreads)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get nThreads option!");
		#endif
		return -1;
	}
	if(nThreads>0) SysUtils::SetOMPThreads(nThreads);

	return 0;

}//close ParseOptions()

