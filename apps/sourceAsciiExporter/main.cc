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
#include <Source.h>
#include <DS9Region.h>
#include <DS9RegionParser.h>
#include <SourceExporter.h>
#include <SourceSelector.h>

#include <ConfigParser.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>
#include <MathUtils.h>
#include <EllipseUtils.h>
#include <Contour.h>
#include <WCSUtils.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace Caesar;

void Usage(char* exeName)
{
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input=[INPUT_FILE] \t Input ROOT file produced by CAESAR containing the source collection to be read and written in a simpler ROOT TTree format"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (ROOT format) where to store source islands TTree (default=sources.dat)"<<endl;
	cout<<"-O, --output-comp=[OUTPUT_FILE] \t Output file name (ROOT format) where to store source components TTree (default=sourceComponents.dat)"<<endl;	
	cout<<"-b, --writeFluxInfoInBrightnessUnits \t Write flux information in brighteness units (Jy/beam) (default=no)"<<endl;
	cout<<"-a, --writeAdditionalSourceInfo \t Write additional source info (spectral index, cross-match objects) (default=no)"<<endl;
	cout<<"-w, --wcsType=[WCS_TYPE] \t WCS type to be used for sky coordinates (0=J2000,1=B1950,2=GALACTIC,3=ECLIPTIC,4=ALTAZ,5=LINEAR) (default=J2000)"<<endl;
	cout<<"-d, --delimiter=[ASCII_DELIMITER] \t Delimiter used to separate ascii columns (1=tab,2=comma,3=pipe) (default=tab)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;

}//close Usage()


static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "output", required_argument, 0, 'o' },
	{ "output-comp", required_argument, 0, 'O' },
	{ "writeFluxInfoInBrightnessUnits", no_argument, 0, 'b' },
	{ "writeAdditionalSourceInfo", no_argument, 0, 'a' },
	{ "wcsType", required_argument, 0, 'w' },
	{ "delimiter", required_argument, 0, 'd' },
	{ "verbosity", required_argument, 0, 'v'},
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string outputFileName= "sources.dat";
std::string outputFileName_comp= "sourceComponents.dat";
int verbosity= 4;//INFO level
bool writeAdditionalSourceInfo= false;
bool convertBrightnessToFlux= true;
int wcsType= eJ2000;
enum AsciiDelimiter 
{
	eTAB_DELIMITER=1,
	eCOMMA_DELIMITER=2,
	ePIPE_DELIMITER=3
};
char asciiDelimiter= '\t';
int delimiterType= eTAB_DELIMITER; 

//Globar vars
std::vector<Source*> m_sources;


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadSourceData(std::string filename);
int Save();

int main(int argc, char *argv[])
{
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
	//== Read source data
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading source data from file "<<fileName<<" ...");
	#endif
	if(ReadSourceData(fileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of source data failed!");
		#endif
		return -1;
	}

	
	//=======================
	//== Save
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving data to file ...");
	#endif
	Save();

	#ifdef LOGGING_ENABLED
		INFO_LOG("End source ascii exporter");
	#endif

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
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hi:o:O:abw:d:v:",options_tab, &option_index)) != -1) {
    
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
    	case 'i':	
			{
				fileName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'O':	
			{
				outputFileName_comp= std::string(optarg);	
				break;	
			}
			case 'a':	
			{
				writeAdditionalSourceInfo= true;
				break;	
			}
			case 'b':	
			{
				convertBrightnessToFlux= false;
				break;
			}
			case 'd':	
			{
				delimiterType= atoi(optarg);	
				break;	
			}
			case 'w':	
			{
				wcsType= atoi(optarg);	
				break;	
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
 
	
	//=======================
	//== Init Logger 
	//=======================
	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	#ifdef LOGGING_ENABLED
		LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	#endif

	//=======================
	//== CHECK ARGS 
	//=======================
	//Check input file name
	if(fileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty input file name given!");
		#endif
		return -1;
	}
	
	//Check output file name
	if(outputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output file name given!");
		#endif
		return -1;
	}

	//Check output file name
	if(outputFileName_comp==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output component file name given!");
		#endif
		return -1;
	}

	//Check delimiter type
	if(delimiterType==eTAB_DELIMITER){
		asciiDelimiter= '\t';
	}
	else if(delimiterType==eCOMMA_DELIMITER){
		asciiDelimiter= ',';
	}
	else if(delimiterType==ePIPE_DELIMITER){
		asciiDelimiter= '|';
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid or not recognized delimiter type given ("<<delimiterType<<")!");
		#endif
		return -1;
	}	


	return 0;

}//close ParseOptions()



int ReadSourceData(std::string filename)
{
	//Open files
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filename<<"!");
		#endif
		return -1;
	}

	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)inputFile->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get access to source tree in file "<<filename<<"!");	
		#endif
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);
	
	//Read sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	#endif
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);
		
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source no. "<<i+1<<"/"<<sourceTree->GetEntries()<<"...");
		#endif

		//Copy source
		Source* source= new Source;
		*source= *aSource;
		
		//Add sources to list
		m_sources.push_back(source);

	}//end loop sources

	return 0;

}//close ReadSourceData()



int Save()
{
	//Return if no sources are found
	if(m_sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source collection to be written is empty, no catalog file will be written!");
		#endif
		return 0;
	}

	//Retrieve source WCS
	WCS* wcs= m_sources[0]->GetWCS(wcsType);
	if(!wcs) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute WCS from sources!");
		#endif
		return -1;
	}	

	//Saving island/blob catalog to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source island catalog to file "<<outputFileName<<" ...");
	#endif
	bool dumpNestedSourceInfo= true;
	int status= SourceExporter::WriteToAscii(outputFileName,m_sources,dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightnessToFlux,asciiDelimiter);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source catalog to file "<<outputFileName<<" failed!");
		#endif
	}
	
	//Saving source fitted components to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source component catalog to file "<<outputFileName_comp<<" ...");
	#endif
	status= SourceExporter::WriteComponentsToAscii(outputFileName_comp,m_sources,dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightnessToFlux,asciiDelimiter);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source fitted component catalog to file "<<outputFileName_comp<<" failed!");
		#endif
	}

	return 0;

}//close Save()



std::string GetStringLogLevel(int verbosity)
{
	std::string slevel= "";
	if(verbosity<=0) slevel= "FATAL";
	else if(verbosity==1) slevel= "FATAL";
	else if(verbosity==2) slevel= "ERROR";
	else if(verbosity==3) slevel= "WARN";
	else if(verbosity==4) slevel= "INFO";
	else if(verbosity>=5) slevel= "DEBUG";
	else slevel= "OFF";

	return slevel;

}//close GetStringLogLevel()

