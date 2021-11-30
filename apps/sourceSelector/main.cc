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
#include <TKey.h>
#include <TROOT.h>

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
	cout<<"-i, --input=[INPUT_FILE] \t Input ROOT file produced by CAESAR containing the source collection to be selected"<<endl;
	cout<<"-c, --cutfile=[CUT_FILE] \t Input ascii file containing the list of cuts to be applied to source collection"<<endl;
	cout<<"-a, --cutor \t Require OR of all cuts listed in cut file (e.g. at lease one passed) (default=require AND)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (ROOT format) where to store selected sources (default=sources.root)"<<endl;
	cout<<"-R, --region-output=[REGION_OUTPUT_FILE] \t Output DS9 region file name where to store selected sources (default=sources.reg)"<<endl;
	cout<<"-C, --catalog-output=[CATALOG_OUTPUT_FILE] \t Output catalog file name where to store selected sources (default=catalog.dat)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;

}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "cutfile", required_argument, 0, 'c' },
	{ "cutor", no_argument, 0, 'a' },
	{ "region-output", required_argument, 0, 'R' },
	{ "catalog-output", required_argument, 0, 'C' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", required_argument, 0, 'o' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string cutFileName= "";
bool requireAllCutsPassed= true;
std::string outputFileName= "sources.root";
std::string regionOutputFileName= "sources.reg";
std::string regionComponentsOutputFileName= "sources_fitcomp.reg";
std::string catalogOutputFileName= "catalog.dat";
std::string catalogComponentsOutputFileName= "catalog_fitcomp.dat";
int ds9WCSType= 0;//use original WCS type to save catalog
int verbosity= 4;//INFO level
TFile* inputFile= 0;
TTree* sourceTree= 0;
TTree* perfTree= 0;
TTree* configTree= 0;

//Globar vars
TFile* outputFile= 0;
TTree* outputTree= 0;

Source* m_source= 0;
std::vector<Source*> m_sources;
std::vector<Source*> m_sources_sel;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadRegionData(std::string filename);
int ReadSourceData(std::string filename);
int SelectSources();
int CloneObjectsInFile(std::vector<std::string> excludedObjNames);
void Save();
int SaveSources();
int SaveDS9Regions();
int SaveCatalog();
int Init();
void ClearData();

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
	//== Init
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initializing data...");
	#endif
	if(Init()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to initialize data!");
		#endif
		ClearData();
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
		ClearData();
		return -1;
	}

	//=======================
	//== Select source data
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Selecting sources located in region ...");
	#endif
	if(SelectSources()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Source selection failed!");
		#endif
		ClearData();
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
		INFO_LOG("End source region selector");
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

	while((c = getopt_long(argc, argv, "hi:c:ao:R:C:v:",options_tab, &option_index)) != -1) {
    
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
			case 'c':	
			{
				cutFileName= std::string(optarg);	
				break;	
			}	
			case 'a':
			{
				requireAllCutsPassed= false;
				break;
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'R':	
			{
				regionOutputFileName= std::string(optarg);	
				break;	
			}
			case 'C':	
			{
				catalogOutputFileName= std::string(optarg);	
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
	
	//Check cut file name
	if(cutFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty cut file name given!");
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

	//Check region output file name
	if(regionOutputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty region output file name given!");
		#endif
		return -1;
	}

	//Check catalog output file name
	if(catalogOutputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty catalog output file name given!");
		#endif
		return -1;
	}

	//Set DS9 region component file name
	regionComponentsOutputFileName= CodeUtils::ExtractFileNameFromPath(regionOutputFileName,true);
	regionComponentsOutputFileName+= "_fitcomp.reg";

	//Set catalog component file name
	catalogComponentsOutputFileName= CodeUtils::ExtractFileNameFromPath(catalogOutputFileName,true);
	catalogComponentsOutputFileName+= "_fitcomp.dat";

	return 0;

}//close ParseOptions()

int Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	//Init output source TTree
	if(!outputTree) outputTree= new TTree("SourceInfo","SourceInfo");
	m_source= 0;
	m_sources.clear();
	outputTree->Branch("Source",&m_source);
	
	return 0;

}//close Init()


void ClearData()
{
	//Delete TTree
	if(outputTree){
		delete outputTree;
		outputTree= 0;
	}

	//Close file
	if(outputFile && outputFile->IsOpen()){
		outputFile->Close();
	}

}//close ClearData()

int SelectSources()
{
	#ifdef LOGGING_ENABLED
		INFO_LOG("Applying cut selection to #"<<m_sources.size()<<" sources (requireAllCutsPassed? "<<requireAllCutsPassed<<") ...");
	#endif
	
	int status= SourceSelector::SelectSources(m_sources_sel,m_sources,cutFileName,requireAllCutsPassed);	
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to apply selection to sources through selector!");
		#endif
		return -1;
	}

	return 0;

}//close SelectSources()



int CloneObjectsInFile(std::vector<std::string> excludedObjNames)
{	
	//Check if input file is open
	if(!inputFile || !inputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Input file is not open!");
		#endif
		return -1;
	}

	//Check if output file is open
	if(!outputFile || !outputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Output file is not open!");
		#endif
		return -1;
	}

	//Loop on all entries present in input file
	TKey* key;
  TIter nextkey(inputFile->GetListOfKeys());
	
	while ((key = (TKey*)nextkey())) 
	{
		std::string className= key->GetClassName();
		std::string keyName= key->GetName();
    TClass* cl = gROOT->GetClass(className.c_str());
    if (!cl) continue;
	
		#ifdef LOGGING_ENABLED
			INFO_LOG("Processing object key "<<keyName<<" (name="<<className<<")...");
		#endif

		bool excludeObj= false;
		for(size_t i=0;i<excludedObjNames.size();i++){
			if(keyName==excludedObjNames[i]){
				excludeObj= true;
				break;
			}
		}
		if(excludeObj) {
			INFO_LOG("Object "<<keyName<<" exluded from the list of objects that will be saved ...");
			continue;
		}

		if (cl->InheritsFrom(TTree::Class())) {
    	TTree* T= (TTree*)inputFile->Get(key->GetName());

			//Write to output file
			outputFile->cd();
      TTree* newT= T->CloneTree(-1,"fast");
      newT->Write();
			
		}//close if TTree object
		else{
			TObject* obj= key->ReadObj();
      outputFile->cd();
      obj->Write();
      delete obj;			
		}

	}//end loop objects

	return 0;

}//close CloneObjectsInFile()


int ReadSourceData(std::string filename)
{
	//Open files
	inputFile= new TFile(filename.c_str(),"READ");
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
	for(int i=0;i<sourceTree->GetEntries();i++)
	{
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Getting entry source no. "<<i+1<<" from tree ...");
		#endif
		sourceTree->GetEntry(i);
		
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source no. "<<i+1<<"/"<<sourceTree->GetEntries()<<"...");
		#endif

		//Copy source
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Reading source no. "<<i+1<<" (name="<<aSource->GetName()<<") ...");
		#endif
		Source* source= new Source;
		*source= *aSource;
		
		//Add sources to list
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Adding source no. "<<i+1<<" (name="<<aSource->GetName()<<") to list ...");
		#endif
		m_sources.push_back(source);

	}//end loop sources

	return 0;

}//close ReadSourceData()


int SaveDS9Regions()
{
	//Save DS9 regions for islands
	bool convertDS9RegionsToWCS= false;
	int status= SourceExporter::WriteToDS9(regionOutputFileName,m_sources_sel,convertDS9RegionsToWCS);
	if(status<0){
		return -1;
	}

	//Save DS9 regions for components
	status= SourceExporter::WriteComponentsToDS9(regionComponentsOutputFileName,m_sources_sel,convertDS9RegionsToWCS);
	if(status<0){
		return -1;
	}

	return 0;

}//close SaveDS9Regions()

int SaveSources()
{
	//Check output file is open
	if(!outputFile || !outputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Output ROOT file not allocated or open!");
		#endif
		return -1;
	}
	if(!outputTree){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Output ROOT source TTree not allocated!");
		#endif
		return -1;
	}

	//Loop over selected sources and write TTree to output file
	for(size_t i=0;i<m_sources_sel.size();i++)
	{
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Fill source "<<i+1<<" (name="<<m_sources_sel[i]->GetName()<<" to TTree ...");
		#endif
		m_source= m_sources_sel[i];
		outputTree->Fill();
	}

	outputTree->Write();

	return 0;

}//close SaveSources()


int SaveCatalog()
{
	//Return if no sources are found
	if(m_sources_sel.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No sources selected, no catalog file will be written!");
		#endif
		return 0;
	}

	//Retrieve source WCS
	WCS* wcs= m_sources_sel[0]->GetWCS(ds9WCSType);
	if(!wcs) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute WCS from sources!");
		#endif
	}	

	//Saving island/blob catalog to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source catalog to file "<<catalogOutputFileName<<" ...");
	#endif
	bool dumpNestedSourceInfo= true;
	int status= SourceExporter::WriteToAscii(catalogOutputFileName,m_sources_sel,dumpNestedSourceInfo,ds9WCSType,wcs);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source catalog to file "<<catalogOutputFileName<<" failed!");
		#endif
	}
	
	//Saving source fitted components to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source catalog to file "<<catalogComponentsOutputFileName<<" ...");
	#endif
	status= SourceExporter::WriteComponentsToAscii(catalogComponentsOutputFileName,m_sources_sel,dumpNestedSourceInfo,ds9WCSType,wcs);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source fitted component catalog to file "<<catalogComponentsOutputFileName<<" failed!");
		#endif
	}

	return 0;

}//close SaveCatalog()

void Save()
{
	//Save TTree to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		
		//Clone all objects present in input file but the Source TTree
		CloneObjectsInFile({"SourceInfo"});

		//Write selected source TTree
		SaveSources();

		//Close file
		outputFile->Close();

	}//close if

	//Save DS9 regions with selected sources	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving DS9 regions with selected sources...");
	#endif
	SaveDS9Regions();

	//Save ascii catalog with selected sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving ascii catalog with selected sources...");
	#endif
	SaveCatalog();

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

