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

#include <ConfigParser.h>
#include <Logger.h>
#include <CodeUtils.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>

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
	cout<<"-i, --input \t Input file name containing sources to be read in ROOT TTree (.root)"<<endl;
	cout<<"-I, --input2 \t Input file name 2 containing sources to be read in ROOT TTree"<<endl;
	cout<<"-o, --output \t Output file name "<<endl;
	cout<<"-t, --threshold \t Fraction of matching pixels to consider sources equal "<<endl;
	cout<<"-v, --verbosity \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "input2", required_argument, 0, 'I' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", optional_argument, 0, 'o' },
	{ "threshold", required_argument, 0, 't'},	
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string fileName2= "";
int verbosity= 4;//INFO level
float matchingThreshold= 0.9;

//Globar vars
TFile* outputFile= 0;
std::string outputFileName= "MatchOutput.root";
TTree* matchedSourceInfo;
std::vector<Source*> sources;
std::vector<Source*> sources2;
int SourceFoundFlag;
std::string SourceName;
int SourceType;
int SourceSimType;
double X0;
double Y0;
double S;
double Smax;
int SourceType2;
double S2;
double Smax2;
double MatchFraction;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int CorrelateSourceCatalogs();
int ReadSourceData();
void Save();
void Init();

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
	//== Init
	//=======================
	INFO_LOG("Initializing data...");
	Init();

	//=======================
	//== Read source data
	//=======================
	INFO_LOG("Reading source data from given files...");
	if(ReadSourceData()<0){
		ERROR_LOG("Reading of source data failed!");
		return -1;
	}

	//=======================
	//== Correlate catalogs
	//=======================
	INFO_LOG("Correlating source catalogs...");
	if(CorrelateSourceCatalogs()<0){
		ERROR_LOG("Correlating source data failed!");
		return -1;
	}

	//=======================
	//== Save
	//=======================
	INFO_LOG("Saving data to file ...");
	Save();

	INFO_LOG("End source catalog correlator");

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

	while((c = getopt_long(argc, argv, "hi:I:o:v:t:",options_tab, &option_index)) != -1) {
    
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
			case 'I':	
			{
				fileName2= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			case 't':
			{
				matchingThreshold= atof(optarg);
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
	LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	
	return 0;

}//close ParseOptions()


std::string GetStringLogLevel(int verbosity){

	std::string slevel= "";
	if(verbosity<=0) slevel= "FATAL";
	else if(verbosity==1) slevel= "FATAL";
	else if(verbosity==2) slevel= "ERROR";
	else if(verbosity==3) slevel= "WARN";
	else if(verbosity==4) slevel= "INFO";
	else if(verbosity>5) slevel= "DEBUG";
	else slevel= "OFF";

	return slevel;

}//close GetStringLogLevel()

int ReadSourceData()
{

	//Open files
	TFile* f= new TFile(fileName.c_str(),"READ");
	if(!f){
		ERROR_LOG("Failed to open file "<<fileName<<"!");
		return -1;
	}
	
	TFile* f2= new TFile(fileName2.c_str(),"READ");
	if(!f2){
		ERROR_LOG("Failed to open file "<<fileName2<<"!");
		return -1;
	}

	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)f->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName<<"!");	
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);

	TTree* sourceTree2= (TTree*)f2->Get("SourceInfo");
	if(!sourceTree2 || sourceTree2->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName2<<"!");	
		return -1;
	}
	sourceTree2->SetBranchAddress("Source",&aSource);

	//Read sources
	INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<fileName<<"...");
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	INFO_LOG("Reading #"<<sourceTree2->GetEntries()<<" sources in file "<<fileName2<<"...");
	for(int i=0;i<sourceTree2->GetEntries();i++){
		sourceTree2->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources2.push_back(source);
	}//end loop sources

	return 0;

}//close ReadSourceData()

int CorrelateSourceCatalogs()
{
	//Check source sizes
	if(sources.empty() || sources2.empty()){
		WARN_LOG("One or both source collections are empty, nothing to correlate...");
		return 0;
	}

	struct MatchingSourceInfo{	
		MatchingSourceInfo(float f,size_t index)
			: matchedPixelFraction(f), sourceIndex(index)
		{}
		double matchedPixelFraction;
		size_t sourceIndex;
	};

	INFO_LOG("Correlating ("<<sources.size()<<","<<sources2.size()<<") sources...");
	long int nFoundSources= 0;

	for(size_t i=0;i<sources.size();i++){
		long int NPix= sources[i]->GetNPixels();
		std::vector<MatchingSourceInfo> matched_info;	

		SourceFoundFlag= 0;
		SourceName= std::string(sources[i]->GetName());
		SourceType= sources[i]->Type;
		SourceSimType= sources[i]->SimType;
		X0= sources[i]->X0;
		Y0= sources[i]->Y0;
		S= sources[i]->GetS();
		Smax= sources[i]->GetSmax();
		S2= -1;
		Smax2= -1;
		SourceType2= -1;
		MatchFraction= 0.;

		for(size_t j=0;j<sources2.size();j++){
			long int NPix2= sources2[j]->GetNPixels();
			std::string SourceName2= std::string(sources2[j]->GetName());
			double X0_2= sources2[j]->X0;
			double Y0_2= sources2[j]->Y0;

			long int NMatchingPixels= sources[i]->GetNMatchingPixels(sources2[j]);
			double matchingPixelFraction= (double)(NMatchingPixels)/(double)(NPix);
			DEBUG_LOG("Source "<<SourceName<<" (X0="<<X0<<", Y0="<<Y0<<", N="<<NPix<<"): finding matching with source "<<SourceName2<<" (X0="<<X0_2<<", Y0="<<Y0_2<<", N="<<NPix2<<"), NMatchingPixels="<<NMatchingPixels<<" f="<<matchingPixelFraction<<" (t="<<matchingThreshold<<")");
			
			if(NMatchingPixels<=0 || matchingPixelFraction<matchingThreshold) continue;

			matched_info.push_back(MatchingSourceInfo(matchingPixelFraction,j));

		}//end loop 2nd collection

		//Fill source matching stats
		if(!matched_info.empty()){
			nFoundSources++;
			SourceFoundFlag= 1;	

			//Find best match in case of multiple matchings
			long int index_best= matched_info[0].sourceIndex;
			float matchFraction_best= matched_info[0].matchedPixelFraction; 
			
			if(matched_info.size()>1){
				for(size_t k=1;k<matched_info.size();k++){
					long int index= matched_info[k].sourceIndex;
					float matchFraction= matched_info[k].matchedPixelFraction; 
					if(matchFraction>matchFraction_best) {
						index_best= index;
						matchFraction_best= matchFraction;
					}
				}//end loop matchings
			}//close if

			SourceType2= sources2[index_best]->Type;
			S2= sources2[index_best]->GetS();
			Smax2= sources2[index_best]->GetSmax();
			MatchFraction= matchFraction_best;

		}//close if has match
	
		matchedSourceInfo->Fill();

	}//end loop collection

	INFO_LOG("#"<<nFoundSources<<"/"<<sources.size()<<" sources found...");

	return 0;

}//close CorrelateSourceCatalogs()


void Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	if(!matchedSourceInfo) matchedSourceInfo= new TTree("data","data");
	matchedSourceInfo->Branch("found",&SourceFoundFlag,"found/I");
	matchedSourceInfo->Branch("name",&SourceName);
	matchedSourceInfo->Branch("type",&SourceType,"type/I");
	matchedSourceInfo->Branch("simtype",&SourceSimType,"simtype/I");
	matchedSourceInfo->Branch("S",&S,"S/D");
	matchedSourceInfo->Branch("Smax",&Smax,"Smax/D");	
	matchedSourceInfo->Branch("X0",&X0,"X0/D");
	matchedSourceInfo->Branch("Y0",&Y0,"Y0/D");
	matchedSourceInfo->Branch("type2",&SourceType2,"type2/I");
	matchedSourceInfo->Branch("S2",&S2,"S2/D");
	matchedSourceInfo->Branch("Smax2",&Smax2,"Smax2/D");
	matchedSourceInfo->Branch("MatchFraction",&MatchFraction,"MatchFraction/D");

}//close Init()

void Save()
{
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		if(matchedSourceInfo) matchedSourceInfo->Write();
		outputFile->Close();
	}

}//close Save()


