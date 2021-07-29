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
#include <StatsUtils.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TCut.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TKey.h>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


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
	cout<<"-I, --interactive \t Run application interactively"<<endl;
	cout<<"-t, --train \t If enabled, train classifier instead of applying it to input data (default=no)"<<endl;
	cout<<"-i, --input=[INPUT_FILE] \t Filename or list of input ROOT file(s) produced by CAESAR containing the source collection to be classified or used for training"<<endl;
	cout<<"-F, --filelist \t Consider input data as a filelist (default=no)"<<endl;
	cout<<"-w, --weights=[CLASSIFIER_WEIGHT_FILE] \t Classifier weight filename (.xml) used to apply classification to input data"<<endl;
	cout<<"-O, --options=[CLASSIFIER_OPTIONS] \t Options given to the classifier"<<endl;
	cout<<"-c, --classifier=[CLASSIFIER] \t Classifier method {Cuts,MLP} (default=MLP)"<<endl;
	cout<<"-n, --sigcut=[CLASSIFIER_SIGNAL_OUTPUT_CUT] \t Classifier signal/bkg cut value (default=0.5)"<<endl;
	cout<<"-e, --cuteff=[CUT_EFFICIENCY] \t Cut classifier signal efficiency (default=0.9)"<<endl;
	cout<<"-E, --eccentricity-diff \t Use eccentricity diff and not ratio (for circular beam cases) (default=no)"<<endl;
	cout<<"-T, --trainsize=[TRAIN_SAMPLE_FRACTION] \t Percentage of input data used for training the classifier (1-fract for testing) (default=0.5)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (ROOT format) where to store selected sources (default=sources.root)"<<endl;
	cout<<"-R, --region-output=[REGION_OUTPUT_FILE] \t Output DS9 region file name where to store selected sources (default=sources.reg)"<<endl;
	cout<<"-C, --catalog-output=[CATALOG_OUTPUT_FILE] \t Output catalog file name where to store selected sources (default=catalog.dat)"<<endl;
	cout<<"-f, --filterByType \t Consider only true sources with given type when searching the match (default=no)"<<endl;
	cout<<"-s, --selectedType=[TYPE] \t True source types to be crossmatched (1=COMPACT, 2=POINT-LIKE, 3=EXTENDED, 4=COMPACT_WITH_EXTENDED) (default=-1)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;

}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "train", no_argument, 0, 't' },
	{ "input", required_argument, 0, 'i' },
	{ "filelist", no_argument, 0, 'F' },
	{ "options", required_argument, 0, 'O' },
	{ "weights", required_argument, 0, 'w' },
	{ "classifier", required_argument, 0, 'c' },
	{ "sigcut", required_argument, 0, 'n' },
	{ "cuteff", required_argument, 0, 'e' },
	{ "eccentricity-diff", no_argument, 0, 'E' },
	{ "trainsize", required_argument, 0, 'T' },
	{ "region-output", required_argument, 0, 'R' },
	{ "catalog-output", required_argument, 0, 'C' },
	{ "filterByType", no_argument, 0, 'f'},
	{ "selectedType", required_argument, 0, 's'},	
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", required_argument, 0, 'o' },
	{ "interactive", no_argument, 0, 'I' },
  {(char*)0, (int)0, (int*)0, (int)0}
};

struct SourceComponentData
{
	int fitComponentIndex;
	int sourceIndex;
	int nestedSourceIndex;
	std::string sname;
	double peakSNR;
	double areaRatio;
	double eccentricityRatio;	
	double eccentricityDiff;
	int flag;
	int classOutput;
	bool isSelected;

	SourceComponentData(){	
		fitComponentIndex= -1;
		sourceIndex= -1;
		nestedSourceIndex= -1;
		sname= "";
		peakSNR= -999;
		areaRatio= -999;
		eccentricityRatio= -999;
		eccentricityDiff= -999;
		flag= 0;
		classOutput= 0;	
		isSelected= true;
	}
	SourceComponentData(int cindex,int sindex,int nindex):
		fitComponentIndex(cindex), sourceIndex(sindex), nestedSourceIndex(nindex)
	{
		sname= "";
		peakSNR= -999;
		areaRatio= -999;
		eccentricityRatio= -999;
		eccentricityDiff= -999;
		flag= 0;
		classOutput= 0;
		isSelected= true;
	}

};//close SourceComponentData()


//Options
bool gRunInteractively= false;
bool gTrainClassifier= false;
bool gIsFileList= false;
bool gApplyLogTransformToVariables= true;
bool gUseEccentricityDiff= false;
std::string fileName= "";
std::string weightFileName= "";
std::string outputFileName_trainData= "output_train.root";
std::string outputFileName= "sources.root";
std::string regionOutputFileName= "sources.reg";
std::string regionComponentsOutputFileName= "sources_fitcomp.reg";
std::string catalogOutputFileName= "catalog.dat";
std::string catalogComponentsOutputFileName= "catalog_fitcomp.dat";
bool selectSourceByType= false;//default=all true sources searched 
std::vector<int> stypes;
int ds9WCSType= 0;//use original WCS type to save catalog
int verbosity= 4;//INFO level
TFile* inputFile= 0;
TTree* sourceTree= 0;
TTree* perfTree= 0;
TTree* configTree= 0;
std::string gClassifierMethod= "MLP";//Possible choices {"MLP","Cuts"}
std::string gClassifierOptions= "";
std::string gCutClassifierOptions= "!H:!V:FitMethod=GA:EffSel";
std::string gNNClassifierOptions= "H:!V:NeuronType=tanh:EstimatorType=CE:VarTransform=None:NCycles=1000:HiddenLayers=N+6,N+3:TrainingMethod=BFGS:UseRegulator:CalculateErrors";

float gTrainSampleFraction= 0.5;
float gSignalCutEff= 0.9;
float gClassSignalOutCut= 0.5;
float gDataNormMin= -1;
float gDataNormMax= 1;
float gPeakSNRMin= 0;
float gPeakSNRMax= 10000;
float gEccentricityRatioMin= 0;
float gEccentricityRatioMax= 10;
float gEccentricityDiffMin= -1;
float gEccentricityDiffMax= 1;
float gSourceToBeamRatioMin= 0;
float gSourceToBeamRatioMax= 100;
float gLogPeakSNRMin= -0.5;
float gLogPeakSNRMax= 4;
float gLogEccentricityRatioMin= -2;
float gLogEccentricityRatioMax= 1;
float gLogEccentricityDiffMin= 0;
float gLogEccentricityDiffMax= 0.5;
float gLogSourceToBeamRatioMin= -4;
float gLogSourceToBeamRatioMax= 2;

//Globar vars
TFile* outputFile= 0;
TFile* outputFile_trainData= 0;
TTree* outputTree= 0;
TApplication* app= 0;
TMVA::Factory* gMVAFactory= 0;
TMVA::Reader* gMVAReader= 0;
TMVA::DataLoader* gDataLoader= 0;
TTree* gDataTree= 0;
TTree* gSignalDataTree_train= 0;
TTree* gBkgDataTree_train= 0;
TTree* gSignalDataTree_test= 0;
TTree* gBkgDataTree_test= 0;
float gPeakSNR;
float gEccentricityPar;
float gEccentricityRatio;
float gEccentricityDiff;
float gAreaRatio;
int gClassOutput;
int gSourceIndex;
int gNestedSourceIndex;
int gFitComponentIndex;
TCanvas* ClassOutPlot= 0;
TH2D* ClassOutPlotBkg= 0;
TH1D* ClassOutHisto_signal= 0;
TH1D* ClassOutHisto_bkg= 0;
TLegend* ClassOutPlotLegend= 0;

std::vector<SourceComponentData*> m_sourceComponentData;
Source* m_source= 0;
std::vector<Source*> m_sources;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadData(std::string filelist);
int ReadSourceData(std::string filename);
int FillSourcePars(std::vector<SourceComponentData*>& pars,Source* aSource,int sourceIndex,int nestedSourceIndex=-1);
int MakeClassifierData();
int TrainClassifier();
int EvaluateMLPClassifier(std::string weightFileName);
int EvaluateCutClassifier(std::string weightFileName);
int ApplyClassifier(std::string weightFileName);
int CloneObjectsInFile(std::vector<std::string> excludedObjNames);
void Save();
int SaveSources();
int SaveDS9Regions();
int SaveCatalog();
int Init();
void ClearData();
void SetStyle();

template<typename T>
int GetClassifierResponse(T& classRes,TMVA::Reader* reader,std::string classMethod)
{
	classRes= T(0);
	try{
		classRes= reader->EvaluateMVA(classMethod.c_str());
	}
	catch(...){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid classifier evaluation on data!");
		#endif
		return -1;
	}

	return 0;

}//close GetClassifierResponse()

int IsClassifiedAsSignal(bool& isSignal,TMVA::Reader* reader,std::string classMethod,double classSignalCut)
{
	isSignal= false;	
	int status= 0;
	if(classMethod=="Cuts") 
	{
		try{
			isSignal= reader->EvaluateMVA(classMethod.c_str());
		}
		catch(...){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Invalid classifier evaluation on data!");
			#endif
			return -1;
		}
	}//close if
	else{
		double classRes= 0;
		if(GetClassifierResponse(classRes,reader,classMethod)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Invalid classifier evaluation on data!");
			#endif
			return -1;
		}
		if(classRes>=classSignalCut) isSignal= true;
	}

	return 0;

}//close IsClassifiedAsSignal()


int main(int argc, char *argv[])
{
	//================================
	//== PARSE CMD LINE OPTIONS
	//================================
	if(ParseOptions(argc,argv)<0){
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Failed to parse command line options!");
		#endif
		return -1;
	}
	
	//=======================
	//== INIT DATA
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
	//== READ DATA
	//=======================
	if(gIsFileList){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading source data files from list "<<fileName<<" ...");
		#endif
		if(ReadData(fileName)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Reading of source data list failed!");
			#endif
			ClearData();
			return -1;
		}
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading source data file "<<fileName<<" ...");
		#endif
		if(ReadSourceData(fileName)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Reading of source data failed!");
			#endif
			ClearData();
			return -1;
		}
	}

	//===========================
	//== MAKE CLASSIFIER DATA
	//===========================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Creating classifier data ...");
	#endif
	if(MakeClassifierData()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Creating classifier data failed!");
		#endif
		ClearData();
		return -1;
	}
	

	//==============================================
	//== TRAIN/APPLY CLASSIFIER
	//==============================================
	if(gTrainClassifier){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Training classifier ...");
		#endif
		if(TrainClassifier()<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Classifier training failed!");
			#endif
			ClearData();
			return -1;
		}

	}//close if train classifier
	else{//apply classifier
		
		#ifdef LOGGING_ENABLED
			INFO_LOG("Applying classifier to input data ...");
		#endif
		if(ApplyClassifier(weightFileName)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Apply classifier failed!");
			#endif
			ClearData();
			return -1;
		}
	}

	//=======================
	//== SAVE DATA TO FILE
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving data to file ...");
	#endif
	Save();

	if(app && gRunInteractively) {
		#ifdef LOGGING_ENABLED
			INFO_LOG("End application run");	
		#endif
		gSystem->ProcessEvents();
		app->Run();
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("End source MVA classifier");
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

	while((c = getopt_long(argc, argv, "hti:Fw:o:R:C:v:fs:IO:c:n:e:ET:",options_tab, &option_index)) != -1) {
    
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
			case 'I':	
			{
				gRunInteractively= true;
				break;	
			}
			case 't':	
			{
				gTrainClassifier= true;
				break;	
			}
    	case 'i':	
			{
				fileName= std::string(optarg);	
				break;	
			}
			case 'F':
			{
				gIsFileList= true;
				break;
			}
			case 'w':	
			{
				weightFileName= std::string(optarg);	
				break;	
			}
			case 'c':	
			{
				gClassifierMethod= std::string(optarg);	
				break;	
			}
			case 'n':	
			{
				gClassSignalOutCut= atof(optarg);	
				break;	
			}	
			case 'e':	
			{
				gSignalCutEff= atof(optarg);	
				break;	
			}
			case 'T':	
			{
				gTrainSampleFraction= atof(optarg);	
				break;	
			}
			case 'O':	
			{
				gClassifierOptions= std::string(optarg);	
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
			case 'f':
			{
				selectSourceByType= true;
				break;
			}
			case 's':	
			{
				int stype= atoi(optarg);
				stypes.push_back(stype);	
				break;	
			}	
			case 'E':
			{
				gUseEccentricityDiff= true;
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

	//Check classifier weights
	if(!gTrainClassifier && weightFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty classifier weight file name given!");
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

	//Set classifier options
	if(gClassifierOptions==""){
		if(gClassifierMethod=="MLP"){
			gClassifierOptions= gNNClassifierOptions;
		}
		else if(gClassifierMethod=="Cuts"){
			gClassifierOptions= gCutClassifierOptions;
		}
		else{
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Invalid/unknown classifier method ("<<gClassifierMethod<<") given!");
			#endif
			return -1;
		}
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("classifierOptions: "<<gClassifierOptions);
	#endif
	

	//=======================
	//== INIT ROOT APP
	//=======================
	if(!app && gRunInteractively){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Initializing ROOT app...");
		#endif
		app= new TApplication("App",&argc,argv);
	}

	return 0;

}//close ParseOptions()

int Init()
{
	//Set draw & graphics style
	SetStyle();

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	//Init output source TTree
	if(!outputTree) outputTree= new TTree("SourceInfo","SourceInfo");
	m_source= 0;
	m_sources.clear();
	outputTree->Branch("Source",&m_source);

	
	//Init MVA 
	TMVA::Tools::Instance();

	//Init output train data file
	if(gTrainClassifier && !gRunInteractively) outputFile_trainData= new TFile(outputFileName_trainData.c_str(),"RECREATE");

	//Init classifier data TTree
	if(!gDataTree) gDataTree= new TTree("data","data");
	gDataTree->Branch("peakSNR",&gPeakSNR);
	gDataTree->Branch("eccentricityPar",&gEccentricityPar);
	//if(gUseEccentricityDiff) gDataTree->Branch("eccentricityPar",&gEccentricityDiff);
	//else gDataTree->Branch("eccentricityPar",&gEccentricityRatio);
	gDataTree->Branch("areaRatio",&gAreaRatio);	
	gDataTree->Branch("classOutput",&gClassOutput);
	gDataTree->Branch("sourceIndex",&gSourceIndex);
	gDataTree->Branch("nestedSourceIndex",&gNestedSourceIndex);
	gDataTree->Branch("fitComponentIndex",&gFitComponentIndex);

	gDataTree->SetBranchAddress("peakSNR",&gPeakSNR);
	gDataTree->SetBranchAddress("eccentricityPar",&gEccentricityPar);
	//if(gUseEccentricityDiff) gDataTree->SetBranchAddress("eccentricityPar",&gEccentricityDiff);
	//else gDataTree->SetBranchAddress("eccentricityPar",&gEccentricityRatio);
	gDataTree->SetBranchAddress("areaRatio",&gAreaRatio);
	gDataTree->SetBranchAddress("classOutput",&gClassOutput);
	gDataTree->SetBranchAddress("sourceIndex",&gSourceIndex);
	gDataTree->SetBranchAddress("nestedSourceIndex",&gNestedSourceIndex);
	gDataTree->SetBranchAddress("fitComponentIndex",&gFitComponentIndex);

	
	if(!gSignalDataTree_train) gSignalDataTree_train= new TTree("signalData_train","signalData_train");
	gSignalDataTree_train->Branch("peakSNR",&gPeakSNR);
	gSignalDataTree_train->Branch("eccentricityPar",&gEccentricityPar);
	//if(gUseEccentricityDiff) gSignalDataTree_train->Branch("eccentricityPar",&gEccentricityDiff);
	//else gSignalDataTree_train->Branch("eccentricityPar",&gEccentricityRatio);
	gSignalDataTree_train->Branch("areaRatio",&gAreaRatio);	
	gSignalDataTree_train->Branch("classOutput",&gClassOutput);
	

	if(!gBkgDataTree_train) gBkgDataTree_train= new TTree("bkgData_train","bkgData_train");
	gBkgDataTree_train->Branch("peakSNR",&gPeakSNR);
	gBkgDataTree_train->Branch("eccentricityPar",&gEccentricityPar);
	//if(gUseEccentricityDiff) gBkgDataTree_train->Branch("eccentricityPar",&gEccentricityDiff);
	//else gBkgDataTree_train->Branch("eccentricityPar",&gEccentricityRatio);
	gBkgDataTree_train->Branch("areaRatio",&gAreaRatio);	
	gBkgDataTree_train->Branch("classOutput",&gClassOutput);

	if(!gSignalDataTree_test) gSignalDataTree_test= new TTree("signalData_test","signalData_test");
	gSignalDataTree_test->Branch("peakSNR",&gPeakSNR);
	gSignalDataTree_test->Branch("eccentricityPar",&gEccentricityPar);
	//if(gUseEccentricityDiff) gSignalDataTree_test->Branch("eccentricityPar",&gEccentricityDiff);
	//else gSignalDataTree_test->Branch("eccentricityPar",&gEccentricityRatio);
	gSignalDataTree_test->Branch("areaRatio",&gAreaRatio);	
	gSignalDataTree_test->Branch("classOutput",&gClassOutput);

	if(!gBkgDataTree_test) gBkgDataTree_test= new TTree("bkgData_test","bkgData_test");
	gBkgDataTree_test->Branch("peakSNR",&gPeakSNR);
	gBkgDataTree_test->Branch("eccentricityPar",&gEccentricityPar);
	//if(gUseEccentricityDiff) gBkgDataTree_test->Branch("eccentricityPar",&gEccentricityDiff);
	//else gBkgDataTree_test->Branch("eccentricityPar",&gEccentricityRatio);
	gBkgDataTree_test->Branch("areaRatio",&gAreaRatio);	
	gBkgDataTree_test->Branch("classOutput",&gClassOutput);
	
	if(!ClassOutHisto_signal) ClassOutHisto_signal= new TH1D("ClassOutHisto_signal","ClassOutHisto_signal",100,-2,2);
	ClassOutHisto_signal->SetFillColor(kRed);
	ClassOutHisto_signal->SetFillStyle(3001);
	ClassOutHisto_signal->SetLineColor(kRed);

	if(!ClassOutHisto_bkg) ClassOutHisto_bkg= new TH1D("ClassOutHisto_bkg","ClassOutHisto_bkg",100,-2,2);
	ClassOutHisto_bkg->SetFillColor(kBlack);
	ClassOutHisto_bkg->SetFillStyle(3004);
	ClassOutHisto_bkg->SetLineColor(kBlack);

	//Set random seed
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);

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

	//Close train data file
	if(outputFile_trainData && outputFile_trainData->IsOpen()){
		outputFile_trainData->Close();
	}

}//close ClearData()



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
			#ifdef LOGGING_ENABLED
				INFO_LOG("Object "<<keyName<<" exluded from the list of objects that will be saved ...");
			#endif
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


int MakeClassifierData()
{
	//Init classifier data
	int nVars= 3;
	std::vector<std::vector<double>> dataVarList;
	std::vector<std::vector<double>> classifierData;
	std::vector<double> dataVars_fixed_min;
	std::vector<double> dataVars_fixed_max;
	if(gUseEccentricityDiff){
		dataVars_fixed_min= {gPeakSNRMin,gEccentricityDiffMin,gSourceToBeamRatioMin};
		dataVars_fixed_max= {gPeakSNRMax,gEccentricityDiffMax,gSourceToBeamRatioMax};
	}
	else{
		dataVars_fixed_min= {gPeakSNRMin,gEccentricityRatioMin,gSourceToBeamRatioMin};
		dataVars_fixed_max= {gPeakSNRMax,gEccentricityRatioMax,gSourceToBeamRatioMax};
	}

	if(gApplyLogTransformToVariables){
		dataVars_fixed_min.clear();
		dataVars_fixed_max.clear();
		if(gUseEccentricityDiff){	
			dataVars_fixed_min= {gLogPeakSNRMin,gLogEccentricityDiffMin,gLogSourceToBeamRatioMin};
			dataVars_fixed_max= {gLogPeakSNRMax,gLogEccentricityDiffMax,gLogSourceToBeamRatioMax};
		}
		else{
			dataVars_fixed_min= {gLogPeakSNRMin,gLogEccentricityRatioMin,gLogSourceToBeamRatioMin};
			dataVars_fixed_max= {gLogPeakSNRMax,gLogEccentricityRatioMax,gLogSourceToBeamRatioMax};
		}
	}

	for(int j=0;j<nVars;j++){
		dataVarList.push_back( std::vector<double>() );
	}

	//Store classifier data for normalization
	for(size_t i=0;i<m_sourceComponentData.size();i++)
	{
		double peakSNR= m_sourceComponentData[i]->peakSNR;
		double eccentricityRatio= m_sourceComponentData[i]->eccentricityRatio;
		double eccentricityDiff= m_sourceComponentData[i]->eccentricityDiff;
		double sourceToBeamRatio= m_sourceComponentData[i]->areaRatio;
		int classOutput= m_sourceComponentData[i]->classOutput;
		int componentIndex= m_sourceComponentData[i]->fitComponentIndex;
		int sindex= m_sourceComponentData[i]->sourceIndex;
		int nestedSourceIndex= m_sourceComponentData[i]->nestedSourceIndex;

		if(gApplyLogTransformToVariables){
			peakSNR= log10(peakSNR);
			eccentricityRatio= log10(eccentricityRatio);
			eccentricityDiff= log10(eccentricityDiff+2);//add +2 to avoid log of negative values (range would be 0,3)
			sourceToBeamRatio= log10(sourceToBeamRatio);
		}
		
		dataVarList[0].push_back(peakSNR);
		if(gUseEccentricityDiff) dataVarList[1].push_back(eccentricityDiff);
		else dataVarList[1].push_back(eccentricityRatio);
		dataVarList[2].push_back(sourceToBeamRatio);
		
		classifierData.push_back( std::vector<double>() );
		classifierData[i].push_back(peakSNR);
		if(gUseEccentricityDiff) classifierData[i].push_back(eccentricityDiff);
		else classifierData[i].push_back(eccentricityRatio);
		classifierData[i].push_back(sourceToBeamRatio);
		classifierData[i].push_back(classOutput);
		classifierData[i].push_back(sindex);
		classifierData[i].push_back(nestedSourceIndex);
		classifierData[i].push_back(componentIndex);

	}//end loop data

	//Compute data var stats
	std::vector<double> dataVars_min;
	std::vector<double> dataVars_max;
	std::vector<Caesar::BoxStats<double>> dataVars_stats;

	for(int j=0;j<nVars;j++){	
		BoxStats<double> stats= StatsUtils::ComputeBoxStats(dataVarList[j],false);
		double wmin= stats.minVal;
		double wmax= stats.maxVal;
		double wmin_fixed= dataVars_fixed_min[j];
		double wmax_fixed= dataVars_fixed_max[j];
		dataVars_min.push_back(wmin);
		dataVars_max.push_back(wmax);

		#ifdef LOGGING_ENABLED
			INFO_LOG("Data var "<<j+1<<": range ("<<wmin<<","<<wmax<<"), norm values ("<<wmin_fixed<<","<<wmax_fixed<<"), norm range ("<<gDataNormMin<<","<<gDataNormMax<<")");
		#endif

		dataVars_stats.push_back(stats);
	}//end loop vars

	// - Normalize input data
	#ifdef LOGGING_ENABLED
		INFO_LOG("Normalizing input data (#"<<classifierData.size()<<" data events) ...");
	#endif

	for(size_t i=0;i<classifierData.size();i++){//loop on events
		for(size_t j=0;j<classifierData[i].size()-4;j++){//loop on data vars
			//double wmin= dataVars_min[j];
			//double wmax= dataVars_max[j];
			double wmin= dataVars_fixed_min[j];
			double wmax= dataVars_fixed_max[j];

			double w= classifierData[i][j];
			double w_norm= gDataNormMin + (gDataNormMax-gDataNormMin)*(w-wmin)/(wmax-wmin);
			classifierData[i][j]= w_norm;
		}//end loop data vars
	}//end loop data

	// - Fill classifier input data TTree
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting classifier input data (#"<<classifierData.size()<<" data entries present) ...");
	#endif
	
	if(gTrainClassifier)
	{
		for(size_t i=0;i<classifierData.size();i++)
		{	
			gPeakSNR= classifierData[i][0];
			gEccentricityPar= classifierData[i][1];
			gAreaRatio= classifierData[i][2];
			gClassOutput= static_cast<int>(classifierData[i][3]);
			gSourceIndex= static_cast<int>(classifierData[i][4]);
			gNestedSourceIndex= static_cast<int>(classifierData[i][5]);
			gFitComponentIndex= static_cast<int>(classifierData[i][6]);

			double rand= gRandom->Uniform(0,1);
			if(rand<=gTrainSampleFraction)
			{
				if(gClassOutput==0) gBkgDataTree_train->Fill();
				else gSignalDataTree_train->Fill();
			}
			else{
				if(gClassOutput==0) gBkgDataTree_test->Fill();
				else gSignalDataTree_test->Fill();
			}
			gDataTree->Fill();

		}//end loop data

		long int nData= gDataTree->GetEntries();
		long int nSignal_train= gSignalDataTree_train->GetEntries();
		long int nSignal_test= gSignalDataTree_test->GetEntries();
		long int nBkg_train= gBkgDataTree_train->GetEntries();
		long int nBkg_test= gBkgDataTree_test->GetEntries();
		long int nTrain= nSignal_train + nBkg_train;
		long int nTest= nSignal_test + nBkg_test;  
	
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<gDataTree->GetEntries()<<" classifier data: TRAIN="<<nTrain<<" (signal="<<nSignal_train<<", bkg="<<nBkg_train<<"), TEST="<<nTrain<<" (signal="<<nSignal_test<<", bkg="<<nBkg_test<<")");
		#endif
	}//close if train mode
	else
	{
		for(size_t i=0;i<classifierData.size();i++)
		{	
			gPeakSNR= classifierData[i][0];
			gEccentricityPar= classifierData[i][1];
			gAreaRatio= classifierData[i][2];
			gClassOutput= static_cast<int>(classifierData[i][3]);
			gSourceIndex= static_cast<int>(classifierData[i][4]);
			gNestedSourceIndex= static_cast<int>(classifierData[i][5]);
			gFitComponentIndex= static_cast<int>(classifierData[i][6]);

			gDataTree->Fill();

		}//end loop data

		long int nData= gDataTree->GetEntries();
		
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<gDataTree->GetEntries()<<" sources will be given to classifier...");
		#endif
	}//close else


	return 0;

}//close MakeClassifierData()


int TrainClassifier()
{
	//######################################
	//##   INITIALIZE MVA
	//######################################
	//Initialize factory
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initializing MVA factory...");
	#endif
	if(!gMVAFactory){
		gMVAFactory= new TMVA::Factory(
			"SourceClassifier", 
			outputFile,
			"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification"
		);
	}

	//######################################
	//##   SET DATA
	//######################################
	//Initialize data loader
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting classifier data in MVA data loader ...");
	#endif
	
	gDataLoader= new TMVA::DataLoader("classifierData");
	gDataLoader->AddVariable("peakSNR", 'F');
	//gDataLoader->AddVariable("eccentricityRatio", 'F');
	gDataLoader->AddVariable("eccentricityPar", 'F');
	gDataLoader->AddVariable("areaRatio", 'F');	
	gDataLoader->AddSpectator("classOutput", 'I');	
	
	//Set the classifier train & test data
	//Assume weights=1
	double sigWeight = 1.0;
	double bkgWeight = 1.0;
	gDataLoader->AddSignalTree(gSignalDataTree_train, sigWeight, TMVA::Types::kTraining);
	gDataLoader->AddBackgroundTree(gBkgDataTree_train, bkgWeight, TMVA::Types::kTraining);
	gDataLoader->AddSignalTree(gSignalDataTree_test, sigWeight, TMVA::Types::kTesting);
	gDataLoader->AddBackgroundTree(gBkgDataTree_test, bkgWeight, TMVA::Types::kTesting);
	
	TCut emptyCuts= "";
	gDataLoader->PrepareTrainingAndTestTree(emptyCuts,"NormMode=None:!V");
	//gDataLoader->PrepareTrainingAndTestTree(emptyCuts,"SplitMode=Random:NormMode=NumEvents:!V");

	//######################################
	//##   SET CLASSIFIERS
	//######################################
	//Select the classifier method
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting classifier method  ...");
	#endif
	if(gClassifierMethod=="MLP"){
		gMVAFactory->BookMethod(gDataLoader,TMVA::Types::kMLP, "MLP", gClassifierOptions.c_str());
	}
	else if(gClassifierMethod=="Cuts"){
		gMVAFactory->BookMethod(gDataLoader,TMVA::Types::kCuts,"Cuts", gClassifierOptions.c_str());
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid or unknown classifier method ("<<gClassifierMethod<<") given!");
		#endif
		return -1;
	}

	//######################################
	//##   TRAIN CLASSIFIERS
	//######################################
	//Train classifier	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Starting classifier training ...");
	#endif
	gMVAFactory->TrainAllMethods();

	//######################################
	//##   TEST CLASSIFIERS
	//######################################
	// Evaluate all MVAs using the set of test events
	#ifdef LOGGING_ENABLED
		INFO_LOG("Evaluating trained classifier using the test sample ...");
	#endif
  gMVAFactory->TestAllMethods();

	// Evaluate and compare performance of all configured MVAs
	#ifdef LOGGING_ENABLED
		INFO_LOG("Evaluating and comparing trained classifier(s) (if more than one booked) ...");
	#endif
  gMVAFactory->EvaluateAllMethods();

	//Evaluate method and draw results
	#ifdef LOGGING_ENABLED
		INFO_LOG("Evaluating and drawing results ...");
	#endif
	int status= 0;
	std::string weightFileName= "";
	if(gClassifierMethod=="MLP"){
		weightFileName= "classifierData/weights/SourceClassifier_MLP.weights.xml";
		status= EvaluateMLPClassifier(weightFileName);
	}
	else if(gClassifierMethod=="Cuts"){
		weightFileName= "classifierData/weights/SourceClassifier_Cuts.weights.xml";
		status= EvaluateCutClassifier(weightFileName);
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid or unknown classifier method ("<<gClassifierMethod<<") given!");
		#endif
		return -1;
	}

	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to evaluate classifier performances!");
		#endif
		return -1;
	}
	

	return 0;

}//close TrainClassifier()


int EvaluateMLPClassifier(std::string weightFileName)
{
	//######################################
	//##   INIT MVA READER
	//######################################
	//Initialize MVA reader
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initializing MVA reader ...");
	#endif
  if(!gMVAReader){
		gMVAReader= new TMVA::Reader("!Color:!Silent");
	}
	gMVAReader->AddVariable("peakSNR",&gPeakSNR);
	//gMVAReader->AddVariable("eccentricityRatio",&gEccentricityRatio);
	gMVAReader->AddVariable("eccentricityPar",&gEccentricityPar);
	gMVAReader->AddVariable("areaRatio",&gAreaRatio);
	gMVAReader->AddSpectator("classOutput",&gClassOutput);	

	//Setting weights in reader
	//TString weightFile= "dataset/weights/SourceNNClassification_MLP.weights.xml";
	TString weightFile= weightFileName.c_str();
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting weights from file "<<weightFile.Data()<<" ...");
	#endif
	gMVAReader->BookMVA("MLP",weightFile);
	

	//######################################
	//##   INIT DATA
	//######################################
	

	int NSig_true= 0;
	int NSig= 0;
	int NSig_right= 0;
	int NSig_wrong= 0;
	int NBkg= 0;
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<gDataTree->GetEntries()<<" data will be read to evaluate classifier performances...");
	#endif
	
	for (int i=0;i<gDataTree->GetEntries();i++) {
		gDataTree->GetEntry(i);

		if(i%1000==0) cout<<"--> Reading "<<i+1<<"/"<<gDataTree->GetEntries()<<" event..."<<endl;

		double NNOut= 0;
		try{
			NNOut= gMVAReader->EvaluateMVA("MLP");
		}
		catch(...){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Invalid classifier evaluation on data "<<i+1<<"!");
			#endif
			continue;
		}

		//Fill NN out histo
		if(gClassOutput==1){
			NSig_true++;
			ClassOutHisto_signal->Fill(NNOut);
		}
		else{
			ClassOutHisto_bkg->Fill(NNOut);
		}

		//Update efficiency
		if(NNOut>=gClassSignalOutCut){//identified as signal
			NSig++;
			if(gClassOutput==1) NSig_right++;//correct identification
			else NSig_wrong++;//wrong identification
		}
		else{
			NBkg++;
		}

	}//end loop events

	//Compute efficiency & purity
	double Eff= (double)NSig_right/(double)NSig_true;
	double Purity= (double)NSig_right/(double)NSig;
	#ifdef LOGGING_ENABLED
		INFO_LOG("NSig_true="<<NSig_true<<", NSig_right="<<NSig_right<<", NSig="<<NSig<<", Efficiency="<<Eff<<", Purity="<<Purity);
	#endif
	
	//Draw NN output histos
	if(outputFile_trainData) outputFile_trainData->cd();

		ClassOutPlot= new TCanvas("ClassOutPlot","ClassOutPlot",600,600);
		ClassOutPlot->cd();

		ClassOutPlotBkg= new TH2D("ClassOutPlotBkg","",100,-0.5,1.5,100,0,1);
		ClassOutPlotBkg->SetXTitle("NNOut");
		ClassOutPlotBkg->SetYTitle("entries");
		ClassOutPlotBkg->SetStats(0);
		ClassOutPlotBkg->Draw();

		gPad->Modified(); 
		gPad->Update();

		ClassOutHisto_signal->DrawNormalized("hist same");

		gPad->Modified(); 
		gPad->Update();

		ClassOutHisto_bkg->DrawNormalized("hist same");

		gPad->Modified(); 
		gPad->Update();

		ClassOutPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
		ClassOutPlotLegend->SetFillColor(0);
		ClassOutPlotLegend->SetTextSize(0.045);
		ClassOutPlotLegend->SetTextFont(52);
		ClassOutPlotLegend->AddEntry(ClassOutHisto_signal,"real","F");
		ClassOutPlotLegend->AddEntry(ClassOutHisto_bkg,"false","F");	
		ClassOutPlotLegend->Draw("same");

		gPad->Modified(); 
		gPad->Update();

		ClassOutPlot->Update();
  	ClassOutPlot->Draw();

	
	if(gRunInteractively) gSystem->ProcessEvents();

	//Clear stuff
	if(gMVAReader) {
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Deleting MVA reader ...");
		#endif
		delete gMVAReader;
		gMVAReader= 0;
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("done");
		#endif
	}

	return 0;

}//close EvaluateMLPClassifier()

int EvaluateCutClassifier(std::string weightFileName)
{
	//######################################
	//##   INIT MVA READER
	//######################################
	//Initialize MVA reader
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initializing MVA reader ...");
	#endif
  if(!gMVAReader){
		gMVAReader= new TMVA::Reader("!Color:!Silent");
	}
	gMVAReader->AddVariable("peakSNR",&gPeakSNR);
	//gMVAReader->AddVariable("eccentricityRatio",&gEccentricityRatio);
	gMVAReader->AddVariable("eccentricityPar",&gEccentricityPar);
	gMVAReader->AddVariable("areaRatio",&gAreaRatio);
	gMVAReader->AddSpectator("classOutput",&gClassOutput);	

	//Setting weights in reader
	TString weightFile= weightFileName.c_str();
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting weights from file "<<weightFile.Data()<<" ...");
	#endif
	gMVAReader->BookMVA("Cuts",weightFile);
	

	//######################################
	//##   INIT DATA
	//######################################
	
	int NSig_true= 0;
	int NSig= 0;
	int NSig_right= 0;
	int NSig_wrong= 0;
	int NBkg= 0;
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<gDataTree->GetEntries()<<" data will be read to evaluate classifier performances...");
	#endif
	
	for (int i=0;i<gDataTree->GetEntries();i++) {
		gDataTree->GetEntry(i);

		if(i%1000==0) cout<<"--> Reading "<<i+1<<"/"<<gDataTree->GetEntries()<<" event..."<<endl;

		bool passed = true;
		try{
			passed= gMVAReader->EvaluateMVA("Cuts",gSignalCutEff);
		}
		catch(...){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Invalid classifier evaluation on data "<<i+1<<"!");
			#endif
			continue;
		}

		//Fill NN out histo
		if(gClassOutput==1){
			NSig_true++;
		}
		
		//Update efficiency
		if(passed){//identified as signal
			NSig++;
			if(gClassOutput==1) NSig_right++;//correct identification
			else NSig_wrong++;//wrong identification
		}
		else{
			NBkg++;
		}

	}//end loop events

	//Compute efficiency & purity
	double Eff= (double)NSig_right/(double)NSig_true;
	double Purity= (double)NSig_right/(double)NSig;
	#ifdef LOGGING_ENABLED
		INFO_LOG("NSig_true="<<NSig_true<<", NSig_right="<<NSig_right<<", NSig="<<NSig<<", Efficiency="<<Eff<<", Purity="<<Purity);
	#endif
	
	//Clear stuff
	if(gMVAReader) {
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Deleting MVA reader ...");
		#endif
		delete gMVAReader;
		gMVAReader= 0;
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("done");
		#endif
	}

	return 0;

}//close EvaluateCutClassifier()



int ApplyClassifier(std::string weightFileName)
{
	//######################################
	//##   INIT MVA READER
	//######################################
	//Initialize MVA reader
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initializing MVA reader ...");
	#endif
  if(!gMVAReader){
		gMVAReader= new TMVA::Reader("!Color:!Silent");
	}
	gMVAReader->AddVariable("peakSNR",&gPeakSNR);
	//gMVAReader->AddVariable("eccentricityRatio",&gEccentricityRatio);	
	gMVAReader->AddVariable("eccentricityPar",&gEccentricityPar);
	gMVAReader->AddVariable("areaRatio",&gAreaRatio);
	gMVAReader->AddSpectator("classOutput",&gClassOutput);	

	//Setting weights in reader
	TString weightFile= weightFileName.c_str();
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting weights from file "<<weightFile.Data()<<" ...");
	#endif
	gMVAReader->BookMVA(gClassifierMethod.c_str(),weightFile);
	

	//#########################################
	//##   APPLY CLASSIFICATION TO SOURCES
	//#########################################
	long int nBkgSources= 0;
	long int nSignalSources= 0;

	for (int i=0;i<gDataTree->GetEntries();i++) 
	{
		gDataTree->GetEntry(i);

		if(i%1000==0) cout<<"--> Reading "<<i+1<<"/"<<gDataTree->GetEntries()<<" event..."<<endl;

		
		//Retrieve source from collection
		Source* source= 0;
		if(gNestedSourceIndex==-1) source= m_sources[gSourceIndex];
		else source= m_sources[gSourceIndex]->GetNestedSource(gNestedSourceIndex);

		//cout<<"--> Reading source "<<source->GetName()<<" (sindex="<<gSourceIndex<<", nindex="<<gNestedSourceIndex<<", cindex="<<gFitComponentIndex<<"), nFitComponents="<<source->GetNFitComponents()<<", nSelFitComponents="<<source->GetNSelFitComponents()<<" ..."<<endl;


		//Skip if no fit info available
		bool hasFitInfo= source->HasFitInfo();
		if(!hasFitInfo || gFitComponentIndex<0) continue;

		//Check is source is classified as signal by classifier
		bool isClassifiedAsSignal= false;
		if(IsClassifiedAsSignal(isClassifiedAsSignal,gMVAReader,gClassifierMethod,gClassSignalOutCut)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to evaluate classifier response on source no. "<<i<<" (sindex="<<gSourceIndex<<", nestedId="<<gNestedSourceIndex<<", compId="<<gFitComponentIndex<<"), skip it...");
			#endif
			continue;
		}
			
		//Set source component flag
		int componentFlag= 0;
		if(isClassifiedAsSignal) {
			componentFlag= eReal;
			nSignalSources++;
		}
		else {
			componentFlag= eFake;
			nBkgSources++;
		}

		if(source->SetFitComponentFlag(gFitComponentIndex,componentFlag)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to set flag of fit component "<<gFitComponentIndex<<"!");
			#endif
			continue;
		}
	}//end loop source data

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<gDataTree->GetEntries()<<" sources processed: #"<<nBkgSources<<" false, #"<<nSignalSources<<" real");
	#endif

	//double mvaErr = reader->GetMVAError();

	return 0;

}//close ApplyClassifier()

int ReadData(std::string filelist)
{
	//Check filenames
	std::ifstream fileStream(filelist);	
	std::string line;
	if (fileStream.fail() || !fileStream.is_open()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filelist<<" for reading...");
		#endif
		return -1;
	}

	//Store filenames present in lists
	std::vector<std::string> fileNames;
	std::string filename= "";
	while (std::getline(fileStream, line)) {
  	std::istringstream iss(line);
    if (!(iss >> filename)) { 
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read line from file "<<filelist<<"!");
			#endif
			return -1; 
		}
    fileNames.push_back(filename);
	}//end file read
				
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<fileNames.size()<<" files present in list...");
	#endif

	//Read source data per each list
	for(size_t i=0;i<fileNames.size();i++){
		if(ReadSourceData(fileNames[i])<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read source data for file no. "<<i+1<<"!");
			#endif
			return -1;
		}
	}//end loop files

	return 0;

}//close ReadData()

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

	
	int sourceIndex= static_cast<int>(m_sourceComponentData.size());//Start from previous list size (in case of multiple source catalog files)

	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);
		int type= aSource->Type;
		
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source no. "<<i+1<<"/"<<sourceTree->GetEntries()<<"...");
		#endif


		//Select source by type?
		if(selectSourceByType){
			bool skipSource= true;
			for(size_t j=0;j<stypes.size();j++){
				if( stypes[j]==-1 || type==stypes[j]) {
					skipSource= false;
					break;
				}
			}
			if(skipSource) continue;
		}


		//Copy source
		Source* source= new Source;
		*source= *aSource;

		//Fill source pars
		if(FillSourcePars(m_sourceComponentData,source,sourceIndex,-1)<0){
			#ifdef LOGGING_ENABLED	
				ERROR_LOG("Failed to fill component data for source no. "<<i+1<<" (name="<<source->GetName()<<")!");
			#endif
			return -1;
		}


		//Add sources to list
		m_sources.push_back(source);
		sourceIndex++;

	}//end loop sources

	
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<m_sources.size()<<" sources read (#"<<m_sourceComponentData.size()<<" source component pars added) ...");
	#endif


	return 0;

}//close ReadSourceData()


int FillSourcePars(std::vector<SourceComponentData*>& pars,Source* aSource,int sourceIndex,int nestedSourceIndex)
{
	
	//Check input source
	if(!aSource) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Nullptr to input source given!");
		#endif
		return -1;	
	}

	//Get source data
	std::string sourceName= std::string(aSource->GetName());
	bool hasFitInfo= aSource->HasFitInfo();
	bool hasNestedSources= aSource->HasNestedSources();
	double beamArea= aSource->GetBeamFluxIntegral();
	double S= aSource->GetS();	
	
	double nPixels= static_cast<double>(aSource->NPix);
	double bkgRMSSum= aSource->GetBkgRMSSum();
	double rmsMean= bkgRMSSum/nPixels;
	if(rmsMean==0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No bkg info stored!");
		#endif
		return -1;
	}


	//## Fill source fit components
	SourceComponentData* cData= 0;

	if(hasFitInfo)
	{
		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= fitPars.GetNComponents();

		for(int k=0;k<nComponents;k++){
			//Skip source component if not selected
			bool isSelected= fitPars.IsSelectedComponent(k);
			//if(!isSelected) continue;

			std::string sname= sourceName + std::string(Form("_fitcomp%d",k+1));

			//Compute source component flag
			int flag= eUnknownSourceFlag;
			fitPars.GetComponentFlag(flag,k);

			int classFlag= -1;
			if(flag==eFake || flag==eCandidate) classFlag= 0;
			else if(flag==eReal) classFlag= 1;
			else{
				if(gTrainClassifier){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Unknown flag for source component "<<sname<<", will not use this in training...");
					#endif		
					continue;
				}	
			}
	
			//Compute peak signal to noise ratio
			double A= fitPars.GetParValue(k,"A");
			double peakSNR= A/rmsMean;

			//Compute eccentricity ratio
			double E= fitPars.GetComponentFitEllipseEccentricity(k);
			double E_beam= fitPars.GetComponentBeamEllipseEccentricity(k);
			double eccentricityRatio= 0;
			if(E_beam>0) eccentricityRatio= E/E_beam; 
			double eccentricityDiff= E-E_beam;

			//Compute area ratio
			double Area= fitPars.GetComponentFitEllipseArea(k);			
			double Area_beam= fitPars.GetComponentBeamEllipseArea(k);
			double areaRatio= 0;
			if(Area_beam>0) areaRatio= Area/Area_beam;

			//Check data integrity
			// - Check source SNR range
			if( peakSNR<gPeakSNRMin || peakSNR>gPeakSNRMax || TMath::IsNaN(peakSNR) || !std::isfinite(peakSNR) ){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Skip source component peakSNR (SNR="<<peakSNR<<") as NaN or outside selected range ("<<gPeakSNRMin<<","<<gPeakSNRMax<<")");
				#endif		
				continue;
			}

			//Check eccentricity par
			if(gUseEccentricityDiff)
			{
				//Check eccentricity diff
				if(eccentricityDiff<gEccentricityDiffMin || eccentricityDiff>gEccentricityDiffMax || TMath::IsNaN(eccentricityDiff) || !std::isfinite(eccentricityDiff) ) {
					#ifdef LOGGING_ENABLED
						WARN_LOG("Skip source component eccentricity diff for source "<<sname<<" (E_diff="<<eccentricityDiff<<") as NaN or outside selected range ("<<gEccentricityDiffMin<<","<<gEccentricityDiffMax<<")");
					#endif
					continue;
				}	
			}//close if
			else
			{
				if(eccentricityRatio<=0 || eccentricityRatio<gEccentricityRatioMin || eccentricityRatio>gEccentricityRatioMax || TMath::IsNaN(eccentricityRatio) || !std::isfinite(eccentricityRatio) ){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Skip source component eccentricity ratio for source "<<sname<<" (E_ratio="<<eccentricityRatio<<") as NaN or outside selected range ("<<gEccentricityRatioMin<<","<<gEccentricityRatioMax<<")");
					#endif
					continue;
				}
			}//close else

			
			//Check source to beam area ratio
			if(areaRatio<=0 || areaRatio<gSourceToBeamRatioMin || areaRatio>gSourceToBeamRatioMax || TMath::IsNaN(areaRatio) || !std::isfinite(areaRatio) ){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Skip source component beam area ratio (A_ratio="<<areaRatio<<") as NaN or outside selected range ("<<gSourceToBeamRatioMin<<","<<gSourceToBeamRatioMax<<")");
				#endif
				continue;
			}


			//Fill component par
			cData= new SourceComponentData(k,sourceIndex,nestedSourceIndex);
			cData->sname= sname;
			cData->flag= flag;		
			cData->classOutput= classFlag;
			cData->peakSNR= peakSNR;
			cData->eccentricityRatio= eccentricityRatio;
			cData->eccentricityDiff= eccentricityDiff;
			cData->areaRatio= areaRatio;
			cData->isSelected= isSelected;

			//Add component to source pars
			pars.push_back(cData);

		}//end loop components

	}//close if has fit info

	
	//## Fill pars for nested sources (if any)
	if(hasNestedSources){
		std::vector<Source*> nestedSources= aSource->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			if(FillSourcePars(pars,nestedSources[j],sourceIndex,j)<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
				#endif
				return -1;
			}
		}//end loop sources
	}//close if nested sources
	
	return 0;

}//close FillSourcePars()

int SaveDS9Regions()
{
	//Save DS9 regions for islands
	bool convertDS9RegionsToWCS= false;
	int status= SourceExporter::WriteToDS9(regionOutputFileName,m_sources,convertDS9RegionsToWCS);
	if(status<0){
		return -1;
	}

	//Save DS9 regions for components
	status= SourceExporter::WriteComponentsToDS9(regionComponentsOutputFileName,m_sources,convertDS9RegionsToWCS);
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
	for(size_t i=0;i<m_sources.size();i++){
		m_source= m_sources[i];
		outputTree->Fill();
	}

	outputTree->Write();

	return 0;

}//close SaveSources()


int SaveCatalog()
{
	//Return if no sources are found
	if(m_sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No sources selected, no catalog file will be written!");
		#endif
		return 0;
	}

	//Retrieve source WCS
	WCS* wcs= m_sources[0]->GetWCS(ds9WCSType);
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
	int status= SourceExporter::WriteToAscii(catalogOutputFileName,m_sources,dumpNestedSourceInfo,ds9WCSType,wcs);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source catalog to file "<<catalogOutputFileName<<" failed!");
		#endif
	}
	
	//Saving source fitted components to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source catalog to file "<<catalogComponentsOutputFileName<<" ...");
	#endif
	status= SourceExporter::WriteComponentsToAscii(catalogComponentsOutputFileName,m_sources,dumpNestedSourceInfo,ds9WCSType,wcs);
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

	//Save train data
	if(!gRunInteractively && outputFile_trainData && outputFile_trainData->IsOpen())
	{
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving train data to file ...");
		#endif
		outputFile_trainData->cd();
		if(gSignalDataTree_train) gSignalDataTree_train->Write();
		if(gSignalDataTree_test) gSignalDataTree_test->Write();
		if(gBkgDataTree_train) gBkgDataTree_train->Write();
		if(gBkgDataTree_test) gBkgDataTree_test->Write();
		if(ClassOutHisto_signal) ClassOutHisto_signal->Write();
		if(ClassOutHisto_bkg) ClassOutHisto_bkg->Write();	
		if(ClassOutPlot) ClassOutPlot->Write();
	}

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


void SetStyle()
{
	TStyle* myStyle= new TStyle("myStyle","myStyle");

	//## CANVAS & PAD
	myStyle->SetCanvasDefH(700); 
  myStyle->SetCanvasDefW(700); 
	myStyle->SetFrameBorderMode(0);
	myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
	myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetPadLeftMargin(0.16);
  myStyle->SetPadRightMargin(0.1);

	//## TITLE
	myStyle->SetOptTitle(0);
	myStyle->SetTitleX(0.1f);
	myStyle->SetTitleW(0.8f);                 
  myStyle->SetTitleXOffset(0.8);
  myStyle->SetTitleYOffset(1.1);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(0);//border size of Title PavelLabel
	myStyle->SetTitleSize(0.06,"X");
  myStyle->SetTitleSize(0.06,"Y");
	myStyle->SetTitleSize(0.06,"Z");
	
	//## STAT
	myStyle->SetOptStat("eMR");
	//myStyle->SetOptStat(1);
  myStyle->SetStatColor(0);
	myStyle->SetStatY(0.975);                
  myStyle->SetStatX(0.95);                
  myStyle->SetStatW(0.35);//0.2                
  myStyle->SetStatH(0.10);//0.15
  myStyle->SetStatBorderSize(1);

	myStyle->SetTitleFont(52,"X");
  myStyle->SetTitleFont(52,"Y");
  myStyle->SetTitleFont(52,"Z");
  myStyle->SetLabelFont(42,"X");
  myStyle->SetLabelFont(42,"Y");
  myStyle->SetLabelFont(42,"Z");   
	
	//## OTHER
  myStyle->SetOptFit(1);
	myStyle->SetOptLogx(0);
	myStyle->SetOptLogy(0);
  //myStyle->SetPalette(1,0);
  myStyle->SetMarkerStyle(8);
  myStyle->SetMarkerSize(0.6);
  myStyle->SetFuncWidth(1.); 
  myStyle->SetErrorX(0.);

	myStyle->SetNumberContours(999);
	//myStyle->SetPalette(kColorPrintableOnGrey);
	//myStyle->SetPalette(kTemperatureMap);
	//myStyle->SetPalette(kBlueGreenYellow);
	myStyle->SetPalette(kRainBow);
	
	gROOT->SetStyle("myStyle");
	gStyle= myStyle;
	myStyle->cd();

}//close SetStyle()

