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
#include <SysUtils.h>
#include <FITSReader.h>
#include <BkgData.h>
#include <BkgFinder.h>
#include <ConfigParser.h>
#include <Consts.h>
#include <Logger.h>


#include <TFile.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <chrono>

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
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
  {(char*)0, (int)0, (int*)0, (int)0}
};



//Options
//--> Main options
//Img* inputImg= 0;
Image* inputImg= 0;
TFile* inputFile= 0;
std::string inputFileName= "";
std::string imageName= "";
TFile* outputFile= 0;	
std::string outputFileName;
//Img* saliencyImg= 0;	
Image* saliencyImg= 0;	
bool saveToFile;
bool saveConfig;
bool saveResidualMap;
bool saveBkgMap;
bool saveNoiseMap;
bool saveSignificanceMap;
bool saveInputMap;
bool saveSaliencyMap;
bool saveSources;
bool saveToFITSFile;
std::string residualMapFITSFile;
std::string inputMapFITSFile;
std::string saliencyMapFITSFile;
std::string bkgMapFITSFile;
std::string noiseMapFITSFile;
std::string significanceMapFITSFile;

//--> bkg options
//BkgData* bkgData= 0;
ImgBkgData* bkgData= 0;
//Img* significanceMap= 0;
Image* significanceMap= 0;
double boxSizeX, boxSizeY;
double gridSizeX, gridSizeY;
int bkgEstimator;
bool useLocalBkg;

//Functions
int ParseOptions(int argc, char *argv[]);
int ReadImage();
int ComputeStats();
int ComputeBkg();
int ApplySmoothing();
int OpenOutputFile();
int FindSaliency();
int Clear();
int Save();


int main(int argc, char *argv[]){

	
	auto t0 = chrono::steady_clock::now();
	
	//================================
	//== Parse command line options
	//================================
	auto t0_parse = chrono::steady_clock::now();
	if(ParseOptions(argc,argv)<0){
		ERROR_LOG("Failed to parse command line options!");
		Clear();
		return -1;
	}
	auto t1_parse = chrono::steady_clock::now();
	double dt_parse= chrono::duration <double, milli> (t1_parse-t0_parse).count();

	
	//=======================
	//== Open out file
	//=======================		
	auto t0_outfile = chrono::steady_clock::now();
	if(OpenOutputFile()<0){
		ERROR_LOG("Failed to open output file!");
		Clear();
		return -1;	
	}
	auto t1_outfile = chrono::steady_clock::now();
	double dt_outfile= chrono::duration <double, milli> (t1_outfile-t0_outfile).count();

	//=======================
	//== Read image
	//=======================	
	auto t0_read = chrono::steady_clock::now();
	if(ReadImage()<0){
		ERROR_LOG("Failed to read image from file!");
		Clear();
		return -1;
	}
	auto t1_read = chrono::steady_clock::now();
	double dt_read= chrono::duration <double, milli> (t1_read-t0_read).count();

	
	//=======================
	//== Apply smoothing
	//=======================
	auto t0_smooth = chrono::steady_clock::now();	
	if(ApplySmoothing()<0){
		ERROR_LOG("Failed to perform image smoothing!");
		Clear();
		return -1;
	}
	auto t1_smooth = chrono::steady_clock::now();
	double dt_smooth= chrono::duration <double, milli> (t1_smooth-t0_smooth).count();

	//=======================
	//== Compute stats
	//=======================
	auto t0_stats = chrono::steady_clock::now();
	if(ComputeStats()<0){
		ERROR_LOG("Failed to read image from file!");		
		Clear();
		return -1;
	}
	auto t1_stats = chrono::steady_clock::now();
	double dt_stats= chrono::duration <double, milli> (t1_stats-t0_stats).count();


	//=======================
	//== Background finder
	//=======================
	auto t0_bkg = chrono::steady_clock::now();	
	if(ComputeBkg()<0){
		ERROR_LOG("Failed to compute bkg!");
		Clear();
		return -1;
	}
	auto t1_bkg = chrono::steady_clock::now();
	double dt_bkg= chrono::duration <double, milli> (t1_bkg-t0_bkg).count();


	//=======================
	//== Saliency filtering
	//=======================
	auto t0_sal = chrono::steady_clock::now();	
	if(FindSaliency()<0){
		ERROR_LOG("Failed to compute saliency map!");
		Clear();
		return -1;
	}		
	auto t1_sal = chrono::steady_clock::now();
	double dt_sal= chrono::duration <double, milli> (t1_sal-t0_sal).count();

	//=======================
	//== Save results 
	//=======================
	auto t0_save = chrono::steady_clock::now();	
	Save();
	auto t1_save = chrono::steady_clock::now();
	double dt_save= chrono::duration <double, milli> (t1_save-t0_save).count();


	//=======================
	//== Clear
	//=======================
	Clear();
	
	//=======================
	//== Print perf stats
	//=======================
	auto t1 = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (t1-t0).count();

	INFO_LOG("===========================");
	INFO_LOG("===   PERFORMANCE INFO  ===");
	INFO_LOG("===========================");
	INFO_LOG("dt(ms)= "<<dt);
	INFO_LOG("dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]");
	INFO_LOG("dt_stats(ms)= "<<dt_stats<<" ["<<dt_stats/dt*100.<<"%]");
	INFO_LOG("dt_bkg(ms)= "<<dt_bkg<<" ["<<dt_bkg/dt*100.<<"%]");
	INFO_LOG("dt_smooth(ms)= "<<dt_smooth<<" ["<<dt_smooth/dt*100.<<"%]");
	INFO_LOG("dt_saliency(ms)= "<<dt_sal<<" ["<<dt_sal/dt*100.<<"%]");
	INFO_LOG("dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]");
	INFO_LOG("===========================");
	
	

	INFO_LOG("End saliency finder");
	
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


int FindSaliency(){
	
	//## Get options
	int saliencyResoMin= 0;
	int saliencyResoMax = 0; 
	int saliencyResoStep= 0;
	if(GET_OPTION_VALUE(saliencyResoMin,saliencyResoMin)<0){
		ERROR_LOG("Failed to get saliencyResoMin option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyResoMax,saliencyResoMax)<0){
		ERROR_LOG("Failed to get saliencyResoMax option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyResoStep,saliencyResoStep)<0){
		ERROR_LOG("Failed to get saliencyResoStep option!");
		return -1;
	}

	double spBeta;
	int spMinArea;
	if(GET_OPTION_VALUE(spBeta,spBeta)<0){
		ERROR_LOG("Failed to get spBeta option!");
		return -1;
	}
	if(GET_OPTION_VALUE(spMinArea,spMinArea)<0){
		ERROR_LOG("Failed to get spMinArea option!");
		return -1;
	}

	double saliencyNNFactor;
	double saliencySpatialRegFactor;
	bool saliencyUseRobustPars;
	double saliencyMultiResoCombThrFactor;
	double saliencyDissExpFalloffPar;
	double saliencySpatialDistRegPar;
	if(GET_OPTION_VALUE(saliencyNNFactor,saliencyNNFactor)<0){
		ERROR_LOG("Failed to get saliencyNNFactor option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyUseRobustPars,saliencyUseRobustPars)<0){
		ERROR_LOG("Failed to get saliencyUseRobustPars option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyMultiResoCombThrFactor,saliencyMultiResoCombThrFactor)<0){
		ERROR_LOG("Failed to get saliencyMultiResoCombThrFactor option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencySpatialRegFactor,saliencySpatialRegFactor)<0){
		ERROR_LOG("Failed to get saliencySpatialRegFactor option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyDissExpFalloffPar,saliencyDissExpFalloffPar)<0){
		ERROR_LOG("Failed to get saliencyDissExpFalloffPar option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencySpatialDistRegPar,saliencySpatialDistRegPar)<0){
		ERROR_LOG("Failed to get saliencySpatialDistRegPar option!");
		return -1;
	}

	bool saliencyUseBkgMap, saliencyUseNoiseMap;
	if(GET_OPTION_VALUE(saliencyUseBkgMap,saliencyUseBkgMap)<0){
		ERROR_LOG("Failed to get saliencyUseBkgMap option!");
		return -1;
	} 
	if(GET_OPTION_VALUE(saliencyUseNoiseMap,saliencyUseNoiseMap)<0){
		ERROR_LOG("Failed to get saliencyUseNoiseMap option!");
		return -1;
	} 

	double saliencyThrFactor, saliencyImgThrFactor;
	if(GET_OPTION_VALUE(saliencyThrFactor,saliencyThrFactor)<0){
		ERROR_LOG("Failed to get saliencyThrFactor option!");
		return -1;
	} 
	if(GET_OPTION_VALUE(saliencyImgThrFactor,saliencyImgThrFactor)<0){
		ERROR_LOG("Failed to get saliencyImgThrFactor option!");
		return -1;
	}
		
	/*
	bool saliencyUseCurvInDiss;	
	if(GET_OPTION_VALUE(saliencyUseCurvInDiss,saliencyUseCurvInDiss)<0){
		ERROR_LOG("Failed to get saliencyUseCurvInDiss option!");
		return -1;
	}
	*/
		
	INFO_LOG("saliencyReso("<<saliencyResoMin<<","<<saliencyResoMax<<","<<saliencyResoStep<<") spBeta="<<spBeta<<". spMinArea="<<spMinArea<<", saliencyNNFactor="<<saliencyNNFactor<<",  saliencyUseRobustPars="<<saliencyUseRobustPars<<" saliencyDissExpFalloffPar="<<saliencyDissExpFalloffPar<<", saliencySpatialDistRegPar="<<saliencySpatialDistRegPar<<", saliencyMultiResoCombThrFactor="<<saliencyMultiResoCombThrFactor<<" saliencyUseBkgMap="<<saliencyUseBkgMap<<" saliencyUseNoiseMap="<<saliencyUseNoiseMap<<" saliencyThrFactor="<<saliencyThrFactor<<", saliencyImgThrFactor="<<saliencyImgThrFactor);

	//## Compute saliency
	saliencyImg= inputImg->GetMultiResoSaliencyMap(
		saliencyResoMin,saliencyResoMax,saliencyResoStep,
		spBeta,spMinArea,saliencyNNFactor,saliencyUseRobustPars,saliencyDissExpFalloffPar,saliencySpatialDistRegPar, 
		saliencyMultiResoCombThrFactor,
  	saliencyUseBkgMap,saliencyUseNoiseMap,bkgData,
		saliencyThrFactor,saliencyImgThrFactor
	);

	if(!saliencyImg){
		ERROR_LOG("Failed to compute saliency map!");
		return -1;
	}

	return 0;

}//FindSaliency()


int ApplySmoothing(){
	
	//## Get options
	bool usePreSmoothing;
	if(GET_OPTION_VALUE(usePreSmoothing,usePreSmoothing)<0){
		ERROR_LOG("Failed to get usePreSmoothing option!");
		return -1;
	}
	if(!usePreSmoothing){
		INFO_LOG("No pre-smoothing selected");
		return 0;
	}

	int smoothFilter;
	int gausFilterKernSize;
	double gausFilterSigma;
	double guidedFilterRadius, guidedFilterColorEps;
	if(GET_OPTION_VALUE(smoothFilter,smoothFilter)<0){
		ERROR_LOG("Failed to get smoothFilter option!");
		return -1;
	}
	if(GET_OPTION_VALUE(gausFilterKernSize,gausFilterKernSize)<0){
		ERROR_LOG("Failed to get gausFilterKernSize option!");
		return -1;
	}
	if(GET_OPTION_VALUE(gausFilterSigma,gausFilterSigma)<0){
		ERROR_LOG("Failed to get gausFilterSigma option!");
		return -1;
	}
	if(GET_OPTION_VALUE(guidedFilterRadius,guidedFilterRadius)<0){
		ERROR_LOG("Failed to get guidedFilterRadius option!");
		return -1;
	}
	if(GET_OPTION_VALUE(guidedFilterColorEps,guidedFilterColorEps)<0){
		ERROR_LOG("Failed to get guidedFilterColorEps option!");
		return -1;
	}	

	//## Apply a smoothing stage?
	//Img* smoothedImg= 0;
	Image* smoothedImg= 0;
	if(smoothFilter==eGaus){
		smoothedImg= inputImg->GetSmoothedImage(gausFilterKernSize,gausFilterKernSize,gausFilterSigma,gausFilterSigma);
	}
	else if(smoothFilter==eGuided){
		smoothedImg= inputImg->GetGuidedFilterImage(guidedFilterRadius,guidedFilterColorEps);
	}
	else{
		ERROR_LOG("Invalid smoothing algo selected!");
		return -1;
	}

	if(!smoothedImg){
		ERROR_LOG("Computation of smoothed map failed!");
		return -1;
	}

	// Compute stats
	/*
	INFO_LOG("Computing input image stats...");
	if(smoothedImg->ComputeStats(true,false,false)<0){
		ERROR_LOG("Stats computing failed!");
		smoothedImg->Delete();
		return -1;
	}
	*/
	
	//## Replace input image with smoothed map
	TString imgName= inputImg->GetName();
	inputImg->Delete();
	smoothedImg->SetNameTitle(imgName,imgName);
	inputImg= smoothedImg;
	
	return 0;

}//close ApplySmoothing()

int ComputeStats(){

	//## Compute stats
	INFO_LOG("Computing input image stats...");
	if(inputImg->ComputeStats(true,false,false)<0){
		ERROR_LOG("Stats computing failed!");
		return -1;
	}
	inputImg->PrintStats();	

	return 0;

}//close ComputeStats()


int ComputeBkg(){

	//## Get bkg options
	//Beam info
	bool useBeamInfoInBkg;
	if(GET_OPTION_VALUE(useBeamInfoInBkg,useBeamInfoInBkg)<0){
		ERROR_LOG("Failed to get useBeamInfoInBkg option!");
		return -1;
	}
	int nPixelsInBeam= 0;
	if(useBeamInfoInBkg && inputImg->HasMetaData()){
		nPixelsInBeam= inputImg->GetMetaData()->GetBeamWidthInPixel();	
	}
		
	//Box size
	if(GET_OPTION_VALUE(boxSizeX,boxSizeX)<0 || GET_OPTION_VALUE(boxSizeY,boxSizeY)<0){
		ERROR_LOG("Failed to get boxSize option!");
		return -1;
	}

	//Grid size
	if(GET_OPTION_VALUE(gridSizeX,gridSizeX)<0 || GET_OPTION_VALUE(gridSizeY,gridSizeY)<0){
		ERROR_LOG("Failed to get gridSize option!");
		return -1;
	}
	
	//Bkg estimator
	if(GET_OPTION_VALUE(bkgEstimator,bkgEstimator)<0){
		ERROR_LOG("Failed to get bkgEstimator option!");
		return -1;
	}

	//Local bkg flag
	bool use2ndPassInLocalBkg;
	bool skipOutliersInLocalBkg;
	int minNPix;
	//double seedBrightThr;
	double seedThr;
	double dilateZThr;
	double mergeThr;

	if(GET_OPTION_VALUE(useLocalBkg,useLocalBkg)<0){
		ERROR_LOG("Failed to get useLocalBkg option!");
		return -1;
	}

	if(GET_OPTION_VALUE(use2ndPassInLocalBkg,use2ndPassInLocalBkg)<0){
		ERROR_LOG("Failed to get use2ndPassInLocalBkg option!");
		return -1;
	}
	
	if(GET_OPTION_VALUE(skipOutliersInLocalBkg,skipOutliersInLocalBkg)<0){
		ERROR_LOG("Failed to get skipOutliersInLocalBkg option!");
		return -1;
	}
	
	if(GET_OPTION_VALUE(minNPix,minNPix)<0){
		ERROR_LOG("Failed to get minNPix option!");
		return -1;
	}
	
	//if(GET_OPTION_VALUE(seedBrightThr,seedBrightThr)<0){
	//	ERROR_LOG("Failed to get seedBrightThr option!");
	//	return -1;
	//}
	
	if(GET_OPTION_VALUE(seedThr,seedThr)<0){
		ERROR_LOG("Failed to get seedThr option!");
		return -1;
	}

	if(GET_OPTION_VALUE(dilateZThr,dilateZThr)<0){
		ERROR_LOG("Failed to get dilateZThr option!");
		return -1;
	}
	
	if(GET_OPTION_VALUE(mergeThr,mergeThr)<0){
		ERROR_LOG("Failed to get mergeThr option!");
		return -1;
	}


	//Compute box size
	//double Nx= static_cast<double>(inputImg->GetNbinsX());
	//double Ny= static_cast<double>(inputImg->GetNbinsY());
	double Nx= static_cast<double>(inputImg->GetNx());
	double Ny= static_cast<double>(inputImg->GetNy());
	if(useBeamInfoInBkg && nPixelsInBeam>0){
		INFO_LOG("Setting bkg boxes as ("<<boxSizeX<<","<<boxSizeY<<") x beam (beam="<<nPixelsInBeam<<" pixels) ...");
		boxSizeX*= nPixelsInBeam;
		boxSizeY*= nPixelsInBeam;
	}
	else{
		WARN_LOG("Beam information is not available or its usage has been turned off, using image fractions...");
		
		boxSizeX*= Nx;
		boxSizeY*= Ny;
	}
	INFO_LOG("Setting bkg boxes to ("<<boxSizeX<<","<<boxSizeY<<") pixels ...");	

	gridSizeX*= boxSizeX;
	gridSizeY*= boxSizeY;
	INFO_LOG("Setting grid size to ("<<gridSizeX<<","<<gridSizeY<<") pixels ...");
	
	//## Check grid & box size
	if(boxSizeX>=Nx || boxSizeY>=Ny || gridSizeX>=Nx || gridSizeY>=Ny){
		ERROR_LOG("Box/grid size are too large compared to image size ("<<Nx<<","<<Ny<<")");
		return -1;
	}

	//Compute bkg & noise maps
	bkgData= inputImg->ComputeBkg(
		bkgEstimator,
		useLocalBkg,boxSizeX,boxSizeY,gridSizeX,gridSizeY,
		use2ndPassInLocalBkg,
		skipOutliersInLocalBkg,seedThr,mergeThr,minNPix
	);
	if(!bkgData){
		ERROR_LOG("Failed to compute bkg data!");
		return -1;
	}

	//Compute significance
	significanceMap= inputImg->GetSignificanceMap(bkgData,useLocalBkg);
	if(!significanceMap){
		ERROR_LOG("Failed to compute significance map!");
		return -1;
	}

	return 0;

}//close ComputeBkg()


int ReadImage(){

	//## Get options
	if(GET_OPTION_VALUE(inputFile,inputFileName)<0){
		ERROR_LOG("Failed to get inputFile option!");
		return -1;
	}	
	if(GET_OPTION_VALUE(inputImage,imageName)<0){
		ERROR_LOG("Failed to get inputImage option!");
		return -1;
	}
	
	//## Check given input file and get info
	Caesar::FileInfo info;
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,false)){
		ERROR_LOG("Invalid input file ("<<inputFileName<<") specified!");
		return -1;
	}
	std::string file_extension= info.extension;
	if(file_extension!= ".fits" && file_extension!=".root") {
		ERROR_LOG("Invalid file extension ("<<file_extension<<")...nothing to be done!");
		return -1;
	}

	//## Read image
	//===== ROOT reading =====
	if(file_extension==".root"){// Read image from ROOT file
		INFO_LOG("Reading ROOT input file "<<inputFileName<<" (image name="<<imageName<<")...");
		inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<inputFileName<<"!");
			return -1;
		}
		//inputImg=  (Img*)inputFile->Get(imageName.c_str());	
		inputImg=  (Image*)inputFile->Get(imageName.c_str());
		if(!inputImg){
			ERROR_LOG("Cannot get image from input file "<<inputFileName<<"!");
			return -1;
		}
	}//close if

	//===== FITS reading =====
	if(file_extension==".fits"){// Read image from FITS file
		INFO_LOG("Reading FITS input file "<<inputFileName<<"...");
		//inputImg= new Caesar::Img; 
		inputImg= new Caesar::Image; 	
		inputImg->SetNameTitle(imageName.c_str(),imageName.c_str());
		if(inputImg->ReadFITS(inputFileName)<0){
			ERROR_LOG("Failed to read image from input file "<<inputFileName<<"!");	
			return -1;
		}
	}//close else if

	if(!inputImg){
		ERROR_LOG("Failed to read image from input file "<<inputFileName<<"!");
		return -1;
	}

	return 0;

}//close ReadImage()

int OpenOutputFile(){

	//Get options
	if(GET_OPTION_VALUE(outputFile,outputFileName)<0){
		ERROR_LOG("Failed to get outputFile option!");
		return -1;
	}	
	if(GET_OPTION_VALUE(saveToFile,saveToFile)<0){
		ERROR_LOG("Failed to get saveToFile option!");
		return -1;
	}	
	if(GET_OPTION_VALUE(saveConfig,saveConfig)<0){
		ERROR_LOG("Failed to get saveConfig option!");
		return -1;
	}	
	if(GET_OPTION_VALUE(saveResidualMap,saveResidualMap)<0){
		ERROR_LOG("Failed to get saveResidualMap option!");
		return -1;
	}	
	if(GET_OPTION_VALUE(saveBkgMap,saveBkgMap)<0){
		ERROR_LOG("Failed to get saveBkgMap option!");
		return -1;
	}	
	if(GET_OPTION_VALUE(saveNoiseMap,saveNoiseMap)<0){
		ERROR_LOG("Failed to get saveNoiseMap option!");
		return -1;
	}		
	if(GET_OPTION_VALUE(saveSignificanceMap,saveSignificanceMap)<0){
		ERROR_LOG("Failed to get saveSignificanceMap option!");
		return -1;
	}	
	if(GET_OPTION_VALUE(saveInputMap,saveInputMap)<0){
		ERROR_LOG("Failed to get saveInputMap option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saveSaliencyMap,saveSaliencyMap)<0){
		ERROR_LOG("Failed to get saveSaliencyMap option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saveSources,saveSources)<0){
		ERROR_LOG("Failed to get saveSources option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saveToFITSFile,saveToFITSFile)<0){
		ERROR_LOG("Failed to get saveToFITSFile option!");
		return -1;
	}

	if(GET_OPTION_VALUE(residualMapFITSFile,residualMapFITSFile)<0){
		ERROR_LOG("Failed to get residualMapFITSFile option!");
		return -1;
	}
	if(GET_OPTION_VALUE(inputMapFITSFile,inputMapFITSFile)<0){
		ERROR_LOG("Failed to get inputMapFITSFile option!");
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyMapFITSFile,saliencyMapFITSFile)<0){
		ERROR_LOG("Failed to get saliencyMapFITSFile option!");
		return -1;
	}
	if(GET_OPTION_VALUE(bkgMapFITSFile,bkgMapFITSFile)<0){
		ERROR_LOG("Failed to get bkgMapFITSFile option!");
		return -1;
	}
	if(GET_OPTION_VALUE(noiseMapFITSFile,noiseMapFITSFile)<0){
		ERROR_LOG("Failed to get noiseMapFITSFile option!");
		return -1;
	}
	if(GET_OPTION_VALUE(significanceMapFITSFile,significanceMapFITSFile)<0){
		ERROR_LOG("Failed to get significanceMapFITSFile option!");
		return -1;
	}

	if(saveToFile){
		outputFile= new TFile(outputFileName.c_str(),"RECREATE");
		if(!outputFile || !outputFile->IsOpen()){
			ERROR_LOG("Failed to open output file!");
			return -1;
		}
	}
	
	return 0;

}//close OpenOutputFile()


int Clear(){

	//Clear input image
	if(inputImg) inputImg->Delete();

	//Clear bkg data
	if(bkgData){ 
		delete bkgData;
		bkgData= 0;
	}

	//Delete significance map
	if(significanceMap) significanceMap->Delete();

	//Clear saliency image
	if(saliencyImg) saliencyImg->Delete();

	return 0;

}//close Clear()

int Save(){

	//## Save to ROOT?
	if(saveToFile && outputFile){
		outputFile->cd();
		if(saveConfig){
			TTree* configTree= ConfigParser::Instance().GetConfigTree();
			if(configTree) configTree->Write();
		}
		if(saveInputMap && inputImg) {
			inputImg->SetName("img");
			inputImg->Write();
		}
		if(saveBkgMap && bkgData->BkgMap) {
			(bkgData->BkgMap)->SetName("img_bkg");
			(bkgData->BkgMap)->Write();
		}
		if(saveNoiseMap && bkgData->NoiseMap) {
			(bkgData->NoiseMap)->SetName("img_rms");
			(bkgData->NoiseMap)->Write();
		}
		if(saveSignificanceMap && significanceMap) {		
			significanceMap->SetName("img_significance");
			significanceMap->Write();
		}
		if(saveSaliencyMap && saliencyImg) {
			saliencyImg->SetName("img_saliency");
			saliencyImg->Write();
		}
		outputFile->Close();
	}

	//## Save to FITS?
	if(saveToFITSFile){
		if(saveInputMap && inputImg) inputImg->WriteFITS(inputMapFITSFile);
		if(saveBkgMap && bkgData->BkgMap) {
			(bkgData->BkgMap)->WriteFITS(bkgMapFITSFile);
		}
		if(saveNoiseMap && bkgData->NoiseMap) {
			(bkgData->NoiseMap)->WriteFITS(noiseMapFITSFile);
		}
		if(saveSignificanceMap && significanceMap) {	
			significanceMap->WriteFITS(significanceMapFITSFile);
		}
		if(saveSaliencyMap && saliencyImg) saliencyImg->WriteFITS(saliencyMapFITSFile);
	}

	return 0;
}//close Save()

