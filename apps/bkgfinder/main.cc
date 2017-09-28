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

//Caesar headers
#include <Image.h>
#include <Logger.h>
#include <SysUtils.h>
#include <FITSReader.h>
#include <BkgData.h>
#include <BkgFinder.h>

//ROOT headers
#include <TFile.h>

//C++ headers
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
	cout<<"-i, --input \t Input file name containing image to be read (only .fits/.root supported)"<<endl;
	cout<<"-b, --boxsize \t Size of sampling box in pixels (often a multiple of the image beam size)"<<endl;
	cout<<"-g, --gridsize \t Granularity size of the interpolation grid in pixels (i.e. half the box size)"<<endl;
	cout<<"-e, --estimator \t Bkg estimator used in the sampling box (1=mean, 2=median, 3=biweight, 4=clipped median)"<<endl;
	cout<<"-f, --fitsout \t Write results in FITS files"<<endl;
	cout<<"-o, --output \t Output file name (1 file for ROOT out and multiple files if FITS out option is selected)"<<endl;
	cout<<"-I, --imgname \t Image name in input ROOT file (if non standard)"<<endl;
	cout<<"-s, --significance \t Compute and store also the significance map (along with bkg and noise maps)"<<endl;
	cout<<"-P, --2ndpass \t If given, perform a 2nd pass in bkg calculation (default=no)"<<endl;
	cout<<"-S, --skipblobs \t If given, skip blobs using a flood-fill algorithm (default=no)"<<endl;
	cout<<"-T, --seedthr \t Seed threshold in flood-fill algorithm in nsigmas significance (typical value=5)"<<endl;
	cout<<"-t, --mergethr \t Merge threshold in flood-fill algorithm in nsigmas significance (typical value=2.6)"<<endl;
	cout<<"-v, --verbosity \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "output", optional_argument, 0, 'o' },
	{ "fitsout", no_argument, 0, 'f' },
	{ "gridsize", required_argument, 0, 'g' },
	{ "boxsize", required_argument, 0, 'b' },
	{ "estimator", required_argument, 0, 'e' },
	{ "significance", no_argument, 0, 's' },
	{ "imgname", optional_argument, 0, 'I' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "2ndpass", no_argument, 0, 'P'},
	{ "skipblobs", no_argument, 0, 'S'},
	{ "seedthr", required_argument, 0, 'T'},
	{ "mergethr", required_argument, 0, 't'},
	{ "minnpixels", required_argument, 0, 'n'},
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
bool useDefaultOutput= true;
bool computeSignificance= false;
int verbosity= 4;//INFO level
bool use2ndPass= false;
bool skipOutliers= false;
double seedThr= 5;
double mergeThr= 2.6;
int minPixels= 10;
bool writeToFITS= false;
double boxSize= 0;
double gridSize= 0;
int bkgEstimator= 0;
Caesar::BkgEstimator estimator;

//Globar vars
TFile* outputFile= 0;
std::string inputFileName= "";
std::string outputFileName= "";
std::string outputFileName_bkg= "";
std::string outputFileName_noise= "";
std::string outputFileName_significance= "";
std::string imageName= "img";
//Caesar::Img* inputImg= 0;
Caesar::Image* inputImg= 0;
//Caesar::BkgData* bkgData= 0;
Caesar::ImgBkgData* bkgData= 0;
Caesar::FileInfo info;


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int OpenOutputFile();
int ReadImage();
int ComputeStats();
int ComputeBkg();
void Clear();
void Save();

//=======================================
//===             MAIN               ====
//=======================================
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
		ERROR_LOG("Failed to compute stats & bkg!");		
		Clear();
		return -1;
	}
	auto t1_bkg = chrono::steady_clock::now();
	double dt_bkg= chrono::duration <double, milli> (t1_bkg-t0_bkg).count();

	//=======================
	//== Save to file
	//=======================
	auto t0_save = chrono::steady_clock::now();	
	Save();
	auto t1_save = chrono::steady_clock::now();
	double dt_save= chrono::duration <double, milli> (t1_save-t0_save).count();

	//=======================
	//== Clear data
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
	INFO_LOG("dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]");
	INFO_LOG("===========================");
	
	
	INFO_LOG("End background finder");
	
	return 0;

}//close main


void Save(){

	//=======================
	//== Output 
	//=======================
	//Caesar::Img* BkgMap= bkgData->BkgMap;
	//Caesar::Img* NoiseMap= bkgData->NoiseMap;
	//Caesar::Img* SignificanceMap= 0;
	Caesar::Image* BkgMap= bkgData->BkgMap;
	Caesar::Image* NoiseMap= bkgData->NoiseMap;
	Caesar::Image* SignificanceMap= 0;
	if(computeSignificance) SignificanceMap= inputImg->GetSignificanceMap(bkgData,true);
	
	if(writeToFITS){
		BkgMap->WriteFITS(outputFileName_bkg);
		NoiseMap->WriteFITS(outputFileName_noise);
		if(computeSignificance && SignificanceMap) SignificanceMap->WriteFITS(outputFileName_significance);
	}
	else{
		if(outputFile && outputFile->IsOpen()){
			outputFile->cd();
			BkgMap->SetName("bkgMap");
			BkgMap->Write();		
			NoiseMap->SetName("rmsMap");
			NoiseMap->Write();
			if(computeSignificance && SignificanceMap) {	
				SignificanceMap->SetName("significanceMap");
				SignificanceMap->Write();
			}
			outputFile->Close();
		}
	}

}//close Save()

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
	while((c = getopt_long(argc, argv, "hfsi:o::b:g:e:I::v:PST:t:n:",options_tab, &option_index)) != -1) {
    
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
				inputFileName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				useDefaultOutput= false;
				break;	
			}
			case 'f':
			{
      	writeToFITS= true;
				break;
			}
			case 'b':
			{
				boxSize= atof(optarg);
				break;
			}
			case 'g':
			{
				gridSize= atof(optarg);
				break;
			}
			case 'e':
			{
				bkgEstimator= atoi(optarg);
				break;
			}
			case 's':	
			{
				computeSignificance= true;
				break;
			}
			case 'I':	
			{
				imageName= std::string(optarg);	
				break;	
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			case 'P':	
			{
				use2ndPass= true;
				break;	
			}
			case 'S':	
			{
				skipOutliers= true;
				break;	
			}
			case 'T':	
			{
				seedThr= atof(optarg);
				break;	
			}
			case 't':	
			{
				mergeThr= atof(optarg);
				break;	
			}	
			case 'n':	
			{
				minPixels= atoi(optarg);
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while

	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	
	//## Set bkg estimator
	estimator= Caesar::eMedianBkg;
	if(bkgEstimator==Caesar::eMedianBkg){
		estimator= Caesar::eMedianBkg;
	}
	else if(bkgEstimator==Caesar::eMeanBkg){
		estimator= Caesar::eMeanBkg;
	}
	else if(bkgEstimator==Caesar::eBiWeightBkg){
		estimator= Caesar::eBiWeightBkg;
	}	
	else if(bkgEstimator==Caesar::eMedianClippedBkg){
		estimator= Caesar::eMedianClippedBkg;
	}
	else{
		ERROR_LOG("Invalid bkg estimator specified!");
		return -1;
	}
	
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


int OpenOutputFile(){

	// Set output filenames
	if(useDefaultOutput){
		std::string basefilename_wext= info.filename_wext;
		if(writeToFITS){//FITS output
			std::string outputFileNamePrefix= basefilename_wext;
			outputFileName_bkg= outputFileNamePrefix + std::string("bkg.fits");
			outputFileName_noise= outputFileNamePrefix + std::string("rms.fits");	
			outputFileName_significance= outputFileNamePrefix + std::string("significance.fits");
		}
		else{//ROOT output
			outputFileName= basefilename_wext + std::string("bkg.root");
			outputFile= new TFile(outputFileName.c_str(),"RECREATE");
		}
	}//close if
	else{
		if(writeToFITS){//FITS output
			std::string outputFileNamePrefix= outputFileName;
			outputFileName_bkg= outputFileNamePrefix + std::string("bkg.fits");
			outputFileName_noise= outputFileNamePrefix + std::string("rms.fits");	
			outputFileName_significance= outputFileNamePrefix + std::string("significance.fits");
		}
		else{
			outputFile= new TFile(outputFileName.c_str(),"RECREATE");
		}
	}//close else

	return 0;

}//close OpenOutputFile()

int ReadImage(){

	// Check given input file and get info
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

	//--> ROOT reading
	if(file_extension==".root"){// Read image from ROOT file
		INFO_LOG("Reading ROOT input file "<<inputFileName<<"...");
		TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<inputFileName<<"!");
			return -1;
		}
		//inputImg=  (Caesar::Img*)inputFile->Get(imageName.c_str());
		inputImg=  (Caesar::Image*)inputFile->Get(imageName.c_str());
		if(!inputImg){
			ERROR_LOG("Cannot get image from input file "<<inputFileName<<"!");
			return -1;
		}
	}//close if

	//--> FITS reading
	if(file_extension==".fits"){// Read image from FITS file
		INFO_LOG("Reading FITS input file "<<inputFileName<<"...");
		//inputImg= new Caesar::Img; 
		//inputImg->SetNameTitle(imageName.c_str(),imageName.c_str());

		inputImg= new Caesar::Image;
		inputImg->SetName(imageName);
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

	//## Compute background 
	INFO_LOG("Starting background finder ...");
	
	//Check grid & box size
	//int Nx= inputImg->GetNbinsX();
	//int Ny= inputImg->GetNbinsY();
	long int Nx= inputImg->GetNx();
	long int Ny= inputImg->GetNy();
	if(boxSize>=min(Nx,Ny) || gridSize>=min(Nx,Ny) ){
		ERROR_LOG("Box/grid size are too large compared to image size ("<<Nx<<","<<Ny<<")");
		//inputImg->Delete();
		return -1;
	}

	//Compute bkg & noise maps
	bool computeLocalBkg= true;
	
	//bkgData= inputImg->ComputeBkg(estimator,true,boxSize,boxSize,gridSize,gridSize,use2ndPass,skipOutliers,seedThr,mergeThr,minPixels);
	bkgData= inputImg->ComputeBkg(estimator,computeLocalBkg,boxSize,boxSize,gridSize,gridSize,use2ndPass,skipOutliers,seedThr,mergeThr,minPixels);
	if(!bkgData){
		ERROR_LOG("Failed to compute bkg data!");
		//inputImg->Delete();
		return -1;
	}

	return 0;

}//close ComputeBkg()


void Clear(){

	//Clear input image
	if(inputImg) inputImg->Delete();

	//Clear bkg data
	if(bkgData){ 
		delete bkgData;
		bkgData= 0;
	}

}//close Clear()

