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

void Usage(char* exeName)
{
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"*** Mandatory options ***"<<endl;
	cout<<"--input=[FILENAME] - Input file name containing image to be read (NB: .fits/.root supported)"<<endl;
	cout<<endl;
	cout<<"*** Optional options ***"<<endl;
  cout<<"-h, --help - Show help message and exit"<<endl;
	cout<<"--boxsize=[SIZE] - Size of sampling box in pixels or expressed as a multiple of the image beam size (if --sizeinbeam option is given) (default=100 pixels)"<<endl;
	cout<<"--sizeinbeam - Consider box size option expressed in multiple of beam size (beam info read from image) (default=no)"<<endl;
	cout<<"--gridsize=[SIZE] - Size of the interpolation grid expressed as fraction of the sampling box (default=0.25)"<<endl;
	cout<<"--estimator=[ESTIMATOR] - Bkg estimator used in the sampling box (1=mean, 2=median, 3=biweight, 4=clipped median) (default=2)"<<endl;
	cout<<"--2ndpass - If given, perform a 2nd pass in bkg calculation (default=no)"<<endl;
	cout<<"--skipblobs - If given, skip blobs using a flood-fill algorithm (default=no)"<<endl;
	cout<<"--seedthr=[NSIGMAS] - Seed threshold in flood-fill algorithm in nsigmas significance (default=5)"<<endl;
	cout<<"--mergethr=[NSIGMAS] - Merge threshold in flood-fill algorithm in nsigmas significance (default=2.6)"<<endl;
	cout<<"--minnpixels=[NPIX] - Minimum number of pixels in a blob (default=5)"<<endl;
	cout<<"--nthreads=[N] - Number of threads to be used for reading (-1=all available threads) (default=1)"<<endl;
	cout<<"--output=[FILENAME] - ROOT file where to save output maps (default=bkg.root)"<<endl;
	cout<<"--output-bkg=[FILENAME] - FITS file where to save bkg map (if --fitsout is given) (default=bkg.fits)"<<endl;	
	cout<<"--output-rms=[FILENAME] - FITS file where to save rms map (if --fitsout is given) (default=rms.fits)"<<endl;
	cout<<"--significance - Save the significance map (along with bkg and noise maps) in output file (default=no)"<<endl;
	cout<<"--output-significance=[FILENAME] - FITS file where to save significance map (if --fitsout is given) (default=significance.fits)"<<endl;
	cout<<"--fitsout - Write results in FITS files (default=no)"<<endl;
	cout<<"--imgname=[NAME] - Image name to be read in input ROOT file (if non standard) (default=img)"<<endl;
	cout<<"--parallel - Use parallel std algorithms for median (default=no)"<<endl;
	cout<<"-v [LEVEL], --verbosity=[LEVEL] - Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;

}//close Usage()

static const struct option options_tab[] = {
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "output", required_argument, 0, 'o' },
	{ "output-bkg", required_argument, 0, 'F' },
	{ "output-rms", required_argument, 0, 'R' },
	{ "significance", no_argument, 0, 's' },
	{ "output-significance", required_argument, 0, 'O' },
	{ "fitsout", no_argument, 0, 'f' },
	{ "gridsize", required_argument, 0, 'g' },
	{ "boxsize", required_argument, 0, 'b' },	
	{ "sizeinbeam", no_argument, 0, 'B' },
	{ "estimator", required_argument, 0, 'e' },
	{ "imgname", required_argument, 0, 'I' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "2ndpass", no_argument, 0, 'P'},
	{ "skipblobs", no_argument, 0, 'S'},
	{ "seedthr", required_argument, 0, 'T'},
	{ "mergethr", required_argument, 0, 't'},
	{ "minnpixels", required_argument, 0, 'm'},
	{ "nthreads", required_argument, 0, 'n' },
	{ "parallel", no_argument, 0, 'p' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
bool useDefaultOutput= true;
bool useDefaultOutput_bkg= true;
bool useDefaultOutput_rms= true;
bool useDefaultOutput_significance= true;
bool computeSignificance= false;
int verbosity= 4;//INFO level
bool use2ndPass= false;
bool skipOutliers= false;
double seedThr= 5;
double mergeThr= 2.6;
int minPixels= 5;
bool writeToFITS= false;
bool boxSizeInBeam= false;
double boxSize= 100;
double gridSize= 0.25;
int bkgEstimator= 2;
Caesar::BkgEstimator estimator;
int nthreads= 1;
bool useParallelVersion= false;

//Globar vars
TFile* outputFile= 0;
std::string inputFileName= "";
std::string outputFileName= "bkg.root";
std::string outputFileName_bkg= "bkg.fits";
std::string outputFileName_rms= "rms.fits";
std::string outputFileName_significance= "significance.fits";
std::string imageName= "img";
Caesar::Image* inputImg= 0;
Caesar::ImgBkgData* bkgData= 0;
Caesar::FileInfo info;

//Functions
int ParseOptions(int argc, char *argv[]);
int CheckOptions();
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
int main(int argc, char *argv[])
{
	//Start timer
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

	//=======================
	//== Read memory usage
	//=======================
	ProcMemInfo memInfo;
	if(SysUtils::GetProcMemoryInfo(memInfo)<0){
		ERROR_LOG("Failed to read process memory info!");		
		Clear();
		return -1;
	}	

	INFO_LOG("===========================");
	INFO_LOG("===   PERFORMANCE INFO  ===");
	INFO_LOG("===========================");
	INFO_LOG("dt(ms)= "<<dt);
	INFO_LOG("dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]");
	INFO_LOG("dt_stats(ms)= "<<dt_stats<<" ["<<dt_stats/dt*100.<<"%]");
	INFO_LOG("dt_bkg(ms)= "<<dt_bkg<<" ["<<dt_bkg/dt*100.<<"%]");
	INFO_LOG("dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]");
	INFO_LOG("real mem(kB)= "<<memInfo.realMem);
	INFO_LOG("real mem peak(kB)= "<<memInfo.realPeakMem);
	INFO_LOG("virt mem(kB)= "<<memInfo.virtMem);
	INFO_LOG("virt mem peak(kB)= "<<memInfo.virtPeakMem);
	INFO_LOG("===========================");
	
	
	INFO_LOG("End background finder");
	
	return 0;

}//close main


void Save(){

	//=======================
	//== Output 
	//=======================
	Caesar::Image* BkgMap= bkgData->BkgMap;
	Caesar::Image* NoiseMap= bkgData->NoiseMap;
	Caesar::Image* SignificanceMap= 0;
	if(computeSignificance) SignificanceMap= inputImg->GetSignificanceMap(bkgData,true);
	
	if(writeToFITS){
		BkgMap->WriteFITS(outputFileName_bkg);
		NoiseMap->WriteFITS(outputFileName_rms);
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
	while((c = getopt_long(argc, argv, "hfsi:o:O:F:R:b:Bg:e:I:v:PST:t:m:n:p",options_tab, &option_index)) != -1) 
	{
    
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
				if(!optarg){
					cerr<<"ERROR: Null string to input file argument given!"<<endl;
					exit(1);
				}
				inputFileName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				if(!optarg){
					cerr<<"ERROR: Null string to output file argument given!"<<endl;
					exit(1);
				}
				outputFileName= std::string(optarg);	
				useDefaultOutput= false;
				break;	
			}
			case 'O':	
			{
				if(!optarg){
					cerr<<"ERROR: Null string to significance output file argument given!"<<endl;
					exit(1);
				}
				outputFileName_significance= std::string(optarg);	
				useDefaultOutput_significance= false;
				break;	
			}
			case 'F':	
			{
				if(!optarg){
					cerr<<"ERROR: Null string to output bkg FITS file argument given!"<<endl;
					exit(1);
				}
				outputFileName_bkg= std::string(optarg);	
				useDefaultOutput_bkg= false;
				break;	
			}
			case 'R':	
			{
				if(!optarg){
					cerr<<"ERROR: Null string to output rms FITS file argument given!"<<endl;
					exit(1);
				}
				outputFileName_rms= std::string(optarg);	
				useDefaultOutput_rms= false;
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
			case 'B':
			{
				boxSizeInBeam= true;
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
				if(!optarg){
					cerr<<"ERROR: Null string to output image name argument given!"<<endl;
					exit(1);
				}
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
			case 'm':	
			{
				minPixels= atoi(optarg);
				break;	
			}
			case 'n':	
			{
				nthreads= atoi(optarg);	
				break;	
			}
			case 'p':	
			{
				useParallelVersion= true;				
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
	
	//## Check options
	if(CheckOptions()<0){
		ERROR_LOG("Invalid program options given, see logs!");
		return -1;
	}

	//## Print options
	INFO_LOG("========= OPTIONS ============");
	INFO_LOG("input file: "<<inputFileName);
	INFO_LOG("image name: "<<imageName);
	INFO_LOG("output file: "<<outputFileName);
	INFO_LOG("boxSize: "<<boxSize);
	INFO_LOG("gridSize: "<<gridSize);
	INFO_LOG("bkgEstimator: "<<bkgEstimator);	
	INFO_LOG("computeSignificance? "<<computeSignificance);
	INFO_LOG("use2ndPass? "<<use2ndPass);
	INFO_LOG("skipOutliers? "<<skipOutliers);
	INFO_LOG("seedThr: "<<seedThr);
	INFO_LOG("mergeThr: "<<mergeThr);
	INFO_LOG("minPixels: "<<minPixels);
	INFO_LOG("nthreads: "<<nthreads<<" (NThreads="<<SysUtils::GetOMPMaxThreads()<<")");
	INFO_LOG("useParallelAlgo? "<<useParallelVersion);	
	INFO_LOG("===============================");
	
	
	return 0;

}//close ParseOptions()

int CheckOptions()
{
	//Check mandatory options
	if(inputFileName==""){
		ERROR_LOG("Missing or empty input file argument!");
		return -1;
	}

	//Check optional options
	//- Output file option
	if(!useDefaultOutput && outputFileName==""){
		ERROR_LOG("Empty output file argument given!");
		return -1;
	}
	if(!useDefaultOutput_bkg && outputFileName_bkg==""){
		ERROR_LOG("Empty fits output file for bkg map argument given!");
		return -1;
	}
	if(!useDefaultOutput_rms && outputFileName_rms==""){
		ERROR_LOG("Empty fits output file for rms map argument given!");
		return -1;
	}	
	if(!useDefaultOutput_significance && outputFileName_significance==""){
		ERROR_LOG("Empty fits output file for significance map argument given!");
		return -1;
	}
	
	//- Detection options
	if(seedThr<=0 || mergeThr<=0 || minPixels<=0){
		ERROR_LOG("Invalid blob detection options given (hint: they must be positive)!");
		return -1;
	}

	//- Sampling box options
	if(boxSize<=0 ){
		ERROR_LOG("Invalid sampling box size option given (hint: must be positive)!");
		return -1;
	}
	if(gridSize<=0 || gridSize>1){
		ERROR_LOG("Invalid sampling grid size option given (hint: must be positive and smaller than 1)!");
		return -1;
	}

	//- Set bkg estimator
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
		ERROR_LOG("Invalid/unknown bkg estimator specified!");
		return -1;
	}

	//- Number of threads
	if(nthreads<=0){
		ERROR_LOG("Invalid number of threads given (hint: must be positive)!");
		return -1;
	}
	if(nthreads>SysUtils::GetOMPMaxThreads()){
		WARN_LOG("A number of threads exceeding detected machine core capacity ("<<SysUtils::GetOMPMaxThreads()<<") was given...");
	}
	SysUtils::SetOMPThreads(nthreads);

	//- Check given input file and get info
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,false)){
		ERROR_LOG("Invalid input file ("<<inputFileName<<") specified!");
		return -1;
	}
	std::string file_extension= info.extension;
	if(file_extension!= ".fits" && file_extension!=".root") {
		ERROR_LOG("Invalid file extension ("<<file_extension<<")...nothing to be done!");
		return -1;
	}

	return 0;

}//close CheckOptions()

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


int OpenOutputFile()
{
	// Set output filenames
	if(!writeToFITS){//ROOT output
		outputFile= new TFile(outputFileName.c_str(),"RECREATE");
	}

	/*
	if(useDefaultOutput){
		std::string basefilename_wext= info.filename_wext;
		INFO_LOG("basefilename_wext="<<basefilename_wext);

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
	*/

	return 0;

}//close OpenOutputFile()

int ReadImage()
{
	//--> ROOT reading
	if(info.extension==".root"){// Read image from ROOT file
		INFO_LOG("Reading ROOT input file "<<inputFileName<<"...");
		TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<inputFileName<<"!");
			return -1;
		}
		inputImg=  (Caesar::Image*)inputFile->Get(imageName.c_str());
		if(!inputImg){
			ERROR_LOG("Cannot get image from input file "<<inputFileName<<"!");
			return -1;
		}
	}//close if

	//--> FITS reading
	if(info.extension==".fits"){// Read image from FITS file
		INFO_LOG("Reading FITS input file "<<inputFileName<<"...");
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
	bool computeRobustStats= true;
	bool useRange= false;
	bool forceRecomputing= false;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	if(inputImg->ComputeStats(computeRobustStats,forceRecomputing,useRange,minThr,maxThr,useParallelVersion)<0){
		ERROR_LOG("Stats computing failed!");
		return -1;
	}
	inputImg->PrintStats();	

	return 0;

}//close ComputeStats()

int ComputeBkg()
{
	//## Compute background 
	INFO_LOG("Starting background finder ...");
	
	//Compute sampling box size in pixels 
	//If box size is given as a beam multiple find first the beam size from image
	double boxSize_pix= boxSize;
	
	if(boxSizeInBeam){
		//Check if image has metadata with beam info
		if(!inputImg->HasMetaData()){
			ERROR_LOG("Requested to use sampling box size as multiple of image beam but input image has no metadata with beam information stored!");
			return -1;
		}
	
		//Compute box size in pixels
		int pixelWidthInBeam= inputImg->GetMetaData()->GetBeamWidthInPixel();	
		INFO_LOG("Read a beam size of "<<pixelWidthInBeam<<" pixels from input image ...");

		boxSize_pix= pixelWidthInBeam*boxSize;
		
	}//close if

	//Compute grid size in pixels
	double gridSize_pix= gridSize*boxSize_pix;

	INFO_LOG("Computing background assuming sampling box of size "<<boxSize_pix<<" pixels and a grid size of "<<gridSize_pix<<" pixels ...");
	
	//Check grid & box size
	long int Nx= inputImg->GetNx();
	long int Ny= inputImg->GetNy();
	if(boxSize_pix>=min(Nx,Ny) || gridSize_pix>=min(Nx,Ny) ){
		ERROR_LOG("Box/grid size are too large compared to image size ("<<Nx<<","<<Ny<<")");
		return -1;
	}

	//Compute bkg & noise maps
	bool computeLocalBkg= true;
	bkgData= inputImg->ComputeBkg(estimator,computeLocalBkg,boxSize_pix,boxSize_pix,gridSize_pix,gridSize_pix,use2ndPass,skipOutliers,seedThr,mergeThr,minPixels);
	if(!bkgData){
		ERROR_LOG("Failed to compute bkg data!");
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

