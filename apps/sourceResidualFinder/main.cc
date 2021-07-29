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
#include <Contour.h>
#include <Source.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <Consts.h>


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
	cout<<"-c, --config=[FILENAME] \t Config file containing option settings"<<endl;	
	cout<<"-i, --inputfile=[FILENAME] \t Filename (fits/root) with input image."<<endl;
	cout<<"-s, --sourcefile=[FILENAME] \t Caesar ROOT file with source list. If provided no sources will be searched."<<endl;
	cout<<"-o, --outputfile=[FILENAME] \t Filename where to store output residual image (default=resmap.fits)"<<endl;
	cout<<"-O, --outputfile_mask=[FILENAME] \t Filename where to store output source mask image (default=smask.fits)"<<endl;
	cout<<"-p, --psSubtractionMethod=[METHOD] - Method used to remove point sources (1=DILATION,2=MODEL SUBTRACTION) (default=1)"<<endl;
	cout<<"-r, --resZThr=[NSIGMAS] - Significance threshold (in sigmas) above which sources are removed (if selected for removal) (default=5)"<<endl;
	cout<<"-R, --resZHighThr=[NSIGMAS] - Significance threshold (in sigmas) above which sources are always removed (even if they have nested or different type) (default=10)"<<endl;
	cout<<"-l, --removedSourceType=[TYPE] - Type of bright sources to be dilated from the input image (-1=ALL,1=COMPACT,2=POINT-LIKE,3=EXTENDED)"<<endl;
	cout<<"-a, --removeNestedSources - If a source has nested sources, remove nested rather than mother source (default=no)"<<endl;
	cout<<"-k, --dilateKernelSize=[SIZE] - Kernel size in pixel used to dilate image around sources (default=9)"<<endl;
	cout<<"-V, --globalbkg - Use global bkg rather than local bkg (default=use local)"<<endl;
	cout<<"-z, --bkgAroundSource - Use bkg computed in a box around source and not from the bkg map (default=use bkg map)"<<endl;
	cout<<"-Z, --bkgBoxThickness=[THICKNESS] - Bkg box thickness in pixels (default=20)"<<endl;
	cout<<"-j, --randomizeBkg - Randomize bkg in dilated pixels (default=no)"<<endl;
	cout<<"-T, --seedthr=[NSIGMAS] - Seed threshold in flood-fill algorithm in nsigmas significance (default=5)"<<endl;
	cout<<"-t, --mergethr=[NSIGMAS] - Merge threshold in flood-fill algorithm in nsigmas significance (default=2.6)"<<endl;
	cout<<"-m, --minnpixels=[NPIX] - Minimum number of pixels in a blob (default=5)"<<endl;
	cout<<"-N, --no-nested - Do not search nested sources (default=search)"<<endl;
	cout<<"-b, --boxsize=[SIZE] - Size of sampling box in pixels or expressed as a multiple of the image beam size (if --sizeinbeam option is given) (default=100 pixels)"<<endl;
	cout<<"-B, --sizeinbeam - Consider box size option expressed in multiple of beam size (beam info read from image) (default=no)"<<endl;
	cout<<"-g, --gridsize=[SIZE] - Size of the interpolation grid expressed as fraction of the sampling box (default=0.25)"<<endl;
	cout<<"-e, --estimator=[ESTIMATOR] - Bkg estimator used in the sampling box (1=mean, 2=median, 3=biweight, 4=clipped median) (default=2)"<<endl;
	cout<<"-P, --2ndpass - If given, perform a 2nd pass in bkg calculation (default=no)"<<endl;
	cout<<"-S, --skipblobs - If given, skip blobs using a flood-fill algorithm (default=no)"<<endl;
	cout<<"-A, --no-selection - Do not select and retag input sources (default=apply selection)"<<endl;
	cout<<"-d, --no-maxnpixcut - Do not apply max n pixel cut (default=apply)"<<endl;
	cout<<"-D, --maxnpix=[MAX_NPIX] - max number of pixel to consider a source as point-like (default=1000)"<<endl;
	cout<<"-f, --no-elongcut - Do not apply elongation cut (default=apply)"<<endl;
	cout<<"-G, --no-circratiocut - Do not apply circular ratio cut (default=apply)"<<endl;
	cout<<"-q, --no-ellipsearearatiocut - Do not apply ellipse area ratio cut (default=apply)"<<endl;
	cout<<"-u, --no-nbeamscut - Do not apply nbeams cut (default=apply)"<<endl;
	cout<<"-U, --maxnbeams=[MAX_NBEAMS] - Max number of beams in source to consider it as point-like (default=10)"<<endl;
	cout<<"-n, --nthreads \t Number of threads to be used (default=1)"<<endl;
	cout<<"-v [LEVEL], --verbosity=[LEVEL] - Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
	{ "inputfile", required_argument, 0, 'i' },
	{ "outputfile", required_argument, 0, 'o' },
	{ "outputfile_mask", required_argument, 0, 'O' },
	{ "sourcefile", required_argument, 0, 's' },
	{ "psSubtractionMethod", required_argument, 0, 'p' },
	{ "resZThr", required_argument, 0, 'r' },
	{ "resZHighThr", required_argument, 0, 'R' },
	{ "removedSourceType", required_argument, 0, 'l' },
	{ "removeNestedSources", no_argument, 0, 'a' },
	{ "dilateKernelSize", required_argument, 0, 'k' },
	{ "globalbkg", no_argument, 0, 'V' },
	{ "bkgAroundSource", no_argument, 0, 'z' },
	{ "bkgBoxThickness", required_argument, 0, 'Z' },
	{ "randomizeBkg", no_argument, 0, 'j' },
	{ "seedthr", required_argument, 0, 'T'},
	{ "mergethr", required_argument, 0, 't'},
	{ "minnpixels", required_argument, 0, 'm'},
	{ "no-nested", no_argument, 0, 'N'},
	{ "gridsize", required_argument, 0, 'g' },
	{ "boxsize", required_argument, 0, 'b' },	
	{ "sizeinbeam", no_argument, 0, 'B' },
	{ "estimator", required_argument, 0, 'e' },
	{ "2ndpass", no_argument, 0, 'P'},
	{ "skipblobs", no_argument, 0, 'S'},
	{ "no-selection", no_argument, 0, 'A'},
	{ "no-maxnpixcut", no_argument, 0, 'd'},
	{ "maxnpix", required_argument, 0, 'D'},
	{ "no-elongcut", no_argument, 0, 'f'},
	{ "no-circratiocut", no_argument, 0, 'G'},
	{ "no-ellipsearearatiocut", no_argument, 0, 'q'},
	{ "no-nbeamscut", no_argument, 0, 'u'},
	{ "maxnbeams", required_argument, 0, 'U'},
	{ "nthreads", required_argument, 0, 'n' },
	{ "verbosity", required_argument, 0, 'v'},
  {(char*)0, (int)0, (int*)0, (int)0}
};



//Options
//--> Main options
int verbosity= 4;//INFO level
int nThreads= 1;
Image* inputImg= 0;
TFile* inputFile= 0;
std::string inputFileName= "";
std::string imageName= "";
std::string sourceFileName= "";
bool findSources= true;
double seedThr= 5;
double mergeThr= 2.6;
int minNPix= 5;
bool searchNestedSources= true;
bool bkgAroundSource= false;
int bkgBoxThickness= 20;

TFile* outputFile= 0;	
std::string outputFileName= "resmap.fits";
Image* residualImg= 0;
std::string outputFileName_mask= "smask.fits";
std::string outputFileName_bkg= "bkg.fits";
Image* smaskImg= 0;
Image* smaskImg_binary= 0;
std::vector<Source*> sources;	


//--> Source residual options
bool removeNestedSources= false;
int dilateKernelSize= 9;
int removedSourceType= 2;
int residualModel= 1;
double residualZHighThr= 10;
double residualZThr= 5;
bool residualModelRandomize= false;
int psSubtractionMethod= 1;

//--> Source selection
bool applySourceSelection= true;
double sourceMinBoundingBox= 2;
bool useCircRatioCut= true;
double psCircRatioThr= 0.4;
bool useElongCut= true;
double psElongThr= 0.7;
bool useEllipseAreaRatioCut= true;
double psEllipseAreaRatioMinThr= 0.6;
double psEllipseAreaRatioMaxThr= 1.4;
bool useMaxNPixCut= true;
int psMaxNPix= 1000;
bool useNBeamsCut= true;
double psNBeamsThr= 10;

//--> Bkg options
ImgBkgData* bkgData= 0;
Image* significanceMap= 0;
bool useLocalBkg= true;
bool boxSizeInBeam= false;
double boxSize= 100;
double boxSizeX= boxSize;
double boxSizeY= boxSize;
double gridSize= 0.25;
double gridSizeX= gridSize;
double gridSizeY= gridSize;
int bkgEstimator= 2;
Caesar::BkgEstimator estimator;
bool use2ndPass= false;
bool skipOutliers= false;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verb);
int SetOptionsFromConfig();
int CheckOptions();
int ReadImage();
int ComputeStats(Image*);
int ComputeBkg(Image*);
//int OpenOutputFile();
int FindSources();
int SelectSources();
int ReadSources();
bool IsGoodSource(Source* aSource);
bool IsPointLikeSource(Source* aSource);
int ComputeSourceResidual();
int Clear();
int Save();

int main(int argc, char *argv[]){

	auto t0 = chrono::steady_clock::now();
	
	//================================
	//== PARSE OPTIONS
	//================================
	auto t0_parse = chrono::steady_clock::now();
	if(ParseOptions(argc,argv)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to parse command line options!");
		#endif
		Clear();
		return -1;
	}
	auto t1_parse = chrono::steady_clock::now();
	double dt_parse= chrono::duration <double, milli> (t1_parse-t0_parse).count();

	//=======================
	//== READ INPUT IMAGE
	//=======================
	auto t0_read = chrono::steady_clock::now();
	if(ReadImage()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read image from file!");		
		#endif
		Clear();
		return -1;
	}
	auto t1_read = chrono::steady_clock::now();
	double dt_read= chrono::duration <double, milli> (t1_read-t0_read).count();
	

	//=======================
	//== EXTRACT SOURCES
	//=======================
	double dt_stats= 0;
	double dt_bkg= 0;
	double dt_sfinder= 0;
	double dt_ssel= 0;

	if(findSources)
	{
		//=======================
		//== COMPUTE STATS
		//=======================
		auto t0_stats = chrono::steady_clock::now();
		if(ComputeStats(inputImg)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute stats!");		
			#endif
			Clear();
			return -1;
		}
		auto t1_stats = chrono::steady_clock::now();
		dt_stats= chrono::duration <double, milli> (t1_stats-t0_stats).count();

		//=======================
		//== FIND BACKGROUND
		//=======================
		auto t0_bkg = chrono::steady_clock::now();	
		if(ComputeBkg(inputImg)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute bkg!");
			#endif
			Clear();
			return -1;
		}
		auto t1_bkg = chrono::steady_clock::now();
		dt_bkg= chrono::duration <double, milli> (t1_bkg-t0_bkg).count();

	
		//=======================
		//== FIND SOURCES
		//=======================
		auto t0_sfinder = chrono::steady_clock::now();	
		if(FindSources()<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to find sources!");
			#endif
			Clear();
			return -1;
		}		
		auto t1_sfinder = chrono::steady_clock::now();
		dt_sfinder= chrono::duration <double, milli> (t1_sfinder-t0_sfinder).count();

		
	}//close if findSources 
	else
	{
		//===========================
		//== READ SOURCES FROM FILE
		//===========================	
		auto t0_sfinder = chrono::steady_clock::now();	
		if(ReadSources()<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read sources from file "<<inputFileName<<"!");	
			#endif
			Clear();
			return -1;
		}
		auto t1_sfinder = chrono::steady_clock::now();
		dt_sfinder= chrono::duration <double, milli> (t1_sfinder-t0_sfinder).count();

	}//close else


	//=======================
	//== SELECT SOURCES
	//=======================
	auto t0_ssel = chrono::steady_clock::now();	
	if(applySourceSelection){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Selecting sources ...");
		#endif
		if(SelectSources()<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to select sources ...");
			#endif
			Clear();
			return -1;
		}
	}
	auto t1_ssel = chrono::steady_clock::now();
	dt_ssel= chrono::duration <double, milli> (t1_ssel-t0_ssel).count();

	
	//====================================
	//== RECOMPUTE BKG EXCLUDING SOURCES
	//====================================
	//- Compute source masks
	#ifdef LOGGING_ENABLED
		INFO_LOG("Mask source pixels in original image ...");	
	#endif
	smaskImg= (Image*)inputImg->GetCloned("smaskImg");
	if(!smaskImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to clone input map!");	
		#endif
		Clear();
		return -1;
	}
	if(smaskImg->MaskSources(sources)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to mask sources in input map!");	
		#endif
		Clear();
		return -1;
	}


	bool isBinary= true;
	bool invertMask= true;
	smaskImg_binary= inputImg->GetSourceMask(sources,isBinary,invertMask);
	if(!smaskImg_binary){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute source binary mask!");	
		#endif
		Clear();
		return -1;
	}


	//Compute stats for source mask
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing source mask stats ...");	
	#endif
	if(ComputeStats(smaskImg)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute stats for masked image!");		
		#endif
		Clear();
		return -1;
	}	

	//- Compute bkg for source mask, first clear existing bkg data
	if(bkgData){ 
		delete bkgData;
		bkgData= 0;
	}
	if(significanceMap) significanceMap->Delete();

	if(ComputeBkg(smaskImg)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute bkg for masked image!");
		#endif
		Clear();
		return -1;
	}


	//==========================
	//== COMPUTE RESIDUAL MAP
	//==========================
	auto t0_res = chrono::steady_clock::now();
	if(ComputeSourceResidual()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute source residual map!");
		#endif
		Clear();
		return -1;
	}
	auto t1_res = chrono::steady_clock::now();
	double dt_res= chrono::duration <double, milli> (t1_res-t0_res).count();

	//=======================
	//== SAVE TO FILE
	//=======================
	auto t0_save = chrono::steady_clock::now();	
	Save();
	auto t1_save = chrono::steady_clock::now();
	double dt_save= chrono::duration <double, milli> (t1_save-t0_save).count();

	//=======================
	//== CLEAR DATA
	//=======================
	Clear();

	//=======================
	//== PRINT STATS
	//=======================
	auto t1 = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (t1-t0).count();

	#ifdef LOGGING_ENABLED
		INFO_LOG("===========================");
		INFO_LOG("===   PERFORMANCE INFO  ===");
		INFO_LOG("===========================");
		INFO_LOG("dt(ms)= "<<dt);
		INFO_LOG("dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]");
		INFO_LOG("dt_stats(ms)= "<<dt_stats<<" ["<<dt_stats/dt*100.<<"%]");
		INFO_LOG("dt_bkg(ms)= "<<dt_bkg<<" ["<<dt_bkg/dt*100.<<"%]");
		INFO_LOG("dt_sfinder(ms)= "<<dt_sfinder<<" ["<<dt_sfinder/dt*100.<<"%]");
		INFO_LOG("dt_ssel(ms)= "<<dt_ssel<<" ["<<dt_ssel/dt*100.<<"%]");
		INFO_LOG("dt_res(ms)= "<<dt_res<<" ["<<dt_res/dt*100.<<"%]");
		INFO_LOG("dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]");
		INFO_LOG("===========================");
	#else
		cout<<"==========================="<<endl;
		cout<<"===   PERFORMANCE INFO  ==="<<endl;
		cout<<"==========================="<<endl;
		cout<<"dt(ms)= "<<dt<<endl;
		cout<<"dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]"<<endl;
		cout<<"dt_stats(ms)= "<<dt_stats<<" ["<<dt_stats/dt*100.<<"%]"<<endl;
		cout<<"dt_bkg(ms)= "<<dt_bkg<<" ["<<dt_bkg/dt*100.<<"%]"<<endl;
		cout<<"dt_sfinder(ms)= "<<dt_sfinder<<" ["<<dt_sfinder/dt*100.<<"%]"<<endl;
		cout<<"dt_ssel(ms)= "<<dt_ssel<<" ["<<dt_ssel/dt*100.<<"%]"<<endl;
		cout<<"dt_res(ms)= "<<dt_res<<" ["<<dt_res/dt*100.<<"%]"<<endl;
		cout<<"dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]"<<endl;
		cout<<"==========================="<<endl;
	#endif
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("End residual computation");
	#endif

	return 0;

}//close main

int ParseOptions(int argc, char *argv[])
{
	//- Check args
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//- Parse options
	std::string configFileName= "";
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hc:i:s:o:O:p:r:R:l:ak:T:t:m:n:Nb:Bg:e:PSAdD:fGquU:zZ:jv:V",options_tab, &option_index)) != -1) {
    
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
				if(!optarg){
					cerr<<"ERROR: Null string to config file argument given!"<<endl;
					exit(1);
				}
				configFileName= std::string(optarg);	
				break;	
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
			case 's':	
			{
				if(!optarg){
					cerr<<"ERROR: Null string to source file argument given!"<<endl;
					exit(1);
				}
				sourceFileName= std::string(optarg);	
				if(sourceFileName!="") findSources= false;
				break;	
			}
			case 'o':	
			{
				if(!optarg){
					cerr<<"ERROR: Null string to output file argument given!"<<endl;
					exit(1);
				}
				outputFileName= std::string(optarg);	
				break;	
			}		
			case 'O':	
			{
				if(!optarg){
					cerr<<"ERROR: Null string to output file mask argument given!"<<endl;
					exit(1);
				}
				outputFileName_mask= std::string(optarg);	
				break;	
			}
			case 'p':	
			{
				psSubtractionMethod= atoi(optarg);
				break;	
			}	
			case 'r':	
			{
				residualZThr= atof(optarg);
				break;	
			}
			case 'R':	
			{
				residualZHighThr= atof(optarg);
				break;	
			}
			case 'l':	
			{
				removedSourceType= atoi(optarg);
				break;	
			}
			case 'a':
			{
				removeNestedSources= true;		
				break;
			}
			case 'k':
			{
				dilateKernelSize= atoi(optarg);		
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
				minNPix= atoi(optarg);
				break;	
			}
			case 'N':	
			{
				searchNestedSources= false;
				break;	
			}
			case 'V':
			{
				useLocalBkg= false;
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
			case 'A':	
			{
				applySourceSelection= false;
				break;	
			}
			case 'd':
			{
				useMaxNPixCut= false;
				break;
			}
			case 'D':
			{
				psMaxNPix= atoi(optarg);
				break;
			}
			case 'f':
			{
				useElongCut= false;
				break;
			}
			case 'G':
			{
				useCircRatioCut= false;
				break;
			}
			case 'q':
			{
				useEllipseAreaRatioCut= false;
				break;
			}
			case 'u':
			{
				useNBeamsCut= false;
				break;
			}
			case 'U':
			{
				psNBeamsThr= atof(optarg);
				break;
			}
			case 'z':
			{
				bkgAroundSource= true;
				break;
			}
			case 'Z':
			{
				bkgBoxThickness= atoi(optarg);
				break;
			}
			case 'j':
			{
				residualModelRandomize= true;
				break;
			}
			case 'n':	
			{
				nThreads= atoi(optarg);	
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
 
	//- Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	#ifdef LOGGING_ENABLED
		LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	#endif

	//- Read config options (if provided) and set options
	if(configFileName!=""){
		if(ConfigParser::Instance().Parse(configFileName)<0){
			cerr<<"ERROR: Failed to parse config options!"<<endl;
			return -1;
		}
		PRINT_OPTIONS();

		SetOptionsFromConfig();
	}

	//- Check options
	if(CheckOptions()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid program options given, see logs!");
		#endif
		return -1;
	}

	return 0;

}//close ParseOptions()



int CheckOptions()
{
	//Check mandatory options
	if(inputFileName==""){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Missing or empty input file argument!");
		#endif
		return -1;
	}

	//- Detection options
	if(seedThr<=0 || mergeThr<=0 || minNPix<=0){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid blob detection options given (hint: they must be positive)!");
		#endif
		return -1;
	}

	//- Sampling box options
	if(boxSize<=0 ){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid sampling box size option given (hint: must be positive)!");
		#endif
		return -1;
	}
	if(gridSize<=0 || gridSize>1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid sampling grid size option given (hint: must be positive and smaller than 1)!");
		#endif
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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unknown bkg estimator specified!");
		#endif
		return -1;
	}

	//- Number of threads
	if(nThreads<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of threads given (hint: must be positive)!");
		#endif
		return -1;
	}
	if(nThreads>SysUtils::GetOMPMaxThreads()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("A number of threads exceeding detected machine core capacity ("<<SysUtils::GetOMPMaxThreads()<<") was given...");
		#endif
	}
	SysUtils::SetOMPThreads(nThreads);

	//- Check given input file and get info
	Caesar::FileInfo info;
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,false)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid input file ("<<inputFileName<<") specified!");
		#endif
		return -1;
	}
	std::string file_extension= info.extension;
	if(file_extension!= ".fits" && file_extension!=".root") {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid file extension ("<<file_extension<<")...nothing to be done!");
		#endif
		return -1;
	}

	return 0;

}//close CheckOptions()


int SetOptionsFromConfig()
{
	//=======================
	//== INIT LOGGER
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
	//== INIT THREAD NO.
	//=======================
	int nThreads;
	if(GET_OPTION_VALUE(nThreads,nThreads)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get nThreads option!");
		#endif
		return -1;
	}
	if(nThreads>0) SysUtils::SetOMPThreads(nThreads);


	//============================
	//== INPUT FILE OPTIONS
	//============================
	if(GET_OPTION_VALUE(inputFile,inputFileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get inputFile option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(inputImage,imageName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get inputImage option!");
		#endif
		return -1;
	}

	//============================
	//== SOURCE RESIDUAL OPTIONS
	//============================
	if(GET_OPTION_VALUE(removeNestedSources,removeNestedSources)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get removeNestedSources option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(dilateKernelSize,dilateKernelSize)<0){		
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get dilateKernelSize option!");
		#endif
		return -1;
	}	
	
	if(GET_OPTION_VALUE(removedSourceType,removedSourceType)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get removedSourceType option!");
		#endif
		return -1;
	}
	
	if(GET_OPTION_VALUE(residualModel,residualModel)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get residualModel option!");
		#endif
		return -1;
	}
	
	if(GET_OPTION_VALUE(residualModelRandomize,residualModelRandomize)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get residualModelRandomize option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(residualZThr,residualZThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get residualZThr option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(residualZHighThr,residualZHighThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get residualZHighThr option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(psSubtractionMethod,psSubtractionMethod)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psSubtractionMethod option!");
		#endif
		return -1;
	}	

	//============================
	//== SOURCE FIND OPTIONS
	//============================
	if(GET_OPTION_VALUE(seedThr,seedThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get seedThr option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(mergeThr,mergeThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get mergeThr option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(minNPix,minNPix)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get minNPix option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(searchNestedSources,searchNestedSources)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get searchNestedSources option!");
		#endif
		return -1;
	}

	//============================
	//== SOURCE SELECTION OPTIONS
	//============================
	if(GET_OPTION_VALUE(applySourceSelection,applySourceSelection)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get applySourceSelection option!");
		#endif
		return -1;
	}
	
	if(GET_OPTION_VALUE(sourceMinBoundingBox,sourceMinBoundingBox)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get sourceMinBoundingBox option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(useCircRatioCut,useCircRatioCut)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get useCircRatioCut option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(psCircRatioThr,psCircRatioThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psCircRatioThr option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(useElongCut,useElongCut)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get useElongCut option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(psElongThr,psElongThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psElongThr option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(useEllipseAreaRatioCut,useEllipseAreaRatioCut)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get useEllipseAreaRatioCut option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(psEllipseAreaRatioMinThr,psEllipseAreaRatioMinThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psEllipseAreaRatioMinThr option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(psEllipseAreaRatioMaxThr,psEllipseAreaRatioMaxThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psEllipseAreaRatioMaxThr option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(useMaxNPixCut,useMaxNPixCut)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psMaxNPix option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(psMaxNPix,psMaxNPix)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psMaxNPix option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(useNBeamsCut,useNBeamsCut)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get useNBeamsCut option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(psNBeamsThr,psNBeamsThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get psNBeamsThr option!");
		#endif
		return -1;
	}

	//============================
	//== BACKGROUND OPTIONS
	//============================
	
	//Box size
	if(GET_OPTION_VALUE(boxSizeX,boxSizeX)<0 || GET_OPTION_VALUE(boxSizeY,boxSizeY)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get boxSize option!");
		#endif
		return -1;
	}

	//Grid size
	if(GET_OPTION_VALUE(gridSizeX,gridSizeX)<0 || GET_OPTION_VALUE(gridSizeY,gridSizeY)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get gridSize option!");
		#endif
		return -1;
	}
	
	//Bkg estimator
	if(GET_OPTION_VALUE(bkgEstimator,bkgEstimator)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get bkgEstimator option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(useLocalBkg,useLocalBkg)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get useLocalBkg option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(use2ndPassInLocalBkg,use2ndPass)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get use2ndPassInLocalBkg option!");
		#endif
		return -1;
	}
	
	if(GET_OPTION_VALUE(skipOutliersInLocalBkg,skipOutliers)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get skipOutliersInLocalBkg option!");
		#endif
		return -1;
	}

	return 0;

}//close SetOptionsFromConfig()


int ComputeSourceResidual()
{
	//- Compute residual
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing source residual map ...");
	#endif
	ImgBkgData* bkgData_curr= nullptr;
	if(!bkgAroundSource) bkgData_curr= bkgData;

	residualImg= inputImg->GetSourceResidual(
		sources,
		dilateKernelSize,
		residualModel,
		removedSourceType,
		removeNestedSources,
		bkgData_curr,useLocalBkg,
		residualModelRandomize,
		residualZThr,residualZHighThr,
		psSubtractionMethod,
		smaskImg_binary,bkgBoxThickness
	);
	if(!residualImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute residual map!");
		#endif
		return -1;
	}

	return 0;

}//close ComputeSourceResidual()


int ReadSources()
{
	//- Read source file
	TFile* sourceFile= new TFile(sourceFileName.c_str(),"READ");
	if(!sourceFile || sourceFile->IsZombie()) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Cannot open file "<<sourceFileName<<"!");
		#endif
		return -1;
	}
	
	//- Get source tree
	TTree* SourceInfo= (TTree*)sourceFile->Get("SourceInfo");	
	if(!SourceInfo){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Cannot get sources from file "<<sourceFileName<<"!");
		#endif
		return -1;
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);

	cout<<"# "<<SourceInfo->GetEntries()<<" sources to be drawn..."<<endl;
	
	//- Store source list 
	sources.clear();
	
	for(int i=0;i<SourceInfo->GetEntries();i++){
		SourceInfo->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);

	}//end loop sources

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sources.size()<<" sources read ...");
	#endif

	return 0;

}//close ReadSources()

int FindSources()
{	
	//- Extract sources
	sources.clear();
	int status= inputImg->FindCompactSource(
		sources,
		significanceMap,bkgData,
		seedThr,mergeThr,minNPix,
		searchNestedSources
	);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Source search failed!");
		#endif
		return -1;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sources.size()<<" sources found...");
	#endif

	return 0;

}//close FindSources()

int SelectSources()
{
	//- Apply source selection
	int nSources= (int)sources.size();
	if(nSources<=0) return 0;
	
	int nSelSources= 0;

	std::vector<Source*> sources_sel;
	for(int i=0;i<nSources;i++)
	{	
		std::string sourceName= sources[i]->GetName();
		int sourceId= sources[i]->Id;
		long int NPix= sources[i]->NPix;
		double X0= sources[i]->X0;
		double Y0= sources[i]->Y0;

		//Is bad source (i.e. line-like blob, etc...)?
		if(!IsGoodSource(sources[i])) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as bad source, skipped!");
			#endif
			sources[i]->SetGoodSourceFlag(false);
			continue;
		}
			
		//Is point-like source?
		if( IsPointLikeSource(sources[i]) ){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as a point-like source ...");
			#endif
			sources[i]->SetType(ePointLike);
		}

		//Tag nested sources
		std::vector<Source*> nestedSources= sources[i]->GetNestedSources();
		for(unsigned int j=0;j<nestedSources.size();j++){
			std::string nestedSourceName= nestedSources[j]->GetName();
			int nestedSourceId= nestedSources[j]->Id;
			long int nestedNPix= nestedSources[j]->NPix;
			double nestedX0= nestedSources[j]->X0;
			double nestedY0= nestedSources[j]->Y0;

			if(!IsGoodSource(nestedSources[j])) {
				#ifdef LOGGING_ENABLED
					INFO_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as bad source, skipped!");
				#endif
				nestedSources[j]->SetGoodSourceFlag(false);
			}
			if( IsPointLikeSource(nestedSources[j]) ){
				#ifdef LOGGING_ENABLED
					INFO_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as a point-like source ...");
				#endif
				nestedSources[j]->SetType(ePointLike);
			}
		}//end loop nested sources
			
		//Add source to the list	
		sources_sel.push_back(sources[i]);
		nSelSources++;
	}//end loop sources

	#ifdef LOGGING_ENABLED
		INFO_LOG("Selected and tagged #"<<nSelSources<<" sources ...");
	#endif

	//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
	sources.clear();
	sources.insert(sources.end(),sources_sel.begin(),sources_sel.end());
	sources_sel.clear();

	return 0;

}//SelectSources()

bool IsGoodSource(Source* aSource)
{	
	if(!aSource) return false;

	//## Check for pixels 	
	if(aSource->NPix<=0 || (aSource->GetPixels()).size()<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No pixels present in this source, cannot perform check!");
		#endif
		return false;
	}

	//## Check for line-like source
	if( (aSource->GetContours()).size()<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contour stored for this source, cannot perform check!");
		#endif
		return true;
	}

	double BoundingBoxMin= ((aSource->GetContours())[0])->BoundingBoxMin;
	if(BoundingBoxMin<sourceMinBoundingBox) {
		#ifdef LOGGING_ENABLED
			INFO_LOG("BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<"<<sourceMinBoundingBox<<")");
		#endif
		return false;
	}

	//## Add other check here ...
	//...
	//...

	return true;

}//close IsGoodSource()

bool IsPointLikeSource(Source* aSource)
{
	if(!aSource) return false;
	if(!aSource->HasParameters()) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("No parameters are available for this source (did you compute them?)...test cannot be performed!");
		#endif
		return true;
	}

	std::string sourceName= aSource->GetName();
	int sourceId= aSource->Id;

	//Loop over contours and check if all of them have circular features
	bool isPointLike= true;
	std::vector<Contour*> contours= aSource->GetContours();

	for(unsigned int i=0;i<contours.size();i++){
		Contour* thisContour= contours[i];

		//Test circularity ratio: 1= circle
		if(useCircRatioCut && thisContour->CircularityRatio<psCircRatioThr) {
			cout<<"SourceFinder::IsCompactSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass CircularityRatio cut (CR="<<thisContour->CircularityRatio<<"<"<<psCircRatioThr<<")"<<endl;
			isPointLike= false;
			break;
		}

		//Test elongation (how symmetrical is the shape): 0=circle,square
		if(useElongCut && thisContour->Elongation>psElongThr) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass Elongation cut (ELONG="<<thisContour->CircularityRatio<<">"<<psElongThr<<")");
			#endif
			isPointLike= false;
			break;	
		}

		//Test ellipse fit
		if(useEllipseAreaRatioCut && (thisContour->EllipseAreaRatio<psEllipseAreaRatioMinThr || thisContour->EllipseAreaRatio>psEllipseAreaRatioMaxThr) ) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass EllipseAreaRatio cut (EAR="<<thisContour->EllipseAreaRatio<<" outside range ["<<psEllipseAreaRatioMinThr<<","<<psEllipseAreaRatioMaxThr<<"])");
			#endif
			isPointLike= false;
			break;	
		}

	}//end contour loop
	

	//Check number of beams contained in source
	double beamArea= aSource->GetBeamFluxIntegral();
	if(useNBeamsCut && beamArea>0){	
		double nBeams= (double)(aSource->NPix)/beamArea;
		if(nBeams>psNBeamsThr){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass nBeams cut (beamArea="<<beamArea<<", NPix="<<aSource->NPix<<", nBeams="<<nBeams<<">"<<psNBeamsThr<<")");
			#endif
			isPointLike= false;
		}
	}

	//Check number of pixels
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") (NPix="<<aSource->NPix<<">"<<psMaxNPix<<")");
	#endif
	if(useMaxNPixCut && aSource->NPix>psMaxNPix){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass nMaxPix cut (NPix="<<aSource->NPix<<">"<<psMaxNPix<<")");
		#endif
		isPointLike= false;
	}

	if(!isPointLike) return false;

	return true;

}//close IsPointLikeSource()

int ComputeStats(Image* img){

	//## Compute stats
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing input image stats...");
	#endif
	bool computeRobustStats= true;
	bool useRange= false;
	bool forceRecomputing= false;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	bool useParallelVersion= false;
	std::vector<float> maskedValues= {0};	
	if(img->ComputeStats(computeRobustStats,forceRecomputing,useRange,minThr,maxThr,useParallelVersion,maskedValues)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Stats computing failed!");
		#endif
		return -1;
	}
	img->PrintStats();	

	return 0;

}//close ComputeStats()

int ComputeBkg(Image* img)
{
	//Get map size
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	
	//Compute sampling box size in pixels 
	//If box size is given as a beam multiple find first the beam size from image
	double boxSizeX_pix= boxSize;
	double boxSizeY_pix= boxSize;
	
	if(boxSizeInBeam)
	{
		//Check if image has metadata with beam info
		if(!img->HasMetaData()){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Requested to use sampling box size as multiple of image beam but input image has no metadata with beam information stored!");
			#endif
			return -1;
		}
	
		//Compute box size in pixels
		int pixelWidthInBeam= img->GetMetaData()->GetBeamWidthInPixel();
		
		#ifdef LOGGING_ENABLED	
			INFO_LOG("Read a beam size of "<<pixelWidthInBeam<<" pixels from input image ...");
		#endif
		if(pixelWidthInBeam>0) {
			boxSizeX_pix= pixelWidthInBeam*boxSizeX;
			boxSizeY_pix= pixelWidthInBeam*boxSizeY;
		}
		else {
			#ifdef LOGGING_ENABLED
				WARN_LOG("Beam information is not available, consider options in pixels ...");
			#endif
		}

	}//close if

	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting bkg boxe sizes to ("<<boxSizeX_pix<<","<<boxSizeY_pix<<") pixels ...");
	#endif
	
	//Compute grid size in pixels
	double gridSizeX_pix= gridSizeX*boxSizeX_pix;
	double gridSizeY_pix= gridSizeY*boxSizeY_pix;
	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting grid size to ("<<gridSizeX_pix<<","<<gridSizeY_pix<<") pixels ...");
	#endif
	
	//Check grid & box size
	if(boxSizeX_pix>=Nx || boxSizeY_pix>=Ny || gridSizeX_pix>=Nx || gridSizeY_pix>=Ny){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Box/grid size are too large compared to image size ("<<Nx<<","<<Ny<<")");
		#endif
		return -1;
	}

	
	/*
	//Beam info
	bool useBeamInfoInBkg;
	if(GET_OPTION_VALUE(useBeamInfoInBkg,useBeamInfoInBkg)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get useBeamInfoInBkg option!");
		#endif
		return -1;
	}
	int nPixelsInBeam= 0;
	if(useBeamInfoInBkg && img->HasMetaData()){
		nPixelsInBeam= img->GetMetaData()->GetBeamWidthInPixel();	
	}
		
	//Compute box size
	double Nx= static_cast<double>(img->GetNx());
	double Ny= static_cast<double>(img->GetNy());
	if(useBeamInfoInBkg && nPixelsInBeam>0){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Setting bkg boxes as ("<<boxSizeX<<","<<boxSizeY<<") x beam (beam="<<nPixelsInBeam<<" pixels) ...");
		#endif
		boxSizeX*= nPixelsInBeam;
		boxSizeY*= nPixelsInBeam;
	}
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Beam information is not available or its usage has been turned off, using image fractions...");
		#endif
		boxSizeX*= Nx;
		boxSizeY*= Ny;
	}
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting bkg boxes to ("<<boxSizeX<<","<<boxSizeY<<") pixels ...");	
	#endif

	gridSizeX*= boxSizeX;
	gridSizeY*= boxSizeY;
	#ifdef LOGGING_ENABLED
		INFO_LOG("Setting grid size to ("<<gridSizeX<<","<<gridSizeY<<") pixels ...");
	#endif

	//## Check grid & box size
	if(boxSizeX>=Nx || boxSizeY>=Ny || gridSizeX>=Nx || gridSizeY>=Ny){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Box/grid size are too large compared to image size ("<<Nx<<","<<Ny<<")");
		#endif
		return -1;
	}
	*/

	//Compute bkg & noise maps
	bkgData= img->ComputeBkg(
		bkgEstimator,
		useLocalBkg,
		boxSizeX_pix,boxSizeY_pix,gridSizeX_pix,gridSizeY_pix,
		use2ndPass,
		skipOutliers,
		seedThr,mergeThr,minNPix
	);

	if(!bkgData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute bkg data!");
		#endif
		return -1;
	}

	//Compute significance
	#ifdef LOGGING_ENABLED
		INFO_LOG("Compute significance map ...");
	#endif
	significanceMap= img->GetSignificanceMap(bkgData,useLocalBkg);
	if(!significanceMap){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute significance map!");
		#endif
		return -1;
	}

	return 0;

}//close ComputeBkg()


int ReadImage()
{
	//## Check given input file and get info
	Caesar::FileInfo info;
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,false)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid input file ("<<inputFileName<<") specified!");
		#endif
		return -1;
	}
	std::string file_extension= info.extension;
	if(file_extension!= ".fits" && file_extension!=".root") {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid file extension ("<<file_extension<<")...nothing to be done!");
		#endif
		return -1;
	}

	//## Read image
	//===== ROOT reading =====
	if(file_extension==".root"){// Read image from ROOT file
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading ROOT input file "<<inputFileName<<"...");
		#endif
		inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Cannot open input file "<<inputFileName<<"!");
			#endif
			return -1;
		}
		inputImg=  (Image*)inputFile->Get(imageName.c_str());
		if(!inputImg){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Cannot get image from input file "<<inputFileName<<"!");
			#endif
			return -1;
		}
	}//close if

	//===== FITS reading =====
	if(file_extension==".fits"){// Read image from FITS file
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading FITS input file "<<inputFileName<<"...");
		#endif
		inputImg= new Caesar::Image; 	
		inputImg->SetNameTitle(imageName.c_str(),imageName.c_str());
		if(inputImg->ReadFITS(inputFileName)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read image from input file "<<inputFileName<<"!");	
			#endif
			return -1;
		}
	}//close else if

	if(!inputImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read image from input file "<<inputFileName<<"!");
		#endif
		return -1;
	}

	return 0;

}//close ReadImage()

/*
int OpenOutputFile()
{
	//Get options
	if(GET_OPTION_VALUE(outputFile,outputFileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get outputFile option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(saveToFile,saveToFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveToFile option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(saveConfig,saveConfig)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveConfig option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(saveResidualMap,saveResidualMap)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveResidualMap option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(saveBkgMap,saveBkgMap)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveBkgMap option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(saveNoiseMap,saveNoiseMap)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveNoiseMap option!");
		#endif
		return -1;
	}		
	if(GET_OPTION_VALUE(saveSignificanceMap,saveSignificanceMap)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveSignificanceMap option!");
		#endif
		return -1;
	}	
	if(GET_OPTION_VALUE(saveInputMap,saveInputMap)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveInputMap option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(saveSaliencyMap,saveSaliencyMap)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveSaliencyMap option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(saveSources,saveSources)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveSources option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(saveToFITSFile,saveToFITSFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saveToFITSFile option!");
		#endif
		return -1;
	}

	if(GET_OPTION_VALUE(residualMapFITSFile,residualMapFITSFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get residualMapFITSFile option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(inputMapFITSFile,inputMapFITSFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("OpenOutputFile(): ERROR: Failed to get inputMapFITSFile option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyMapFITSFile,saliencyMapFITSFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get saliencyMapFITSFile option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(bkgMapFITSFile,bkgMapFITSFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get bkgMapFITSFile option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(noiseMapFITSFile,noiseMapFITSFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get noiseMapFITSFile option!");
		#endif
		return -1;
	}
	if(GET_OPTION_VALUE(significanceMapFITSFile,significanceMapFITSFile)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get significanceMapFITSFile option!");
		#endif
		return -1;
	}

	if(saveToFile){
		outputFile= new TFile(outputFileName.c_str(),"RECREATE");
		if(!outputFile || !outputFile->IsOpen()){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to open output file!");
			#endif
			return -1;
		}
	}
	
	return 0;

}//close OpenOutputFile()
*/

int Clear()
{
	//Clear input image
	if(inputImg) inputImg->Delete();

	//Clear bkg data
	if(bkgData){ 
		delete bkgData;
		bkgData= 0;
	}

	//Delete significance map
	if(significanceMap) significanceMap->Delete();

	//Delete source mask
	if(smaskImg) smaskImg->Delete();

	//Delete sources
	for(size_t i=0;i<sources.size();i++){
		if(sources[i]){
			delete sources[i];
			sources[i]= 0;
		}
	}
	sources.clear();

	//Delete residual image
	if(residualImg) residualImg->Delete();

	return 0;

}//close Clear()

int Save()
{
	/*
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
		if(saveResidualMap && residualImg) {
			residualImg->SetName("img_residual");
			residualImg->Write();	
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
		if(saveResidualMap && residualImg) residualImg->WriteFITS(residualMapFITSFile);
	}
	*/

	//Save bkg map
	if(bkgData->BkgMap){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving bkg map to file "<<outputFileName_bkg);
		#endif	
		(bkgData->BkgMap)->WriteFITS(outputFileName_bkg);
	}

	//Save mask map
	if(smaskImg){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving source mask map to file "<<outputFileName_mask);
		#endif	
		smaskImg->WriteFITS(outputFileName_mask);
	}

	//Save residual map
	if(residualImg){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving residual map to file "<<outputFileName);
		#endif	
		residualImg->WriteFITS(outputFileName);
	}

	return 0;

}//close Save()

std::string GetStringLogLevel(int verb)
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

