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
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <SysUtils.h>
#include <FITSReader.h>


//ROOT headers
#include <TFile.h>
#include <TTree.h>

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
	cout<<"-x, --xmin \t Minimum x pixel id to be read"<<endl;
	cout<<"-X, --xmax \t Maximum x pixel id to be read"<<endl;
	cout<<"-y, --ymin \t Minimum y pixel id to be read"<<endl;
	cout<<"-Y, --ymax \t Maximum y pixel id to be read"<<endl;
	cout<<"-n, --nthreads \t Number of threads to be used for reading (-1=all available threads)"<<endl;
	cout<<"-p, --parallel \t Use parallel std algorithms for median "<<endl;
	cout<<"-o, --output \t Output file name (.root)"<<endl;
	cout<<"-s, --save \t Save perf data to file (.root)"<<endl;
	cout<<"-v, --verbosity \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "xmin", optional_argument, 0, 'x' },
	{ "xmax", optional_argument, 0, 'X' },
	{ "ymin", optional_argument, 0, 'y' },
	{ "ymax", optional_argument, 0, 'Y' },
	{ "nthreads", required_argument, 0, 'n' },
	{ "parallel", no_argument, 0, 'p' },
	{ "output", required_argument, 0, 'o' },
	{ "save", no_argument, 0, 's' },
	{ "verbosity", required_argument, 0, 'v'},
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string inputFileName= "";
int verbosity= 4;//INFO level
int nthreads= -1;
bool useParallelVersion= false;
bool saveToFile= false;
std::string outputFileName= "Output.root";
bool readFullImage= true;
long int minx= -1;
long int maxx= -1;
long int miny= -1;
long int maxy= -1;

//Globar vars
Caesar::Image* inputImg= 0;
std::string imageName= "img";
Caesar::FileInfo info;
TFile* outputFile= 0;
TTree* PerfInfo= 0;
double dt= 0;
double dt_init= 0;
double dt_stats= 0;
double dt_read= 0;
double nx= 0;
double ny= 0;
double NThreads= 0;
double virtMem= 0;
double realMem= 0;
double virtPeakMem= 0;
double realPeakMem= 0;


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int OpenOutputFile();
int ReadImage();
int ComputeStats();
void Clear();
void Save();

//=======================================
//===             MAIN               ====
//=======================================
int main(int argc, char *argv[]){

	auto t0 = chrono::steady_clock::now();
	
	//================================
	//==   INIT
	//================================
	auto t0_init = chrono::steady_clock::now();

	//## Parse command line options
	if(ParseOptions(argc,argv)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to parse command line options!");
		#endif
		Clear();
		return -1;
	}
	
	//## Open output file
	if(saveToFile && OpenOutputFile()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open output file!");
		#endif
		Clear();
		return -1;	
	}

	auto t1_init = chrono::steady_clock::now();
	dt_init= chrono::duration <double, milli> (t1_init-t0_init).count();

	//=======================
	//== Read image
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
	dt_read= chrono::duration <double, milli> (t1_read-t0_read).count();

	//=======================
	//== Compute stats
	//=======================
	auto t0_stats = chrono::steady_clock::now();
	if(ComputeStats()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read image from file!");		
		#endif
		Clear();
		return -1;
	}
	auto t1_stats = chrono::steady_clock::now();
	dt_stats= chrono::duration <double, milli> (t1_stats-t0_stats).count();
	
	//=======================
	//== Read memory usage
	//=======================
	ProcMemInfo memInfo;
	if(SysUtils::GetProcMemoryInfo(memInfo)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read process memory info!");		
		#endif
		Clear();
		return -1;
	}	
	
	virtMem= memInfo.virtMem;
	realMem= memInfo.realMem;
	virtPeakMem= memInfo.virtPeakMem;
	realPeakMem= memInfo.realPeakMem;
	
	//=======================
	//== Fill stats
	//=======================
	auto t1 = chrono::steady_clock::now();
	dt= chrono::duration <double, milli> (t1-t0).count();
	if(PerfInfo) PerfInfo->Fill();

	
	//=======================
	//== Save to file
	//=======================
	auto t0_save = chrono::steady_clock::now();	
	if(saveToFile){
		Save();
	}
	auto t1_save = chrono::steady_clock::now();
	double dt_save= chrono::duration <double, milli> (t1_save-t0_save).count();

	//=======================
	//== Clear data
	//=======================
	Clear();

	#ifdef LOGGING_ENABLED
		INFO_LOG("===========================");
		INFO_LOG("===   PERFORMANCE INFO  ===");
		INFO_LOG("===========================");
		INFO_LOG("dt(ms)= "<<dt);
		INFO_LOG("dt_init(ms)= "<<dt_init<<" ["<<dt_init/dt*100.<<"%]");
		INFO_LOG("dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]");
		INFO_LOG("dt_stats(ms)= "<<dt_stats<<" ["<<dt_stats/dt*100.<<"%]");
		INFO_LOG("dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]");
		INFO_LOG("real mem(kB)= "<<realMem);
		INFO_LOG("real mem peak(kB)= "<<realPeakMem);
		INFO_LOG("virt mem(kB)= "<<virtMem);
		INFO_LOG("virt mem peak(kB)= "<<virtPeakMem);
		INFO_LOG("===========================");
	#else
		cout<<"==========================="<<endl;
		cout<<"===   PERFORMANCE INFO  ==="<<endl;
		cout<<"==========================="<<endl;
		cout<<"dt(ms)= "<<dt<<endl;
		cout<<"dt_init(ms)= "<<dt_init<<" ["<<dt_init/dt*100.<<"%]"<<endl;
		cout<<"dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]"<<endl;
		cout<<"dt_stats(ms)= "<<dt_stats<<" ["<<dt_stats/dt*100.<<"%]"<<endl;
		cout<<"dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]"<<endl;
		cout<<"real mem(kB)= "<<realMem<<endl;
		cout<<"real mem peak(kB)= "<<realPeakMem<<endl;
		cout<<"virt mem(kB)= "<<virtMem<<endl;
		cout<<"virt mem peak(kB)= "<<virtPeakMem<<endl;
		cout<<"==========================="<<endl;
	#endif
	
	#ifdef LOGGING_ENABLED
		INFO_LOG("End image read");
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
	while((c = getopt_long(argc, argv, "hi:n:psv:o:x::X::y::Y:",options_tab, &option_index)) != -1) {
    
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
			case 'x':	
			{
				minx= atol(optarg);	
				readFullImage= false;
				break;	
			}
			case 'X':	
			{
				maxx= atol(optarg);	
				readFullImage= false;
				break;	
			}
			case 'y':	
			{
				miny= atol(optarg);	
				readFullImage= false;
				break;	
			}
			case 'Y':	
			{
				maxy= atol(optarg);	
				readFullImage= false;
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
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
			case 's':	
			{
				saveToFile= true;				
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

	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	#ifdef LOGGING_ENABLED
		LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	#endif

	//## Check coords range in case 
	if(!readFullImage) {
		if(minx>=maxx || miny>=maxy){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Invalid coord range selected (x["<<minx<<","<<maxx<<"] y["<<miny<<","<<maxy<<"])");
			#endif
			return -1;
		}
	}

	//## Set number of threads
	if(nthreads>0) SysUtils::SetOMPThreads(nthreads);
	NThreads= SysUtils::GetOMPMaxThreads();

	//## Print options
	#ifdef LOGGING_ENABLED
		INFO_LOG("========= OPTIONS ============");
		INFO_LOG("input file: "<<inputFileName);
		if(!readFullImage) INFO_LOG("x["<<minx<<","<<maxx<<"] y["<<miny<<","<<maxy<<"]");
		INFO_LOG("image name: "<<imageName);
		INFO_LOG("nthreads: "<<nthreads<<" (NThreads="<<NThreads<<")");
		INFO_LOG("useParallelAlgo? "<<useParallelVersion);
		INFO_LOG("output file: "<<outputFileName);
		INFO_LOG("===============================");
	#endif

	return 0;

}//close ParseOptions()

std::string GetStringLogLevel(int verbosity){

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



int ReadImage(){

	// Check given input file and get info
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

	//--> ROOT reading
	if(file_extension==".root"){// Read image from ROOT file
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading ROOT input file "<<inputFileName<<"...");
		#endif
		TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Cannot open input file "<<inputFileName<<"!");
			#endif
			return -1;
		}
			
		Image* img=  (Caesar::Image*)inputFile->Get(imageName.c_str());
		if(!img){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Cannot get image from input file "<<inputFileName<<"!");
			#endif
			return -1;
		}

		// Get sub-image?
		if(readFullImage){
			inputImg= img;
		}
		else{
			inputImg= img->GetTile(minx,maxx,miny,maxy);
			if(!inputImg){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to read subimage!");
				#endif
				delete img;
				img= 0;
				return -1;	
			}
			delete img;
			img= 0;
		}

	}//close if

	//--> FITS reading
	if(file_extension==".fits"){// Read image from FITS file
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading FITS input file "<<inputFileName<<"...");
		#endif
		inputImg= new Caesar::Image;
		inputImg->SetName(imageName);
		
		int hdu_id= 1;
		int status= 0;
		if(readFullImage) status= inputImg->ReadFITS(inputFileName,hdu_id);
		else inputImg->ReadFITS(inputFileName,hdu_id,minx,maxx,miny,maxy);

		if(status<0){
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

	nx= inputImg->GetNx();
	ny= inputImg->GetNy();

	return 0;

}//close ReadImage()

int ComputeStats(){

	//## Compute stats
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing input image stats...");	
	#endif
	bool computeRobustStats= true;
	bool useRange= false;
	bool forceRecomputing= false;	
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	if(inputImg->ComputeStats(computeRobustStats,forceRecomputing,useRange,minThr,maxThr,useParallelVersion)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Stats computing failed!");
		#endif
		return -1;
	}
	inputImg->PrintStats();	

	return 0;

}//close ComputeStats()


void Clear(){

	//Clear input image
	if(inputImg) inputImg->Delete();
	
	
}//close Clear()


int OpenOutputFile(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"UPDATE");

	//Try to retrieve the TTree from the file
	PerfInfo= (TTree*)outputFile->Get("PerfInfo");

	//Create perf tree
	if(PerfInfo){
		PerfInfo->SetBranchAddress("Nx",&nx);	
		PerfInfo->SetBranchAddress("Ny",&ny);	
		PerfInfo->SetBranchAddress("dt",&dt);	
		PerfInfo->SetBranchAddress("NThreads",&NThreads);	
		PerfInfo->SetBranchAddress("dt_init",&dt_init);
		PerfInfo->SetBranchAddress("dt_read",&dt_read);
		PerfInfo->SetBranchAddress("dt_stats",&dt_stats);	
		PerfInfo->SetBranchAddress("virtMem",&virtMem);
		PerfInfo->SetBranchAddress("realMem",&realMem);
		PerfInfo->SetBranchAddress("virtPeakMem",&virtPeakMem);
		PerfInfo->SetBranchAddress("realPeakMem",&realPeakMem);
	}
	else{//create the TTree as not existing in file
		PerfInfo= new TTree("PerfInfo","PerfInfo");
		PerfInfo->Branch("Nx",&nx,"Nx/D");	
		PerfInfo->Branch("Ny",&ny,"Ny/D");
		PerfInfo->Branch("NThreads",&NThreads,"NThreads/D");	
		PerfInfo->Branch("dt",&dt,"dt/D");	
		PerfInfo->Branch("dt_init",&dt_init,"dt_init/D");
		PerfInfo->Branch("dt_read",&dt_read,"dt_read/D");
		PerfInfo->Branch("dt_stats",&dt_stats,"dt_stats/D");
		PerfInfo->Branch("virtMem",&virtMem,"virtMem/D");
		PerfInfo->Branch("realMem",&realMem,"realMem/D");
		PerfInfo->Branch("virtPeakMem",&virtPeakMem,"virtPeakMem/D");
		PerfInfo->Branch("realPeakMem",&realPeakMem,"realPeakMem/D");
	}

	return 0;

}//close OpenOutputFile()

void Save(){

	//Save to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();
		if(PerfInfo) PerfInfo->Write("PerfInfo");			
		outputFile->Close();
	}

}//close Save()

