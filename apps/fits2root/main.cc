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
#include <FITSReader.h>
#include <SysUtils.h>
#include <Image.h>
#include <Logger.h>

#include <TFITS.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TROOT.h>
#include <TLine.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TMinuit.h>
#include <TMultiDimFit.h>
#include <TCut.h>
#include <TEntryList.h>
#include <TVector3.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TApplication.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TRandom3.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <stdexcept>
#include <chrono>
#include <vector>

#include <getopt.h>

using namespace std;
using namespace Caesar;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input \t Input file name containing image to be read (.fits)"<<endl;
	cout<<"-x, --xmin \t Minimum x pixel id to be read"<<endl;
	cout<<"-X, --xmax \t Maximum x pixel id to be read"<<endl;
	cout<<"-y, --ymin \t Minimum y pixel id to be read"<<endl;
	cout<<"-Y, --ymax \t Maximum y pixel id to be read"<<endl;
	cout<<"-o, --output \t Output file name (ROOT format)"<<endl;
	cout<<"-I, --imgname \t Image name in input ROOT file (if non standard)"<<endl;
	cout<<"-v, --verbosity \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "xmin", optional_argument, 0, 'x' },
	{ "xmax", optional_argument, 0, 'X' },
	{ "ymin", optional_argument, 0, 'y' },
	{ "ymax", optional_argument, 0, 'Y' },
	{ "imgname", optional_argument, 0, 'I' },
	{ "output", optional_argument, 0, 'o' },
	{ "verbosity", required_argument, 0, 'v'},
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string inputFileName= "";
std::string outputFileName= "";
bool readFullImage= true;
bool useDefaultOutput= true;
int minx= -1;
int maxx= -1;
int miny= -1;
int maxy= -1;
std::string imageName= "img";
int verbosity= 4;

//Variables
Caesar::FileInfo info;
Image* inputImg= 0;
//Img* inputImg= 0;
TFile* outputFile;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int OpenOutputFile();
void Clear();
int ReadImage();
void Save();


int main(int argc, char **argv){

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
	INFO_LOG("dt_parse(ms)= "<<dt_parse<<" ["<<dt_parse/dt*100.<<"%]");
	INFO_LOG("dt_outfile(ms)= "<<dt_outfile<<" ["<<dt_outfile/dt*100.<<"%]");
	INFO_LOG("dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]");
	INFO_LOG("dt_save(ms)= "<<dt_save<<" ["<<dt_save/dt*100.<<"%]");
	INFO_LOG("===========================");
	
	
	INFO_LOG("End FITS2ROOT application");


	return 0; 

}//close macro


int ParseOptions(int argc, char *argv[])
{
	//## Check args given
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//## Get command args
	int c = 0;
  int option_index = 0;
	
	
	while((c = getopt_long(argc, argv, "hi:x::X::y::Y::o::I::v:",options_tab, &option_index)) != -1) {
    
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
				minx= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'X':	
			{
				maxx= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'y':	
			{
				miny= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'Y':	
			{
				maxy= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				useDefaultOutput= false;
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

	//## Check coords range in case 
	if(!readFullImage) {
		if(minx>=maxx || miny>=maxy){
			ERROR_LOG("Invalid coord range selected (x["<<minx<<","<<maxx<<"] y["<<miny<<","<<maxy<<"])");
			return -1;
		}
	}

	//## Check given input file and get info
	INFO_LOG("Check input file name "<<inputFileName<<" ...");
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,true,".fits")){
		ERROR_LOG("Invalid input file ("<<inputFileName<<") specified!");
		return -1;
	}

	if(useDefaultOutput){
		std::string basefilename_wext= info.filename_wext;
		outputFileName= basefilename_wext + std::string(".root");
	}

	//## Print options
	INFO_LOG("========= OPTIONS ============");
	INFO_LOG("input file: "<<inputFileName);
	if(!readFullImage) INFO_LOG("x["<<minx<<","<<maxx<<"] y["<<miny<<","<<maxy<<"]");
	INFO_LOG("image name: "<<imageName);
	INFO_LOG("output file: "<<outputFileName);
	INFO_LOG("===============================");
	

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

	//Find output file name
	if(useDefaultOutput){
		std::string basefilename_wext= info.filename_wext;
		outputFileName= basefilename_wext + std::string(".root");
	}

	//Open output file
	outputFile= new TFile(outputFileName.c_str(),"RECREATE");	
	
	return 0;

}//close OpenOutputFile()


int ReadImage(){

	//## Create image
	//inputImg= new Caesar::Img; 
	inputImg= new Caesar::Image; 
	
	//## Create image
	int status= 0;
	if(readFullImage) status= inputImg->ReadFITS(inputFileName);
	else status= inputImg->ReadFITS(inputFileName,minx,maxx,miny,maxy);
	
	if(status<0){
		ERROR_LOG("Failed to read image from FITS file!");
		return -1;
	}

	return 0;

}//close ReadImage()


void Save(){
	
	//Save image to ROOT file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();
		if(inputImg) {
			inputImg->SetName(imageName);
			inputImg->Write();
		}
		outputFile->Close();
	}

}//close Save()

void Clear(){

	//Clear input image
	if(inputImg) inputImg->Delete();

	//Close output file if open			
	if(outputFile && outputFile->IsOpen()) outputFile->Close();
	
}//close Clear()


