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


/**
* @file Image.cc
* @class Image
* @brief Image
*
* Image class
* @author S. Riggi
* @date 20/01/2015
*/

#include <Image.h>
#include <ImgStats.h>
#include <ImgMetaData.h>

//== Img headers
#include <ImgStats.h>

//== IO headers
#include <FITSReader.h>
#include <FITSWriter.h>

//== Utils headers
#include <StatsUtils.h>
#include <GraphicsUtils.h>
#include <CodeUtils.h>
#include <SysUtils.h>
#include <MathUtils.h>
#include <WCSUtils.h>

//== Img processing headers
#include <BkgFinder.h>
#include <BkgData.h>
#include <Source.h>
#include <BlobFinder.h>
#include <ChanVeseSegmenter.h>
#include <Pixel.h>

//== Filter headers
#include <GuidedFilter.h>
#include <WTFilter.h>
#include <KirschFilter.h>
#include <LoGFilter.h>
#include <GradientFilter.h>
#include <MorphFilter.h>
#include <SaliencyFilter.h>
#include <GausFilter.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

//ROOT headers
#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TExec.h>
#include <TF1.h>
#include <TRandom3.h>

//OpenCV headers
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

//VTK headers
#ifdef VTK_ENABLED
	#include <vtkSmartPointer.h>
	#include <vtkImageData.h>
#endif

//C++ headers
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>
#include <chrono>

using namespace std;


ClassImp(Caesar::Image)
ClassImp(Caesar::ImgRange)
ClassImp(Caesar::ImgPeak)
ClassImp(Caesar::ImgBkgPars)


namespace Caesar {

//Default constructor
Image::Image()
	: TNamed()
{
	//Do not allocate memory in default constructor otherwise you will have a memory leak 
	//(see https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html#the-default-constructor)
  Init();
}//close costructor

Image::Image(long int nbinsx,long int nbinsy,std::string name) 
	: TNamed(name.c_str(),name.c_str())
{

	//Check mismatch between pixels size and dimx/dimy
	if(nbinsx<=0 || nbinsy<=0){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) given!";
		#ifdef LOGGING_ENABLED
			ERROR_LOG(ss.str());
		#endif
		throw std::out_of_range(ss.str().c_str());
	}
	
	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,0,0);

}//close constructor

Image::Image(long int nbinsx,long int nbinsy,float xlow,float ylow,std::string name) 
	: TNamed(name.c_str(),name.c_str())
{

	//Check mismatch between pixels size and dimx/dimy
	if(nbinsx<=0 || nbinsy<=0){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) given!";
		#ifdef LOGGING_ENABLED
			ERROR_LOG(ss.str());
		#endif
		throw std::out_of_range(ss.str().c_str());
	}
	
	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,xlow,ylow);

}//close constructor

Image::Image(long int nbinsx,long int nbinsy,std::vector<float>const& pixels,std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Check mismatch between pixels size and dimx/dimy
	long int npixels= (long int)(pixels.size());
	if(nbinsx<=0 || nbinsy<=0 || npixels<=0 || npixels!=(nbinsx*nbinsy)){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) or pixel size given (should be equal to Nx*Ny)!";
		#ifdef LOGGING_ENABLED
			ERROR_LOG(ss.str());
		#endif
		throw std::out_of_range(ss.str().c_str());
	}

	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,0,0);

	//Fill pixels	
	bool useNegativePixInStats= true;
	for(size_t i=0;i<pixels.size();i++){
		double w= pixels[i];
		if(FillPixel(i,w,useNegativePixInStats)<0) continue;
	}

}//close constructor

Image::Image(long int nbinsx,long int nbinsy,std::vector<float>const& pixels,float xlow,float ylow,std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Check mismatch between pixels size and dimx/dimy
	long int npixels= (long int)(pixels.size());
	if(nbinsx<=0 || nbinsy<=0 || npixels<=0 || npixels!=(nbinsx*nbinsy)){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) or pixel size given (should be equal to Nx*Ny)!";
		#ifdef LOGGING_ENABLED
			ERROR_LOG(ss.str());
		#endif
		throw std::out_of_range(ss.str().c_str());
	}

	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,xlow,ylow);

	//Fill pixels	
	bool useNegativePixInStats= true;
	for(size_t i=0;i<pixels.size();i++){
		double w= pixels[i];
		if(FillPixel(i,w,useNegativePixInStats)<0) continue;
	}

}//close constructor

Image::Image(long int nbinsx,long int nbinsy,float w,std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Check mismatch between pixels size and dimx/dimy
	if(nbinsx<=0 || nbinsy<=0){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) or pixel size given (should be equal to Nx*Ny)!";
		#ifdef LOGGING_ENABLED
			ERROR_LOG(ss.str());
		#endif
		throw std::out_of_range(ss.str().c_str());
	}
	long int npixels= nbinsx*nbinsy;

	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,0,0);

	//Fill pixels	
	bool useNegativePixInStats= true;
	for(long int i=0;i<npixels;i++){
		if(FillPixel(i,w,useNegativePixInStats)<0) continue;
	}

}//close constructor
		
Image::Image(long int nbinsx,long int nbinsy,float w,float xlow,float ylow,std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Check mismatch between pixels size and dimx/dimy
	if(nbinsx<=0 || nbinsy<=0){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) or pixel size given (should be equal to Nx*Ny)!";
		#ifdef LOGGING_ENABLED
			ERROR_LOG(ss.str());
		#endif
		throw std::out_of_range(ss.str().c_str());
	}
	long int npixels= nbinsx*nbinsy;

	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,xlow,ylow);

	//Fill pixels	
	bool useNegativePixInStats= true;
	for(long int i=0;i<npixels;i++){
		if(FillPixel(i,w,useNegativePixInStats)<0) continue;
	}

}//close constructor


Image::Image(const TH2& histo,std::string name)
{
	//Get histo size/offset/name
	long int nbinsx= histo.GetNbinsX();
	long int nbinsy= histo.GetNbinsY();
	float xlow= histo.GetXaxis()->GetBinCenter(1);
	float ylow= histo.GetYaxis()->GetBinCenter(1);
	std::string histoName= std::string(histo.GetName());
	
	if(nbinsx<=0 || nbinsy<=0){
		std::stringstream ss;
		ss<<"Invalid histo size (<=0) given (hint: histo shall have dimensions)!";
		#ifdef LOGGING_ENABLED
			ERROR_LOG(ss.str());
		#endif
		throw std::out_of_range(ss.str().c_str());
	}
		
	//Init pars
  Init();

	//Set image name
	if(name=="") m_name= histoName;
	else m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,xlow,ylow);

	//Fill from TH2
	bool useNegativePixInStats= true;
	FillFromTH2(histo,useNegativePixInStats);	

}//close constructor from TH2


Image::Image(const Image &img)
{
	// Copy constructor
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy constuctor called...");
	#endif
  ((Image&)img).Copy(*this);
}//close copy constructor



void Image::Copy(TObject& obj) const
{
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying parent TNamed...");
	#endif
	TNamed::Copy((Image&)obj);

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying main vars...");
	#endif
	((Image&)obj).m_name= m_name;
	((Image&)obj).m_Nx= m_Nx;
	((Image&)obj).m_Ny= m_Ny;
	((Image&)obj).m_Xmin= m_Xmin;
	((Image&)obj).m_Ymin= m_Ymin;
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying pixel collection...");
	#endif
	((Image&)obj).m_pixels= m_pixels;
		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying stats vars...");
	#endif
	((Image&)obj).m_HasMetaData= m_HasMetaData;
	((Image&)obj).m_HasStats= m_HasStats;
	((Image&)obj).m_StatMoments= m_StatMoments;
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying meta data...");
	#endif
	if(m_MetaData){
		((Image&)obj).m_MetaData= new ImgMetaData;
		*((Image&)obj).m_MetaData = *m_MetaData;
	}
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying stats data...");
	#endif
	((Image&)obj).m_StatMoments= m_StatMoments;
	if(m_Stats){
		((Image&)obj).m_Stats= new ImgStats;
		*((Image&)obj).m_Stats = *m_Stats;
	}
	
}//close Copy()


Image::~Image()
{
	if(m_Stats){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Deleting stats...");	
		#endif
		delete m_Stats;
		m_Stats= 0;
	}

	if(m_MetaData){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Deleting meta-data...");
		#endif
		delete m_MetaData;
		m_MetaData= 0;	
	}	
 	
}//close destructor


Image& Image::operator=(const Image &img) { 
	// Operator =
  if (this != &img) ((Image&)img).Copy(*this);
  return *this;
}

//================================================================
//===    INITIALIZATION METHODS
//================================================================
void Image::Init(){
		
	//Set sizes
	m_Nx= m_Ny= 0;
	m_pixels.clear();
	m_Xmin= 0;
	m_Ymin= 0;

	//Meta-Data
	m_HasMetaData= false;
	m_MetaData= 0;
 
	//Stats
	m_HasStats= false;
	m_Stats= 0;
	m_StatMoments.Reset();
	ResetImgStats(true);//Reset stats & moments

}//close Init()


int Image::SetSize(long int Nx,long int Ny,float xlow,float ylow){
			
	//Check dimensions
	if(Nx<=0 || Ny<=0) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid (empty or negative) sizes given, nothing will be done!");
		#endif
		return -1;
	}

	//Check range
	//bool hasValidOffset= (xlow>=0 && ylow>=0); 
	//if(!hasValidOffset){
	//	ERROR_LOG("Invalid xlow/ylow given (should be >=0)!");
	//	return -1;
	//}

	//Reset & delete stats, reset moments
	ResetImgStats(true,true);
	
	//Allocate new space
	long int N= Nx*Ny;
	try {
		m_pixels.resize(N);
	}
	catch(std::exception& e){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("C++ exception occurred while allocating memory (N="<<N<<") for pixel vector!");
		#endif
		return -1;
	}

	//Init sizes
	m_Nx= Nx;
	m_Ny= Ny;

	//Init min coordinates
	m_Xmin= xlow;
	m_Ymin= ylow;

	return 0;

}//close SetSize()

//================================================================
//===    READ/WRITE METHODS
//================================================================
int Image::ReadFITS(std::string filename,int hdu_id,int ix_min,int ix_max,int iy_min,int iy_max){

	//Start timer		
	auto start = chrono::steady_clock::now();

	//Read image or tile
	bool check_file= true;
	Caesar::FITSFileInfo fits_info;
	int status= FITSReader::Read(*this,fits_info,filename,hdu_id,ix_min,ix_max,iy_min,iy_max,check_file);
	
	//Stop timer and print
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Read FITS image "<<filename<<" in "<<dt<<" ms");	
	#endif

	//Check error status
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to fill image from FITS file!");
		#endif
		return -1;
	}

	return 0;

}//close ReadFITS()


int Image::WriteFITS(std::string outfilename){
		
	//## Write fits to file
	if(FITSWriter::WriteFITS(this,outfilename)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to write image to FITS file!");
		#endif
		return -1;
	}
	
	return 0;

}//close WriteFITS()


int Image::ReadFile(std::string filename,bool invert){

	//## Detect file extension
	std::string extension= filename.substr(filename.find_last_of(".") + 1);
	if(extension!= "png" && extension!="jpg" && extension!="bmp" && extension!="gif" ) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Unknown file extension detected: ext="<<extension<<" (valid ones are png/jpg/bmp/gif)!");
		#endif
		return -1;
	}
	
	//## Load image from file and set a matrix
	//cv::Mat mat = cv::imread(filename.c_str(), CV_LOAD_IMAGE_COLOR);//deprecated C API
	cv::Mat mat = cv::imread(filename.c_str(), cv::IMREAD_COLOR);

	//## Convert to gray scale
	cv::Mat mat_gray;
  //cvtColor( mat, mat_gray, CV_RGB2GRAY );//deprecated C API
	cvtColor( mat, mat_gray, cv::COLOR_RGB2GRAY);

	//## Fill an image
	int Nx= mat.cols;
	int Ny= mat.rows;
	
	this->SetSize(Nx,Ny);
	
	if(invert){
		for(int j=0;j<mat_gray.rows;j++){
			int rowId= Ny-1-j;
			for(int i=0;i<mat_gray.cols;i++){
				int colId= Nx-1-i;
				unsigned int matrixElement= mat_gray.at<uchar>(j,i);	
				long int ix= colId ;
				long int iy= rowId ;
				this->FillPixel(ix,iy,matrixElement);
			}
		}
	}//close if
	else{
		for(int j=0;j<mat_gray.rows;j++){
			int rowId= Ny-1-j;
			for(int i=0;i<mat_gray.cols;i++){
				int colId= i;
				unsigned int matrixElement= mat_gray.at<uchar>(j,i);	
				long int ix= colId ;
				long int iy= rowId ;
				this->FillPixel(ix,iy,matrixElement);
			}
		}
	}

	return 0;

}//close ReadFile()

Image* Image::GetTile(long int ix_min,long int ix_max,long int iy_min,long int iy_max,std::string imgname)
{
	//Extract pixel rectangular selection
	//NB: Include all pixels otherwise the image size would be not the one expected from ix_min/ix_max/iy_min/iy_max
	std::vector<float> tile_pixels;
	bool useRange= false;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	std::vector<float> maskedValues= {};//include also 0!
	bool requireFiniteValues= false;//include also inf/nan

	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,useRange,minThr,maxThr,maskedValues,requireFiniteValues)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to extract pixel rectangular data, returning nullptr!");
		#endif
		return nullptr;
	}

	long int TileSizeX= ix_max - ix_min + 1;
	long int TileSizeY= iy_max - iy_min + 1;
	float xlow= ix_min;
	float ylow= iy_min;
	std::string name= imgname;
	if(name=="") name= "tileimg";

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("ix min/max="<<ix_min<<"/"<<ix_max<<", iy min/max="<<iy_min<<"/"<<iy_max<<" TileSizeX="<<TileSizeX<<", TileSizeY="<<TileSizeY<<", pixels size="<<tile_pixels.size());
	#endif

	Image* tile= new Image(TileSizeX,TileSizeY,tile_pixels,xlow,ylow,name);
	tile->CopyMetaData(this->m_MetaData);
	
	return tile;

}//close GetTile()

//================================================================
//===    FILLING METHODS
//================================================================
int Image::CheckFillPixel(long int gbin,double w){

	//Check bin & value
	if(gbin<0 || gbin>=(long int)(m_pixels.size())) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid pixel bin ("<<gbin<<") given to be filled (hint: check if image size is not initialized or given bin exceed size)!");
		#endif
		return -1;
	}
	if(TMath::IsNaN(w) || fabs(w)==TMath::Infinity()) {
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Given value (w="<<w<<") for bin "<<gbin<<" is NaN or inf, skipping!");
		#endif
		return 0;
	}

	//Compute binx & biny
	long int binx = GetBinX(gbin);
 	long int biny = GetBinY(gbin);

	//Check if pixel value has been already filled
	double w_old= m_pixels[gbin];
	if( w_old!=0 ){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Pixel "<<gbin<<" ("<<binx<<","<<biny<<") has been already filled (w="<<w_old<<"), skipping...");
		#endif
		return -1;
	}
	
	return 0;

}//close CheckFillPixel()


int Image::FillPixel(long int gbin,double w,bool useNegativePixInStats){

	//Check bin & value
	if(CheckFillPixel(gbin,w)<0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid bin given (gbin="<<gbin<<") or or pixel already filled!");
		#endif
		return -1;
	}

	//Set pixel value
	m_pixels[gbin]= w;

	//Update moments
	if(w>=0 || (w<0 && useNegativePixInStats)) {
		Caesar::StatsUtils::UpdateMoments(m_StatMoments,w);
	}

	return 0;

}//close FillPixel()


int Image::FillPixel(long int ix,long int iy,double w,bool useNegativePixInStats){

	//Compute global bin
	long int gbin= GetBin(ix,iy);
	
	//Fill pixel
	return FillPixel(gbin,w,useNegativePixInStats);

}//close FillPixel()


int Image::Fill(double x,double y,double w,bool useNegativePixInStats)
{
	//Find global bin
	long int gbin= FindBin(x,y);
	if(gbin<0 || gbin>=(long int)(m_pixels.size())) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Underflow/overflow bin ("<<gbin<<") corresponding to (x,y)=("<<x<<","<<y<<") given to be filled (hint: check if image size is not initialized)!");
		#endif
		return -1;
	}
	
	//Fill pixel
	return FillPixel(gbin,w,useNegativePixInStats);

}//close FillPixel()


#ifdef OPENMP_ENABLED
int Image::FillPixelMT(Caesar::StatMoments<double>& moments,long int gbin,double w,bool useNegativePixInStats)
{
	//Check bin & value
	if(CheckFillPixel(gbin,w)<0) {
		#ifdef LOGGING_ENABLED			
			WARN_LOG("Invalid bin given (gbin="<<gbin<<") or nan/inf value given (w="<<w<<") or pixel already filled!");
		#endif
		return -1;
	}

	//Set pixel value
	m_pixels[gbin]= w;

	//Update moments
	if(w>=0 || (w<0 && useNegativePixInStats)) {
		Caesar::StatsUtils::UpdateMoments(moments,w);
	}

	return 0;

}//close FillPixelMT()
#endif

#ifdef OPENMP_ENABLED
int Image::FillPixelMT(Caesar::StatMoments<double>& moments,long int ix,long int iy,double w,bool useNegativePixInStats)
{
	//Compute global bin
	long int gbin= GetBin(ix,iy);
	
	//Fill pixel
	return FillPixelMT(moments,gbin,w,useNegativePixInStats);

}//close FillPixelMT()
#endif

#ifdef OPENMP_ENABLED
int Image::FillMT(Caesar::StatMoments<double>& moments,double x,double y,double w,bool useNegativePixInStats)
{
	//Find global bin
	long int gbin= FindBin(x,y);
	if(gbin<0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Underflow/overflow bin requested to be filled!");	
		#endif
	}

	//Fill pixel
	return FillPixelMT(moments,gbin,w,useNegativePixInStats);

}//close FillPixelMT()
#endif


int Image::FillFromMat(cv::Mat& mat,bool useNegativePixInStats)
{
	//Get image size
	long int nRows = mat.rows;
  long int nCols = mat.cols;
	long int Ny= this->GetNy();
	long int Nx= this->GetNx();
	if(nRows<=0 || nCols<=0){
		#ifdef LOGGING_ENABLED		
			ERROR_LOG("Matrix has zero rows and/or cols!");
		#endif
		return -1;
	}
	if(nCols!=Nx || nRows!=Ny){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid nrows/ncols (hint: should be equal to image size "<<Nx<<" x "<<Ny<<")!");
		#endif
		return -1;
	}


	//Reset this image
	Reset();
		
	#ifdef OPENMP_ENABLED
		
		Caesar::StatMoments<double> moments_t;		
		std::vector<Caesar::StatMoments<double>> parallel_moments;

		#pragma omp parallel private(moments_t)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			#pragma omp for
			for(long int i=0;i<nRows;++i) {
				long int rowId= i;
				long int iy= Ny-1-rowId;
				double* p = mat.ptr<double>(i);
    		for (long int j=0;j<nCols;++j){
					long int colId= j;
					long int ix= colId;
					double w= p[j];
					this->FillPixelMT(moments_t,ix,iy,w,useNegativePixInStats);
    		}//end loop cols
  		}//end loop rows
		
			//Fill parallel moments per thread
			parallel_moments[thread_id]= moments_t;
			
		}//close parallel section

		//Update moments from parallel estimates
		Caesar::StatMoments<double> moments;
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(moments,parallel_moments)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
			#endif
			return -1;
		}	
		this->SetMoments(moments);

	#else	
		for(long int i=0;i<nRows;++i) {
			long int rowId= i;
			long int iy= Ny-1-rowId;
			double* p = mat.ptr<double>(i);
    	for (long int j=0;j<nCols;++j){
				long int colId= j;
				long int ix= colId;
				double w= p[j];
				this->FillPixel(ix,iy,w,useNegativePixInStats);
    	}//end loop cols
  	}//end loop rows
	#endif

	return 0;

}//close FillFromMat()


int Image::FillFromTMatrix(TMatrixD& mat,bool useNegativePixInStats){

	//Get image size
	long int nRows = mat.GetNrows();
  long int nCols = mat.GetNcols();
	long int Ny= this->GetNy();
	long int Nx= this->GetNx();
	if(nRows<=0 || nCols<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Matrix has zero rows and/or cols!");
		#endif
		return -1;
	}
	if(nCols!=Nx || nRows!=Ny){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid nrows/ncols (hint: should be equal to image size "<<Nx<<" x "<<Ny<<")!");
		#endif
		return -1;
	}

	//Reset this image
	Reset();

	//Fill image from TMatrixD
	#ifdef OPENMP_ENABLED
		Caesar::StatMoments<double> moments_t;		
		std::vector<Caesar::StatMoments<double>> parallel_moments;

		#pragma omp parallel private(moments_t)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			#pragma omp for
			for(long int i=0;i<nRows;++i) {
				long int rowId= i;
				long int iy= Ny-1-rowId;
    		for (long int j=0;j<nCols;++j){
					long int colId= j;
					long int ix= colId;
					double w= mat(i,j);
					this->FillPixelMT(moments_t,ix,iy,w,useNegativePixInStats);
    		}//end loop cols
  		}//end loop rows
		
			//Fill parallel moments per thread
			parallel_moments[thread_id]= moments_t;
			
		}//close parallel section

		//Update moments from parallel estimates
		Caesar::StatMoments<double> moments;
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(moments,parallel_moments)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
			#endif
			return -1;
		}	
		this->SetMoments(moments);
		
	#else
		for(long int i=0;i<nRows;++i) {
			long int rowId= i;
			long int iy= Ny-1-rowId;
    	for (long int j=0;j<nCols;++j){
				long int colId= j;
				long int ix= colId;
				double w= mat(i,j);
				this->FillPixel(ix,iy,w,useNegativePixInStats);
    	}//end loop cols
  	}//end loop rows
	#endif

	return 0;

}//close FillFromTMatrix()


int Image::FillFromTH2(const TH2& histo,bool useNegativePixInStats){

	//Get image size
	long int nBinsX = histo.GetNbinsX();
  long int nBinsY = histo.GetNbinsY();
	long int Ny= this->GetNy();
	long int Nx= this->GetNx();
	if(nBinsX<=0 || nBinsY<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Histo has zero bins in one/both dimensions!");
		#endif
		return -1;
	}
	if(nBinsX!=Nx || nBinsY!=Ny){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid bins given (hint: should be equal to image size "<<Nx<<" x "<<Ny<<")!");
		#endif
		return -1;
	}

	//Reset this image
	Reset();

	//Fill image from TMatrixD
	#ifdef OPENMP_ENABLED
		Caesar::StatMoments<double> moments_t;		
		std::vector<Caesar::StatMoments<double>> parallel_moments;

		#pragma omp parallel private(moments_t)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			#pragma omp for
			for(long int i=0;i<nBinsX;i++) {
				long int ix= i;
    		for (long int j=0;j<nBinsY;j++){
					long int iy= j;
					double w= histo.GetBinContent(i+1,j+1);
					this->FillPixelMT(moments_t,ix,iy,w,useNegativePixInStats);
    		}//end loop cols
  		}//end loop rows
		
			//Fill parallel moments per thread
			parallel_moments[thread_id]= moments_t;
			
		}//close parallel section

		//Update moments from parallel estimates
		Caesar::StatMoments<double> moments;
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(moments,parallel_moments)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
			#endif
			return -1;
		}	
		this->SetMoments(moments);
		
	#else
		for(long int i=0;i<nBinsX;i++) {
			long int ix= i;
    	for (long int j=0;j<nBinsY;j++){
				long int iy= j;
				double w= histo.GetBinContent(i+1,j+1);
				this->FillPixel(ix,iy,w,useNegativePixInStats);
    	}//end loop cols
  	}//end loop rows
	#endif

	return 0;

}//close FillFromTH2()


#ifdef VTK_ENABLED
int Image::FillFromVtkImage(vtkSmartPointer<vtkImageData> imageData,bool useNegativePixInStats)
{
	//Check image data
	if(!imageData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input VTK image data given!");
		#endif
		return -1;
	}

	//Check dims
	int* dims = imageData->GetDimensions();
	if(dims[2]!=1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid z dimension (must be =1) in VTK image data!");
		#endif
		return -1;
	}

	long int Nx= dims[1];
	long int Ny= dims[0];

	//Set image size
	this->SetSize(Nx,Ny);

	//Reset image data
	this->Reset();	
	
	//Loop over pixel data and fill
	#ifdef OPENMP_ENABLED
		Caesar::StatMoments<double> moments_t;		
		std::vector<Caesar::StatMoments<double>> parallel_moments;

		#pragma omp parallel private(moments_t)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			#pragma omp for
			for(int row = 0; row < dims[0]; ++row) {
				long int iy= Ny-1-row;
    		for(int col = 0; col < dims[1]; ++col) {
					long int ix= col;
					double* pixel = static_cast<double*>(imageData->GetScalarPointer(row, col, 0));
					double w= pixel[0];
					this->FillPixelMT(moments_t,ix,iy,w,useNegativePixInStats);
				}//end loop columns
			}//end loop rows

			//Fill parallel moments per thread
			parallel_moments[thread_id]= moments_t;
			
		}//close parallel section

		//Update moments from parallel estimates
		Caesar::StatMoments<double> moments;
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(moments,parallel_moments)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
			#endif
			return -1;
		}	
		this->SetMoments(moments);

	#else
		for(int row = 0; row < dims[0]; ++row) {
			long int iy= Ny-1-row;
    	for(int col = 0; col < dims[1]; ++col) {
				long int ix= col;
				double* pixel = static_cast<double*>(imageData->GetScalarPointer(row, col, 0));
				double w= pixel[0];
				this->FillPixel(ix,iy,w,useNegativePixInStats);
			}//end loop columns
		}//end loop rows
	#endif

	return 0;

}//close FillFromVtkImage()
#endif

//================================================================
//===    STATS METHODS
//================================================================
void Image::ResetImgStats(bool resetMoments,bool clearStats){

	//Reset stats
	if(m_Stats) m_Stats->Reset();

	//Reset moments
	if(resetMoments){
		m_StatMoments.Reset();
	}

	//Delete current stats data
	if(clearStats){
		ClearImgStats();
	}

}//close ResetImgStats()

int Image::ComputeMoments(bool useRange,double minThr,double maxThr,std::vector<float> maskedValues)
{
	#ifdef LOGGING_ENABLED	
		DEBUG_LOG("m_StatMoments min/max (before): "<<m_StatMoments.minVal<<"/"<<m_StatMoments.maxVal);
	#endif
	
	//## Recompute stat moments
	//## NB: If OMP is enabled this is done in parallel and moments are aggregated to return the correct cumulative estimate 
	//bool maskNanInfValues= true;
	int status= Caesar::StatsUtils::ComputeStatsMoments(m_StatMoments,m_pixels,useRange,minThr,maxThr,maskedValues);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute stat moments!");
		#endif
	}

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("m_StatMoments min/max (after): "<<m_StatMoments.minVal<<"/"<<m_StatMoments.maxVal);
	#endif

	return status;

}//close ComputeMoments()



void Image::ComputeStatsParams(bool computeRobustStats,bool useRange,double minThr,double maxThr,bool useParallelVersion,std::vector<float> maskedValues)
{
	//## Reset previous stats params (Reset only Stats not moments!)
	ResetImgStats(false);

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Npix="<<m_StatMoments.N<<", M1="<<m_StatMoments.M1<<", M2="<<m_StatMoments.M2<<", M3="<<m_StatMoments.M3<<", M4="<<m_StatMoments.M4<<", pixel sizes="<<m_pixels.size());
	#endif

	//##########################################
	//##  COMPUTE STANDARD STATS PARAMETERS
	//##########################################
	auto start_stats = chrono::steady_clock::now();
	
	m_Stats->n= m_StatMoments.N;
	m_Stats->min= m_StatMoments.minVal;
	m_Stats->max= m_StatMoments.maxVal;
	m_Stats->mean= m_StatMoments.M1;
	m_Stats->rms= 0;
	if(m_StatMoments.N>2) m_Stats->rms= sqrt(m_StatMoments.M2/(m_StatMoments.N-1));
	m_Stats->skewness= 0;
	m_Stats->kurtosis= 0;
	if(m_StatMoments.M2!=0) {
  	m_Stats->skewness= sqrt(m_StatMoments.N)*m_StatMoments.M3/pow(m_StatMoments.M2,1.5);//need to adjust for finite population?
		m_Stats->kurtosis= m_StatMoments.N*m_StatMoments.M4/(m_StatMoments.M2*m_StatMoments.M2)-3;
	}

	//## Compute stats param errors
	m_Stats->meanErr= 0;
	if(m_StatMoments.N>0) m_Stats->meanErr= (m_Stats->rms)/sqrt(m_StatMoments.N);
	double varianceErr= 0;
	m_Stats->rmsErr= 0;
	if(m_StatMoments.N>1) {
		varianceErr= (m_StatMoments.M4-(m_StatMoments.N-3)/(m_StatMoments.N-1)*pow(m_Stats->rms,4))/m_StatMoments.N;
		m_Stats->rmsErr= varianceErr/(2*m_Stats->rms);
	}	 
	m_Stats->skewnessErr= 0;
	m_Stats->kurtosisErr= 0;
	if(m_StatMoments.N>2) m_Stats->skewnessErr= sqrt(6.*m_StatMoments.N*(m_StatMoments.N-1)/((m_StatMoments.N-2)*(m_StatMoments.N+1)*(m_StatMoments.N+3)));//approximate for normal distribution
	
	auto end_stats = chrono::steady_clock::now();
	double dt_stats= chrono::duration <double, milli> (end_stats-start_stats).count();

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image standard stats computed in "<<dt_stats<<" ms");
  #endif

	//## End if no robust stats are to be computed
	if(!computeRobustStats) return;

	//##########################################
	//##  COMPUTE ROBUST STATS PARAMETERS
	//##########################################
	auto start_robuststats = chrono::steady_clock::now();

	//## Remove negative values? 
	//## NB: Copy vector otherwise it is modified by sorting operation inside median and other robust estimators
	std::vector<float> pixels;
	for(size_t i=0;i<m_pixels.size();i++){
		float w= m_pixels[i];
		//bool isNormalValue= std::isnormal(w);//exclude 0 and inf/nan
		bool isFiniteValue= std::isfinite(w);
		bool isMaskedValue= false;
		for(size_t l=0;l<maskedValues.size();l++){
			if(w==maskedValues[l]){
				isMaskedValue= true;
				break;
			}
		}
		
		//if( w==0 || (useRange && (w<=minThr || w>=maxThr)) ) continue;
		if( !isFiniteValue || isMaskedValue || (useRange && (w<=minThr || w>=maxThr)) ) continue;
		pixels.push_back(w);
	}
	
	//## Compute robust stats (median, MAD, ...)	
	//Sort and compute median for all image	
	auto start_median = chrono::steady_clock::now();
	
	float median= Caesar::StatsUtils::GetMedianFast<float>(pixels,useParallelVersion);
	m_Stats->median= median;

	//Compute MAD = median(|x_i-median|)
	float medianMAD= Caesar::StatsUtils::GetMADFast(pixels,median);	
	double medianRMS= medianMAD*1.4826;//0.6744888;
	m_Stats->medianRMS= medianRMS;
	auto end_median = chrono::steady_clock::now();
	double dt_median= chrono::duration <double, milli> (end_median-start_median).count();

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image median pars computed in "<<dt_median<<" ms");
  #endif

	//## Compute biweight robust estimators	
	/*
	auto start_biweights = chrono::steady_clock::now();
	double C= 6.;
	double tol= 0.0001;
	double nmaxIter= 10;
	std::pair<float,float> biweightEstimators= Caesar::StatsUtils::GetBiWeightEstimators<float>(pixels,median,medianRMS,C,tol,nmaxIter);
	m_Stats->bwLocation= biweightEstimators.first;
	m_Stats->bwScale= biweightEstimators.second;
	auto end_biweights = chrono::steady_clock::now();
	double dt_biweights= chrono::duration <double, milli> (end_biweights-start_biweights).count();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image biweight pars computed in "<<dt_biweights<<" ms");
	#endif
  */

	//## Compute clipped estimators
	auto start_clipped = chrono::steady_clock::now();
	double clipSigma= 3;
	int clipMaxIter= 100;
	double clipTolerance= 0.1;
	float mean= m_Stats->mean;
	float rms= m_Stats->rms;
	ClippedStats<float> clipped_stats;
	Caesar::StatsUtils::GetClippedEstimators(clipped_stats,pixels,median,mean,rms,clipSigma,clipMaxIter,clipTolerance,useParallelVersion);
	m_Stats->clippedMedian= clipped_stats.median;
	m_Stats->clippedRMS= clipped_stats.stddev;
	auto end_clipped = chrono::steady_clock::now();
	double dt_clipped= chrono::duration <double, milli> (end_clipped-start_clipped).count();

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image clipped stats pars computed in "<<dt_clipped<<" ms");
	#endif

	auto end_robuststats = chrono::steady_clock::now();
	double dt_robuststats= chrono::duration <double, milli> (end_robuststats-start_robuststats).count();

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image robust stats computed in "<<dt_robuststats<<" ms");	
	#endif

}//close ComputeStatsParams()



int Image::ComputeStats(bool computeRobustStats,bool forceRecomputing,bool useRange,double minThr,double maxThr,bool useParallelVersion,std::vector<float> maskedValues)
{
	//## Start timer
	auto start = chrono::steady_clock::now();

	//## Check if image has already stats computed
	if(!HasStats()){
		m_Stats= new Caesar::ImgStats;
	}
	else{	
		#ifdef LOGGING_ENABLED	
			DEBUG_LOG("Image has stats already computed...");
		#endif
	}

	//## If recomputing is not requested (i.e. some pixels has been reset by the user, just set the stats params!
	if(!forceRecomputing){
		ComputeStatsParams(computeRobustStats,useRange,minThr,maxThr,useParallelVersion,maskedValues);
		m_HasStats= true;
		
		auto stop = chrono::steady_clock::now();
		double elapsed_time= chrono::duration <double, milli> (stop-start).count();
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Image stats computed in "<<elapsed_time<<" ms");
		#endif
		return 0;
	}

	
	//## Recompute the moments and stats params
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Recomputing image stats...");
	#endif

	//--> Reset stats
	ResetImgStats(true);

	//--> Recompute moments
	ComputeMoments(useRange,minThr,maxThr,maskedValues);

	//--> Recompute stats params
	ComputeStatsParams(computeRobustStats,useRange,minThr,maxThr,useParallelVersion,maskedValues);

	m_HasStats= true;

	//Stop timer
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image stats recomputed in "<<dt<<" ms");
	#endif

	return 0;

}//close ComputeStats()


int Image::GetTilePixels(std::vector<float>& pixels,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool useRange,double minThr,double maxThr,std::vector<float> maskedValues,bool requireFinitePixValues)
{
	//Check range given
	if(ix_min<0 || ix_max<0 || ix_max>=m_Nx || ix_min>=ix_max){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid x range given!");
		#endif
		return -1;
	}
	if(iy_min<0 || iy_max<0 || iy_max>=m_Ny || iy_min>=iy_max){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid y range given!");
		#endif
		return -1;
	}

	//Check thr
	if(minThr>=maxThr){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid min/max threshold given (hint: max thr must be larger then min thr)!");
		#endif
		return -1;
	}
	
	//Extract pixel sub vector
	pixels.clear();

	if(useRange)
	{
		for(long int j=iy_min;j<=iy_max;j++)
		{
			//Get row start/end iterators
			long int gBin_min= GetBin(ix_min,j);		
			long int gBin_max= GetBin(ix_max,j);

			std::copy_if (
				m_pixels.begin() + gBin_min, 
				m_pixels.begin() + gBin_max + 1, 
				std::back_inserter(pixels), 
				[&maskedValues,&minThr,&maxThr,&requireFinitePixValues](float w){
					bool goodValue= true;
					if(requireFinitePixValues) goodValue= std::isfinite(w);
					bool inRange= (w>minThr && w<maxThr);
					bool notMasked= true;
					for(size_t l=0;l<maskedValues.size();l++){
						if(w==maskedValues[l]){
							notMasked= false;
							break;
						}
					}
					bool selected= goodValue && inRange && notMasked;
					
					return selected;
				}
			);
		}//end loop rows
	}//close if
	else
	{
		for(long int j=iy_min;j<=iy_max;j++)
		{
			//Get row start/end iterators
			long int gBin_min= GetBin(ix_min,j);		
			long int gBin_max= GetBin(ix_max,j);		

			//Extract row and append to vector
			//pixels.insert(pixels.end(), m_pixels.begin() + gBin_min, m_pixels.begin() + gBin_max + 1);

			std::copy_if (
				m_pixels.begin() + gBin_min, 
				m_pixels.begin() + gBin_max + 1, 
				std::back_inserter(pixels), 
				[&maskedValues,&minThr,&maxThr,&requireFinitePixValues](float w){
					bool goodValue= true;
					if(requireFinitePixValues) goodValue= std::isfinite(w);
					bool notMasked= true;
					for(size_t l=0;l<maskedValues.size();l++){
						if(w==maskedValues[l]){
							notMasked= false;
							break;
						}
					}
					bool selected= goodValue && notMasked;
					return selected;
				}
			);

		}//end loop row
	}//close else
	
	return 0;

}//close GetTilePixels()


int Image::GetTileMeanStats(float& mean,float& stddev,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool useRange,double minThr,double maxThr,std::vector<float> maskedValues)
{
	//Init
	mean= 0;
	stddev= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	bool requireFiniteValues= true;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,useRange,minThr,maxThr,maskedValues,requireFiniteValues)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to extract pixel rectangular data!");
		#endif
		return -1;
	}
	
	//Compute mean/rms
	Caesar::StatsUtils::ComputeMeanAndRMS(mean,stddev,tile_pixels);
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileMeanStats()
		

int Image::GetTileMedianStats(float& median,float& mad_rms,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool useRange,double minThr,double maxThr,std::vector<float> maskedValues)
{
	//Init
	median= 0;
	mad_rms= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	bool requireFiniteValues= true;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,useRange,minThr,maxThr,maskedValues,requireFiniteValues)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to extract pixel rectangular data!");
		#endif
		return -1;
	}
	
	//Compute median/mad
	bool useParallelVersion= false;
	median= Caesar::StatsUtils::GetMedianFast(tile_pixels,useParallelVersion);
	mad_rms= Caesar::StatsUtils::GetMADFast(tile_pixels,median,useParallelVersion)*1.4826;
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileMedianStats()


int Image::GetTileClippedStats(ClippedStats<float>& clipped_stats,long int& npix,double clipSigma,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool useRange,double minThr,double maxThr,std::vector<float> maskedValues)
{
	//Init
	npix= 0;
	
	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	bool requireFiniteValues= true;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,useRange,minThr,maxThr,maskedValues,requireFiniteValues)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to extract pixel rectangular data!");
		#endif
		return -1;
	}
	
	//Compute mean, rms, median
	float mean= 0;
	float rms= 0;
	Caesar::StatsUtils::ComputeMeanAndRMS(mean,rms,tile_pixels);
	float median= Caesar::StatsUtils::GetMedianFast<float>(tile_pixels);
	
	//Compute clipped stats
	int clipMaxIter= 100;
	double clipTolerance= 0.1;
	bool useParallelVersion= false;
	Caesar::StatsUtils::GetClippedEstimators(clipped_stats,tile_pixels,median,mean,rms,clipSigma,clipMaxIter,clipTolerance,useParallelVersion);
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileClippedStats()
		

int Image::GetTileBiWeightStats(float& bwLocation,float& bwScale,long int& npix,double C,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool useRange,double minThr,double maxThr,std::vector<float> maskedValues)
{
	//Init
	bwLocation= 0;
	bwScale= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	bool requireFiniteValues= true;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,useRange,minThr,maxThr,maskedValues,requireFiniteValues)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to extract pixel rectangular data!");
		#endif
		return -1;
	}
	
	//Compute median & MAD
	bool useParallelVersion= false;
	float median= Caesar::StatsUtils::GetMedianFast<float>(tile_pixels);
	float medianRMS= Caesar::StatsUtils::GetMADFast(tile_pixels,median,useParallelVersion)*1.4826;
	
	//## Compute biweight robust estimators
	double tol= 0.0001;
	double nmaxIter= 10;
	std::pair<float,float> biweightEstimators= Caesar::StatsUtils::GetBiWeightEstimators<float>(tile_pixels,median,medianRMS,C,tol,nmaxIter);
	bwLocation= biweightEstimators.first;
	bwScale= biweightEstimators.second;
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileBiWeightStats()

//================================================================
//===    BACKGROUND CALCULATION
//================================================================
ImgBkgData* Image::ComputeBkg(int estimator,bool computeLocalBkg,int boxSizeX,int boxSizeY, double gridStepSizeX,double gridStepSizeY,bool use2ndPass,bool skipOutliers,double seedThr,double mergeThr,int minPixels,bool useRange,double minThr,double maxThr)
{	
	//## Compute bkg data
	ImgBkgData* bkgData= BkgFinder::FindBkg(
		this,estimator,
		computeLocalBkg, boxSizeX, boxSizeY, gridStepSizeX,gridStepSizeY,
		use2ndPass,
		skipOutliers,seedThr,mergeThr,minPixels,
		useRange,minThr,maxThr
	);	
	if(!bkgData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Computation of local background failed for this image!");
		#endif
		return 0;
	}
	return bkgData;
	
}//close ComputeBkg()


ImgBkgData* Image::ComputeBkg(ImgBkgPars pars)
{	
	//## Compute bkg data
	ImgBkgData* bkgData= BkgFinder::FindBkg(
		this,pars.estimator,
		pars.computeLocalBkg, pars.boxSizeX, pars.boxSizeY, pars.gridStepSizeX, pars.gridStepSizeY,
		pars.use2ndPass,
		pars.skipOutliers,pars.seedThr,pars.mergeThr,pars.minPixels,
		pars.useRange,pars.minThr,pars.maxThr
	);	
	if(!bkgData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Computation of local background failed for this image!");
		#endif
		return 0;
	}
	return bkgData;
	
}//close ComputeBkg()


Image* Image::GetSignificanceMap(ImgBkgData* bkgData,bool useLocalBkg)
{
	//Check image
	if(!bkgData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to bkg data!");
		#endif
		return 0;
	}

	//Integrity check for local bkg
	long int Nx= this->GetNx();
	long int Ny= this->GetNy();
	if(useLocalBkg){
		if(!bkgData->HasLocalBkg()){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Local bkg option requested but no local bkg data are available!");
			#endif
			return 0;
		}
		if( Nx!=(bkgData->BkgMap)->GetNx() || Ny!=(bkgData->BkgMap)->GetNy() ||
				Nx!=(bkgData->NoiseMap)->GetNx() || Ny!=(bkgData->NoiseMap)->GetNy()
		){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Bkg/Noise maps have different size!");		
			#endif
			return 0;
		}
	}//close if
		
	//Clone this image and reset content
	TString imgName= Form("%s_significance",m_name.c_str());
	Image* significanceMap= this->GetCloned(std::string(imgName),true,true);
	significanceMap->Reset();

	
	#ifdef OPENMP_ENABLED
		Caesar::StatMoments<double> moments_t;	
		std::vector<Caesar::StatMoments<double>> parallel_moments;
	
		//#pragma omp declare reduction (merge : std::vector<Caesar::StatMoments<double>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
		#pragma omp parallel private(moments_t) //reduction(merge: parallel_moments)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			if(useLocalBkg){
				#pragma omp for collapse(2)
				for(long int i=0;i<Nx;i++){		
					for(long int j=0;j<Ny;j++){	
						long int gBin= GetBin(i,j);
						double w= m_pixels[gBin];											
						double bkgLevel= (bkgData->BkgMap)->GetBinContent(i,j);
						double bkgRMS= (bkgData->NoiseMap)->GetBinContent(i,j);
						if( w==0 || bkgRMS<=0) {
							#ifdef LOGGING_ENABLED
								DEBUG_LOG("Empty pixel or invalid bkg  ("<<i<<","<<j<<") skip it!");
							#endif
							continue;
						}
						double Z= (w-bkgLevel)/bkgRMS;
						significanceMap->FillPixelMT(moments_t,i,j,Z);
					}//end loop
				}//end loop 

				//parallel_moments.push_back(moments_t);
				parallel_moments[thread_id]= moments_t;
			}//close if local bkg

			else{//global bkg
				double bkgLevel= bkgData->gBkg;
				double bkgRMS= bkgData->gNoise;
						
				#pragma omp for collapse(2)
				for(long int i=0;i<Nx;i++){		
					for(long int j=0;j<Ny;j++){	
						long int gBin= GetBin(i,j);
						double w= m_pixels[gBin];										
						if( w==0 || bkgRMS<=0) continue;
				
						double Z= (w-bkgLevel)/bkgRMS;
						significanceMap->FillPixelMT(moments_t,i,j,Z);
					}//end loop
				}//end loop

				//parallel_moments.push_back(moments_t);
				parallel_moments[thread_id]= moments_t;

			}//close else 
			
		}//close parallel section

		//Update moments from parallel estimates
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(m_StatMoments,parallel_moments)<0){	
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");	
			#endif
			if(significanceMap){
				delete significanceMap;
				significanceMap= 0;
			}
			return nullptr;
		}

	#else
			
		if(useLocalBkg){
			for(long int i=0;i<Nx;i++){		
				for(long int j=0;j<Ny;j++){	
					long int gBin= GetBin(i,j);
					double w= m_pixels[gBin];											
					double bkgLevel= (bkgData->BkgMap)->GetBinContent(i,j);
					double bkgRMS= (bkgData->NoiseMap)->GetBinContent(i,j);
					if( w==0 || bkgRMS<=0) {	
						#ifdef LOGGING_ENABLED
							DEBUG_LOG("Empty pixel or invalid bkg  ("<<i<<","<<j<<") skip it!");
						#endif
						continue;
					}
					double Z= (w-bkgLevel)/bkgRMS;
					significanceMap->FillPixel(i,j,Z);
				}//end loop
			}//end loop 
		}//close if useLocalBkg
	
		else{
			double bkgLevel= bkgData->gBkg;
			double bkgRMS= bkgData->gNoise;
			for(long int i=0;i<Nx;i++){		
				for(long int j=0;j<Ny;j++){	
					long int gBin= GetBin(i,j);
					double w= m_pixels[gBin];										
					if( w==0 || bkgRMS<=0) continue;
				
					double Z= (w-bkgLevel)/bkgRMS;
					significanceMap->FillPixel(i,j,Z);
				}//end loop
			}//end loop
		}//close else

	#endif

	return significanceMap;

}//close GetSignificanceMap()

//=======================================================
//==          SOURCE EXTRACTION
//=======================================================
int Image::FindCompactSource(std::vector<Source*>& sources,double thr,int minPixels){

	//Find sources by simple thresholding using the same image as significance map and no bkgdata
	bool findNestedSources= false;
	if(this->FindCompactSource(sources,this,0,thr,thr,minPixels,findNestedSources)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Compact source finder failed!");	
		#endif
		return -1;
	}
	return 0;

}//close FindCompactSource()


int Image::FindCompactSource(std::vector<Source*>& sources,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNestedSources,Image* blobMask,double minNestedMotherDist,double maxMatchingPixFraction,long int nPixThrToSearchNested)
{
	//Check blob mask (if nested source search is required)
	if(findNestedSources && !blobMask){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to blob mask given (hint: you must provide a binary blob mask to search nested sources!");
		#endif
		return -1;
	}

	//Find sources
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	int status= BlobFinder::FindBlobs(this,sources,floodImg,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Blob finder failed!");
		#endif
		CodeUtils::DeletePtrCollection<Source>(sources);
		return -1;
	}

	//Find nested sources?
	if(findNestedSources && sources.size()>0){
		int status_nested= FindNestedSource(sources,blobMask,bkgData,minPixels,minNestedMotherDist,maxMatchingPixFraction,nPixThrToSearchNested);
		if(status_nested<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Nested source search failed!");
			#endif
		}
	}//close if

	return 0;

}//close FindCompactSource()



int Image::FindNestedSource(std::vector<Source*>& sources,Image* blobMask,ImgBkgData* bkgData,int minPixels,double minNestedMotherDist,double maxMatchingPixFraction,long int nPixThrToSearchNested)
{
	//Check if given mother source list is empty
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty source list given!");
		#endif
		return 0;
	}

	//Check blob mask given
	if(!blobMask){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to blob mask given!");
		#endif
		return -1;
	}

	//Find image mask of found sources
	Image* sourceMask= this->GetSourceMask(sources,false);
	if(!sourceMask){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to computed source mask!");
		#endif
		return -1;
	}

	//Find blob+source mask
	Image* sourcePlusBlobMask= sourceMask->GetMask(blobMask,true);
	if(!sourcePlusBlobMask){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute (source+blob) mask!");
		#endif
		CodeUtils::DeletePtrCollection<Image>({sourceMask});	
		return -1;
	}

	//Find nested blobs 
	std::vector<Source*> NestedSources;
	double fgValue= 1;
	int status= BlobFinder::FindBlobs(this,NestedSources,sourcePlusBlobMask,bkgData,fgValue,fgValue,minPixels,false,false);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Nested blob finder failed!");
		#endif
		CodeUtils::DeletePtrCollection<Image>({sourceMask,sourcePlusBlobMask});	
		CodeUtils::DeletePtrCollection<Source>(NestedSources);	
		return -1;
	}

	//Clear maps
	CodeUtils::DeletePtrCollection<Image>({sourceMask,sourcePlusBlobMask});	

	//Add nested sources to mother source
	int nNestedSources= static_cast<int>(NestedSources.size());
	
	if(nNestedSources>=0){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("#"<<nNestedSources<<" blobs found in curvature map ...");
		#endif

		//## Init mother-nested association list
		std::vector<int> nestedSourcesToBeRemoved;
		std::vector< std::vector<int> > MotherNestedAssociationList;
		for(int i=0;i<nSources;i++) MotherNestedAssociationList.push_back( std::vector<int>() );

		//## Find matching between mother and nested sources
		for(int j=0;j<nNestedSources;j++){
			bool isMotherFound= false;
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Finding matching for nested source no. "<<j);
			#endif

			for(int i=0;i<nSources;i++){
				int sourceId= sources[i]->Id;
				bool isOverlapping= NestedSources[j]->IsAdjacentSource(sources[i]);
				
				if(isOverlapping){
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Nested source no. "<<j<<" added to source id="<<sourceId<<" ...");
					#endif
					MotherNestedAssociationList[i].push_back(j);
					isMotherFound= true;
					break;
				}
				
			}//end loop mother sources

			//If nested is not associated to any mother source, mark for removal
			if(!isMotherFound){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Cannot find mother source for nested source no. "<<j<<", will remove it from the list of nested sources...");
				#endif
				nestedSourcesToBeRemoved.push_back(j);
			}			

		}//end loop nested sources	

		//## Select nested
		int nSelNestedSources= 0;
		for(size_t i=0;i<MotherNestedAssociationList.size();i++){
			long int NPix= sources[i]->GetNPixels(); 
			double beamArea= sources[i]->GetBeamFluxIntegral();
			bool isMotherSourceSplittableInNestedComponents= (NPix>nPixThrToSearchNested);
			int nComponents= static_cast<int>(MotherNestedAssociationList[i].size());
			if(nComponents<=0 || !isMotherSourceSplittableInNestedComponents) continue;

			//Get mother source stats
			if(!sources[i]->HasStats()){
				sources[i]->ComputeStats();
			}
			
			//If only one component is present select it if:
			//  1) mother and nested distance is > thr (e.g. 
			//  2) mother and nested pix superposition is <thr (e.g. 50%)
			for(int j=0;j<nComponents;j++){
				int nestedIndex= MotherNestedAssociationList[i][j];

				//Set pars
				NestedSources[nestedIndex]->SetBeamFluxIntegral(beamArea);

				//Compute nested source stats & pars
				NestedSources[nestedIndex]->ComputeStats();
				NestedSources[nestedIndex]->ComputeMorphologyParams();

				if(nComponents==1){
					//Compute centroid distances
					float centroidDistX= fabs(sources[i]->X0-NestedSources[nestedIndex]->X0);
					float centroidDistY= fabs(sources[i]->Y0-NestedSources[nestedIndex]->Y0);
					
					//Compute nmatching pixels
					long int nMatchingPixels= sources[i]->GetNMatchingPixels(NestedSources[nestedIndex]);
					float matchingPixFraction= (float)(nMatchingPixels)/(float)(NPix);

					//Select nested?
					bool areOffset= (centroidDistX>minNestedMotherDist || centroidDistY>minNestedMotherDist);
					bool isNestedSmaller= (matchingPixFraction<maxMatchingPixFraction);
					
					if( areOffset || isNestedSmaller){//Add nested to mother source
						#ifdef LOGGING_ENABLED
							DEBUG_LOG("Adding nested source to mother source (name="<<sources[i]->GetName()<<", id="<<sources[i]->Id<<", nPix="<<NPix<<"): areOffset? "<<areOffset<<" (dist_x="<<centroidDistX<<", dist_y="<<centroidDistY<<"), isNestedSmaller?"<<isNestedSmaller<<" (matchingPixFraction="<<matchingPixFraction<<", maxMatchingPixFraction="<<maxMatchingPixFraction<<", nPixThrToSearchNested="<<nPixThrToSearchNested<<")");					
						#endif

						sources[i]->AddNestedSource(NestedSources[nestedIndex]);
						nSelNestedSources++;
					}
					else{//do not select nested!
						nestedSourcesToBeRemoved.push_back(nestedIndex);
					}
				}//close if
				else{
					//Add nested to mother source
					sources[i]->AddNestedSource(NestedSources[nestedIndex]);
					nSelNestedSources++;
				}
			}//end loop nested components
		}//end loop mother sources

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("#"<<nSelNestedSources<<" nested sources found and added to mother sources ...");
		#endif

		//Delete nested source selected for removal
		for(size_t k=0;k<nestedSourcesToBeRemoved.size();k++){
			int nestedSourceIndex= nestedSourcesToBeRemoved[k];
			if(NestedSources[nestedSourceIndex]){
				delete NestedSources[nestedSourceIndex];
				NestedSources[nestedSourceIndex]= 0;
			}
		}
		NestedSources.clear();
						
	}//close nNestedBlobs>0

	return 0;

}//close FindNestedSources()


int Image::FindBlendedSources(std::vector<Source*>& deblendedSources,std::vector<ImgPeak>& deblendedPeaks,double sigmaMin,double sigmaMax,double sigmaStep,int minBlobSize,double thrFactor,int kernelFactor)
{
	//Find blended sources
	int status= BlobFinder::FindBlendedBlobs(
		deblendedSources,deblendedPeaks,
		this,
		sigmaMin,sigmaMax,sigmaStep,
		minBlobSize,
		thrFactor,kernelFactor
	);

	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to find blended sources inside image!");
		#endif
		return -1;
	}

	return 0;

}//close FindBlendedSources()


int Image::FindExtendedSource_CV(std::vector<Source*>& sources,Image* initSegmImg,ImgBkgData* bkgData,int minPixels,bool findNegativeExcess,double dt,double h,double lambda1,double lambda2,double mu,double nu,double p,int niters,double tol,int niters_inner,int niters_reinit)
{
	//## Compute segmented image
	Image* segmentedImg= ChanVeseSegmenter::FindSegmentation(this,initSegmImg,false,dt,h,lambda1,lambda2,mu,nu,p,niters,tol,niters_inner,niters_reinit);
	if(!segmentedImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute ChanVese image segmentation!");
		#endif
		return -1;
	}
	
	//## Finding blobs in masked image
	double fgValue= 1;	
	bool findNestedSources= false;
	int status= this->FindCompactSource(sources,segmentedImg,bkgData,fgValue,fgValue,minPixels,findNestedSources);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Finding sources in Chan-Vese segmented mask failed!");
		#endif
		return -1;
	}
	
	//## Clear segmented image
	if(segmentedImg){
		delete segmentedImg;
		segmentedImg= 0;
	}

	//## Remove sources of negative excess 
	if(!findNegativeExcess){

		//Compute image stats if not already present
		if(!this->HasStats()){
			this->ComputeStats(true,false,false);
		}	
		if(!m_Stats){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute image stats, clearing all detected sources and returning empty list!");
			#endif
			for(size_t k=0;k<sources.size();k++){
				if(sources[k]){
					delete sources[k];
					sources[k]= 0;
				}	
			}//end loop sources
			sources.clear();
			return -1;
		}//close if

		//Find and remove sources with flux below the input image median (THIS COULD BE IMPROVED!)
		std::vector<size_t> sourcesToBeRemoved;		
		int imgMedian= m_Stats->median;
		for(size_t k=0;k<sources.size();k++){
			//Tag sources as extended
			sources[k]->SetType(eExtended);

			//Check if source is a "negative excess"
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeletePtrItems(sources, sourcesToBeRemoved);
	}//close if

	return 0;
	
}//close FindExtendedSources_CV()

//==================================================
//==       FILTERING METHODS
//==================================================

Image* Image::GetMask(Image* mask,bool isBinary)
{

	//## Check input mask
	if(!mask) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image mask!");
		#endif
		return 0;
	}
		
	//## Check mask bins
	long int Nx= mask->GetNx();
	long int Ny= mask->GetNy();
	if(Nx!=this->GetNx() || Ny!=this->GetNy()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Mask binning is different than current image!");
		#endif
		return 0;
	}	

	//## Clone map
	TString imgName= Form("%s_Mask",m_name.c_str());	
	Image* maskedImage= this->GetCloned(std::string(imgName),true,true);
	maskedImage->Reset();

	//## Loop over mask	
	if(isBinary){
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2)
		#endif
		for(long int i=0;i<Nx;i++){	
			for(long int j=0;j<Ny;j++){
				double binContent= this->GetPixelValue(i,j);
				if(binContent==0) continue;
				double maskContent= mask->GetPixelValue(i,j);
				if(maskContent!=0) maskedImage->SetPixelValue(i,j,1);
			}//end loop bins Y
		}//end loop bins X
	}//close if
	else{
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2)
		#endif
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				double binContent= this->GetPixelValue(i,j);	
				if(binContent==0) continue;
				double maskContent= mask->GetPixelValue(i,j);
				if(maskContent!=0) maskedImage->SetPixelValue(i,j,binContent);
			}//end loop bins Y
		}//end loop bins X
	}//close else

	//Force re-computation of stats (in parallel computation moments are wrong)
	maskedImage->ComputeStats(true,false,true);

	return maskedImage;

}//close GetMask()


int Image::SubtractFittedSources(const std::vector<Source*>& sources,bool subtractNested,bool recomputeStats,double nsigmaTrunc)
{
	//Nothing to be done is source collection is empty
	if(sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Given source collection is empty, nothing to be done!");	
		#endif
		return 0;
	}

	//Get source fit model image
	Image* modelImg= GetSourceFitModelImage(sources,subtractNested,nsigmaTrunc);
	if(!modelImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute source fit model image!");
		#endif
		return -1;
	}
	
	//Subtract model image from current image
	if(this->Add(modelImg,-1,recomputeStats)<0){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to subtract model image from this image!");
		#endif
		delete modelImg;
		modelImg= 0;
		return -1;
	}

	//Clear data
	if(modelImg){
		delete modelImg;
		modelImg= 0;
	}

	return 0;

}//close SubtractFittedSources()

int Image::SubtractFittedSource(Source* source,bool subtractNested,bool recomputeStats,double nsigmaTrunc)
{
	//Nothing to be done is source collection is empty
	if(!source){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input source given!");	
		#endif
		return -1;
	}

	//Get source fit model image	
	std::vector<Source*> sources {source};
	Image* modelImg= GetSourceFitModelImage(sources,subtractNested,nsigmaTrunc);
	if(!modelImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute source fit model image!");
		#endif
		return -1;
	}
	
	//Subtract model image from current image
	if(this->Add(modelImg,-1,recomputeStats)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to subtract model image from this image!");
		#endif
		delete modelImg;
		modelImg= 0;
		return -1;
	}

	//Clear data
	if(modelImg){
		delete modelImg;
		modelImg= 0;
	}

	return 0;

}//close SubtractFittedSource()


int Image::FillFromSourceFitModel(Source* source,bool useNested,double nsigmaTrunc)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to source given!");
		#endif
		return -1;
	}

	//Check if has fit info
	bool hasFitInfo= source->HasFitInfo();
	if(hasFitInfo) {
		//Loop over source components and fill image
		SourceFitPars fitPars= source->GetFitPars();
		for(int k=0;k<fitPars.GetNComponents();k++){
			double X0= fitPars.GetParValue(k,"x0");	
			double Y0= fitPars.GetParValue(k,"y0");
			double A= fitPars.GetParValue(k,"A");	
			double sigmaX= fitPars.GetParValue(k,"sigmaX");
			double sigmaY= fitPars.GetParValue(k,"sigmaY");
			double theta= fitPars.GetParValue(k,"theta")*TMath::DegToRad();
			long int boxHalfWidthX= std::ceil(nsigmaTrunc*sigmaX);
			long int boxHalfWidthY= std::ceil(nsigmaTrunc*sigmaY);
			long int binId_centroid= this->FindBin(X0,Y0);
			if(binId_centroid<0) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Source fitted centroid ("<<X0<<","<<Y0<<") outside image range, skip source...");
				#endif
				continue;
			}
		
			long int ix_centroid= GetBinX(binId_centroid);
			long int iy_centroid= GetBinY(binId_centroid);

			//Loop over model box
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Source "<<source->GetName()<<", component "<<k+1<<": A="<<A<<", pos("<<X0<<","<<Y0<<"), sigma("<<sigmaX<<","<<sigmaY<<"), theta="<<theta<<", box half widths=("<<boxHalfWidthX<<","<<boxHalfWidthY<<"), centroid bin("<<ix_centroid<<","<<iy_centroid<<")");
			#endif

			for(long int ix=-boxHalfWidthX;ix<=boxHalfWidthX;ix++){
				for(long int iy=-boxHalfWidthY;iy<=boxHalfWidthY;iy++){
					//Check bin exists
					long int ix_pix= ix_centroid + ix;
					long int iy_pix= iy_centroid + iy;
					if(!this->HasBin(ix_pix,iy_pix)) continue;
					double x_pix= this->GetX(ix_pix);
					double y_pix= this->GetY(iy_pix);

					//Evaluate source fit model in bin center
					double w= MathUtils::EvalGaus2D(x_pix,y_pix,A,X0,Y0,sigmaX,sigmaY,theta);

					//Add to model bin
					double w_old= this->GetBinContent(ix_pix,iy_pix);
					this->SetBinContent(ix_pix,iy_pix,w+w_old);
				}//end loop y box
			}//end loop x box		
		}//end loop components
	}//close if
	else{	
		if(useNested){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Source "<<source->GetName()<<" without fit info, check nested sources...");	
			#endif
			std::vector<Source*> nestedSources= source->GetNestedSources();
			for(size_t k=0;k<nestedSources.size();k++){
				this->FillFromSourceFitModel(nestedSources[k],nsigmaTrunc);
			}
		}
	}//close else 	

	return 0;

}//close FillFromSourceFitModel()

Image* Image::GetSourceFitModelImage(const std::vector<Source*>& sources,bool useNested,double nsigmaTrunc)
{
	// Clone map
	Image* modelImg= this->GetCloned("",true,true);
	modelImg->Reset();
	
	// Check if source collection is empty
	if(sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source collection is empty, no model can be created, returning empty image!");	
		#endif
		return modelImg;
	}

	//Loop over sources and fill from source fit model (including nested sources)
	for(size_t i=0;i<sources.size();i++){		
		modelImg->FillFromSourceFitModel(sources[i],useNested,nsigmaTrunc);
	}

	return modelImg;

}//close GetSourceFitModelImage()

int Image::MaskSources(std::vector<Source*>const& sources,float maskValue)
{
	//## Check source list
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source list is empty, nothing to be done...");
		#endif
		return 0;	
	}

	//## Mask sources 
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(int k=0;k<nSources;k++){
		for(int l=0;l<sources[k]->GetNPixels();l++){
			long int id= (sources[k]->GetPixel(l))->id;
			this->SetPixelValue(id,maskValue);
		}//end loop pixels
	}//end loop sources		
		
	//Force re-computation of stats after masks
	bool computeRobustStats= true;
	bool useRange= false;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	std::vector<float> maskedValues {0,maskValue};
	bool forceRecomputing= true;
	bool useParallelVersion= false;
	this->ComputeStats(computeRobustStats,forceRecomputing,useRange,minThr,maxThr,useParallelVersion,maskedValues);

	return 0;

}//close MaskSources()


int Image::GetBkgInfoAroundSource(BkgSampleData& bkgSampleData,Source* source,int boxThickness,int bkgEstimator,Image* mask,bool useParallelVersion,std::vector<float> maskedValues)
{
	//Check source
	if(!source){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input source given!");
		#endif
		return -1;
	}

	//Check box size
	if(boxThickness<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Box size must be >0!");	
		#endif
		return -1;
	}

	//Get source bounding box
	long int ix_min= -1;
	long int ix_max= -1;
	long int iy_min= -1;
	long int iy_max= -1;
	source->GetSourcePixelRange(ix_min,ix_max,iy_min,iy_max);

	if(ix_min<0 || ix_min>=ix_max || ix_max>=m_Nx){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid source x range ("<<ix_min<<","<<ix_max<<")!");
		#endif
		return -1;
	}
	if(iy_min<0 || iy_min>=iy_max || iy_max>=m_Ny){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid source y range ("<<iy_min<<","<<iy_max<<")!");
		#endif
		return -1;
	}

	//Exclude pixels belonging to source
	long int nPixels= this->GetNPixels();
	std::vector<bool> use_pixels(nPixels,true);

	for(int l=0;l<source->GetNPixels();l++){
		long int id= (source->GetPixel(l))->id;
		if(!this->HasBin(id)) continue;
		use_pixels[id]= false;
	}//end loop pixels
	
	
	//Select pixels to compute bkg info
	std::vector<double> pixels;
	for(long int i=ix_min-boxThickness;i<=ix_max+boxThickness;i++){
		for(long int j=iy_min-boxThickness;j<=iy_max+boxThickness;j++){
			long int id= this->GetBin(i,j);
			if( id<0 || !use_pixels[id] || (mask && mask->GetPixelValue(id)<=0) ) continue;
			double S= this->GetPixelValue(id);	
			bool isFiniteValue= std::isfinite(S);
			bool isMaskedValue= false;
			for(size_t k=0;k<maskedValues.size();k++){
				if(S==maskedValues[k]){
					isMaskedValue= true;
					break;
				}
			}
			bool isGoodValue= (isFiniteValue && !isMaskedValue);
			if(isGoodValue) pixels.push_back(S);
		}//end loop y
	}//end loop x

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<pixels.size()<<" pixels selected to measure box around source "<<source->GetName()<<" ...");
	#endif

	//Compute bkg info
	double bkgLevel= 0;
	double bkgRMS= 0;
	if(bkgEstimator==eMedianBkg)
	{
		//Compute median
		double median= Caesar::StatsUtils::GetMedianFast<double>(pixels,useParallelVersion);
		
		//Compute MAD = median(|x_i-median|)
		double medianMAD= Caesar::StatsUtils::GetMADFast(pixels,median);	
		double medianRMS= medianMAD*1.4826;//0.6744888;

		bkgLevel= median;
		bkgRMS= medianRMS;
	}
	else if(bkgEstimator==eMedianClippedBkg)
	{
		//Compute mean & rms
		double mean= 0;
		double rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS<double>(mean,rms,pixels);
		
		//Compute median
		double median= Caesar::StatsUtils::GetMedianFast<double>(pixels,useParallelVersion);
		
		//Compute clipped estimators
		double clipSigma= 3;
		int clipMaxIter= 100;
		double clipTolerance= 0.1;
		ClippedStats<double> clipped_stats;
		Caesar::StatsUtils::GetClippedEstimators(clipped_stats,pixels,median,mean,rms,clipSigma,clipMaxIter,clipTolerance,useParallelVersion);
		
		bkgLevel= clipped_stats.median;
		bkgRMS= clipped_stats.stddev;
	}
	else if(bkgEstimator==eMeanBkg)
	{
		//Compute mean & rms
		double mean= 0;
		double rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS<double>(mean,rms,pixels);
		
		bkgLevel= mean;
		bkgRMS= rms;
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid bkg estimator ("<<bkgEstimator<<") given!");
		#endif
		return -1;
	}

	bkgSampleData.ix_min= ix_min-boxThickness;
	bkgSampleData.ix_max= ix_max+boxThickness;
	bkgSampleData.iy_min= iy_min-boxThickness; 
	bkgSampleData.iy_max= iy_max+boxThickness;
	bkgSampleData.npix= pixels.size();
	bkgSampleData.isReliable= true;
	bkgSampleData.bkgLevel= bkgLevel;
	bkgSampleData.bkgRMS= bkgRMS;

	return 0;

}//close GetBkgInfoAroundSource()


Image* Image::GetSourceMask(std::vector<Source*>const& sources,bool isBinary,bool invert,bool searchSourceCoords)
{
	//## Clone map
	long int Nx= this->GetNx();
	long int Ny= this->GetNy();
	bool copyMetaData= true;
	bool resetStats= true;
	TString imgName= Form("%s_SourceMask",m_name.c_str());
	//Image* maskedImage= this->GetCloned(std::string(imgName),copyMetaData,resetStats);
	Image* maskedImage= this->GetCloned("",copyMetaData,resetStats);
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("pixel size="<<maskedImage->GetPixelDataSize()<<" (original size="<<this->GetPixelDataSize()<<")");
	#endif

	//## Check source list
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source list is empty, returning same image!");
		#endif
		return maskedImage;	
	}

	//cout<<"Image::GetSourceMask(): pto 1"<<endl;

	//## Search source pixel coordinates in image?
	//## If sources are extracted from the same map there is no need to do that
	//## however for sources extracted from submaps or tiles it is needed
	std::vector<long int> masked_pix_ids;
	if(searchSourceCoords) {
		
		for(int k=0;k<nSources;k++){
			for(int l=0;l<sources[k]->GetNPixels();l++){
				Pixel* pixel= sources[k]->GetPixel(l);
				double x= pixel->x;
				double y= pixel->y;
				long int id= this->FindBin(x,y);
				if(id<0) continue;
				masked_pix_ids.push_back(id);	
			}
		}
	}//close if
	else{
		
		for(int k=0;k<nSources;k++){
			for(int l=0;l<sources[k]->GetNPixels();l++){
				Pixel* pixel= sources[k]->GetPixel(l);
				long int id= pixel->id;
				masked_pix_ids.push_back(id);	
			}
		}
	}//close else

	
	//## Create mask
	if(invert)
	{
		if(isBinary)
		{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t i=0;i<masked_pix_ids.size();i++){
				long int id= masked_pix_ids[i];
				maskedImage->SetPixelValue(id,0);
			}
			
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(long int i=0;i<Nx;i++){
				for(long int j=0;j<Ny;j++){
					double w= maskedImage->GetPixelValue(i,j);		
					if(w==0) continue;
					maskedImage->SetPixelValue(i,j,1);				
				}
			}
		}//close if
		else{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t i=0;i<masked_pix_ids.size();i++){
				long int id= masked_pix_ids[i];
				maskedImage->SetPixelValue(id,0);
			}
		}//close else

		//Force re-computation of stats
		//maskedImage->ComputeStats(true,false,true);

	}//close if invert
	else	
	{
		//Reset map and loop over sources
		maskedImage->Reset();

		if(isBinary)
		{	
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif	
			for(size_t i=0;i<masked_pix_ids.size();i++){
				long int id= masked_pix_ids[i];
				maskedImage->SetPixelValue(id,1);
			}	
	
		}//close if
		else
		{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t i=0;i<masked_pix_ids.size();i++){
				long int id= masked_pix_ids[i];
				double w= this->GetPixelValue(id);
				maskedImage->SetPixelValue(id,w);
			}
		}//close else

		//Force re-computation of stats
		//maskedImage->ComputeStats(true,false,true);

	}//close else

	//# Force re-computation of stats
	double minThr=-std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	bool useRange= true;
	std::vector<float> maskedValues {};
	if(!isBinary) maskedValues= {0};
	bool computeRobustStats= true;
	//bool forceRecomputing= false;
	bool forceRecomputing= true;
	bool useParallelVersion= false;

	maskedImage->ComputeStats(computeRobustStats,forceRecomputing,useRange,minThr,maxThr,useParallelVersion,maskedValues);

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("(after) pixel size="<<maskedImage->GetPixelDataSize()<<" (original size="<<this->GetPixelDataSize()<<")");
	#endif

	return maskedImage;

}//close GetSourceMask()


Image* Image::GetSourceResidual(std::vector<Source*>const& sources,int kernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr,double zBrightThr,int psSubtractionMethod,Image* mask,int bkgBoxThickness)
{
	//Check bkg data
	if(dilateModel==eDilateWithBkg){
		/*
	 	if(!bkgData){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Selected to use bkg dilation but null ptr to bkg data!");	
			#endif
			return nullptr;
		}
		*/
		if(useLocalBkg && bkgData && !bkgData->HasLocalBkg()){		
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Selected to use local bkg but no local bkg data are available!");
			#endif
			return nullptr;
		}
	}//close if

	
	//Check thresholds
	if(zBrightThr<zThr){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid z threshold given (hint: zBrightThr shall be >= than zThr)!");
		#endif
		return nullptr;
	}

	//Clone input image	
	TString imgName= Form("%s_Residual",m_name.c_str());	
	Image* residualImg= this->GetCloned(std::string(imgName),true,true);
	if(!residualImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to clone input image, returning nullptr!");
		#endif
		return nullptr;
	}

	//Check source list
	if(sources.size()<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source list empty, nothing to be dilated!");
		#endif
		return residualImg;
	}


	//Fill list of sources to be subtracted
	std::vector<Source*> dilatedSources;
	std::vector<Source*> psSources;
	for(size_t i=0;i<sources.size();i++){
		Source* source= sources[i];
		int sourceType= source->Type;
		bool hasFitInfo= source->HasFitInfo();
		bool hasNestedSources= source->HasNestedSources();
		bool isSelectedSource= (dilateSourceType==-1 || sourceType==dilateSourceType);
		bool addToPSList= (
			hasFitInfo && 
			psSubtractionMethod==ePS_MODELSUBTRACTION &&
			(dilateSourceType==-1 || dilateSourceType==ePointLike || dilateSourceType==eCompact)
		);
		
		//Check is source stats are available
		bool hasStats= source->HasStats();
		if(!hasStats){
			#ifdef LOGGING_ENABLED
				WARN_LOG("No stats computed for input source "<<source->GetName()<<"...computing!");
			#endif
			source->ComputeStats(true,true);
		}

		//Check if bright source
		Pixel* seedPixel= source->GetSeedPixel();
		if(!seedPixel){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to find seed pixel for source "<<source->GetName()<<", skip to next!");
			#endif
			continue;
		}
		std::pair<double,double> bkgInfo= seedPixel->GetBkg();
		double bkgLevel= bkgInfo.first;
		double noiseLevel= bkgInfo.second;
		double S= seedPixel->S;
		double Z= (S-bkgLevel)/noiseLevel;
		bool isDilatedSource= false;
		bool isBrightSource= false;
		if(fabs(Z)>zThr) isDilatedSource= true;
		if(fabs(Z)>zBrightThr) isBrightSource= true;
		
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Source "<<source->GetName()<<" sourceType="<<sourceType<<", hasFitInfo?"<<hasFitInfo<<", hasNestedSources?"<<hasNestedSources<<", isSelectedSource?"<<isSelectedSource<<", addToPSList? "<<addToPSList<<" (psSubtractionMethod="<<psSubtractionMethod<<", dilateSourceType="<<dilateSourceType<<"), isDilatedSource? "<<isDilatedSource<<", isBrightSource? "<<isBrightSource);
		#endif

		//Select source
		if(!isDilatedSource) continue;

		//If bright source add to dilate list (if not a fitted point-source), otherwise check nested
		if(isBrightSource){
			if(addToPSList) psSources.push_back(source);
			else dilatedSources.push_back(source);
		}
		else{
			if(hasNestedSources && skipToNested){
				std::vector<Source*> nestedSources= source->GetNestedSources();
				for(size_t k=0;k<nestedSources.size();k++){
					int nestedSourceType= nestedSources[k]->Type;
					bool hasFitInfo_nested= source->HasFitInfo();
					bool isSelectedSource_nested= (dilateSourceType==-1 || nestedSourceType==sourceType);
					if(!isSelectedSource_nested) continue;
					bool addToPSList_nested= (
						hasFitInfo_nested && 
						psSubtractionMethod==ePS_MODELSUBTRACTION &&
						(dilateSourceType==-1 || dilateSourceType==ePointLike || dilateSourceType==eCompact)
					);
					
					if(addToPSList_nested) psSources.push_back(nestedSources[k]);
					else dilatedSources.push_back(nestedSources[k]);

				}//end loop nested sources
			}//close if has nested
			else{
				if(!isSelectedSource) continue;
				if(addToPSList) psSources.push_back(source);
				else dilatedSources.push_back(source);			
			}//close else mother source
		}//close else bright source
	}//end loop sources

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<psSources.size()<<" point-sources selected for subtraction...");
	#endif

	//Subtract point-source using fit information?
	if( psSubtractionMethod==ePS_MODELSUBTRACTION ){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Subtracting "<<psSources.size()<<" point-source from input map...");
		#endif
		bool subtractNested= false;
		bool recomputeStats= false;
		if(residualImg->SubtractFittedSources(psSources,subtractNested,recomputeStats)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to subtract point-sources from input map!");
			#endif
			if(residualImg) {
				delete residualImg;
				residualImg= 0;
			}
			return nullptr;
		}
	}//close if point source subtraction

	//Dilate sources from input map (consider only non-point source if PS_MODELSUBTRACTION flag is used)
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<dilatedSources.size()<<" sources to be dilated from input map...");
	#endif
	for(size_t i=0;i<dilatedSources.size();i++){
		//Dilate source from image
		int status= MorphFilter::DilateAroundSource(
			residualImg,dilatedSources[i],
			kernSize,dilateModel,
			bkgData,useLocalBkg,
			randomize,
			mask,bkgBoxThickness
		);

		if(status<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to dilate source "<<dilatedSources[i]->GetName()<<"!");
			#endif
			if(residualImg) {
				delete residualImg;
				residualImg= 0;
			}
			return nullptr;		
		}//close if
	}//end loop dilated sources

	//Re-compute stats of residual image
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Recomputing stats of residual image...");
	#endif
	residualImg->ComputeStats(true,true,false);

	return residualImg;

}//close GetSourceResidual()


Image* Image::GetNormalizedImage(std::string normScale,int normmin,int normmax,bool skipEmptyBins)
{
	//Check if image has data
	if(m_pixels.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Image has no data stored!");
		#endif
		return nullptr;
	}

	//Get image content min/max
	double wmin= (this->m_StatMoments).minVal;
	double wmax= (this->m_StatMoments).maxVal;
	if(TMath::IsNaN(wmin) || fabs(wmin)==TMath::Infinity() || TMath::IsNaN(wmax) || fabs(wmax)==TMath::Infinity()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Min ("<<wmin<<") or max ("<<wmax<<") image values are nan or inf, this should not occur (hint: check moment calculation!)");
		#endif
		return nullptr;
	}
	
	//Create normalized image
	TString imgName= Form("%s_Normalized",m_name.c_str());	
	Image* norm_img= this->GetCloned(std::string(imgName),true,true);
	norm_img->Reset();
	
	//Fill norm image
	if(normScale=="LINEAR"){

		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(size_t i=0;i<m_pixels.size();i++){
			double w= m_pixels[i];
			if(skipEmptyBins && w==0) continue;
			double w_norm= normmin + (normmax-normmin)*(w-wmin)/(wmax-wmin);
			norm_img->SetPixelValue(i,w_norm);
		}
	}//close if

	else if(normScale=="LOG"){
		double safemin= 1;
		double safemax= 256;

		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(size_t i=0;i<m_pixels.size();i++){
			double w= m_pixels[i];
			if(skipEmptyBins && w==0) continue;
			double w_norm= safemin + (safemax-safemin)*(w-wmin)/(wmax-wmin);
			double w_log= normmin + log10(w_norm/safemin)/log10(safemax/safemin) * (normmax-normmin);
			norm_img->SetPixelValue(i,w_log);
		}

	}//close else if
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid norm scale option selected ("<<normScale<<") no transform applied to original image!");
		#endif
		return norm_img;
	}

	//Force recomputation of stats if present, otherwise recompute only moments
	bool useRange= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	std::vector<float> maskedValues= {0};
	if(this->HasStats()) status= norm_img->ComputeStats(computeRobustStats,forceRecomputing,useRange);
	else status= norm_img->ComputeMoments(useRange,minThr,maxThr,maskedValues);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to re-compute moments/stats for normalized image!");
		#endif
	}
	
	return norm_img;

}//close GetNormalizedImage()

Image* Image::GetGradientImage(bool invert)
{
	//Compute gradient filtered image
	Image* gradImg= GradientFilter::GetGradientFilter(this);
	if(invert) gradImg->Scale(-1);

	return gradImg;

}//close GetGradientImage()

Image* Image::GetKirschImage()
{
	
	//Compute Kirsh filtered image
	Image* kirschImg= KirschFilter::GetKirschFilter(this);

	return kirschImg;

}//close GetKirschImage()

Image* Image::GetLoGImage(bool invert)
{
	//Compute LoG filtered image
	Image* LoGImg= LoGFilter::GetLoGFilter(this);
	if(LoGImg && invert) LoGImg->Scale(-1);
	
	return LoGImg;

}//close GetLoGImage()

Image* Image::GetNormLoGImage(int size,double scale,bool invert)
{
	//Compute normalized LoG filtered image
	Image* normLoGImg= LoGFilter::GetNormLoGFilter(this,size,scale);
	if(normLoGImg && invert) normLoGImg->Scale(-1);

	return normLoGImg;

}//close GetNormLoGImage()

Image* Image::GetBeamConvolvedImage(double bmaj,double bmin,double bpa,int nsigmas,double scale)
{
	//Compute image convolved with an elliptical gaussian beam
	Image* convImg= GausFilter::GetGausFilter(this,bmaj,bmin,bpa,nsigmas,scale);
	return convImg;

}//close GetBeamConvolvedImage()

Image* Image::GetLaplacianImage(bool invert)
{

	//Compute laplacian filtered image
	Image* laplImg= GradientFilter::GetLaplaceFilter(this);
	if(invert) laplImg->Scale(-1);
	
	return laplImg;

}//close GetLaplacianImage()


Image* Image::GetGuidedFilterImage(int radius,double eps)
{

	//Normalize img
	Image* img_norm= this->GetNormalizedImage("LINEAR",1,256);
	img_norm->SetName("tmpImg");
	if(!img_norm) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get normalized image!");
		#endif
		return 0;
	}

	//## Convert image to OpenCV mat
	cv::Mat I= img_norm->GetOpenCVMat("64");
	cv::Mat p= I;
	eps *= 255*255;   // Because the intensity range of our images is [0, 255]

	//## Run guided filter
	cv::Mat dst = Caesar::guidedFilter(I, p, radius, eps);

	//## Fill filtered image
	TString imgName= Form("%s_GuidedFilter",m_name.c_str());
	Image* FilterImg= this->GetCloned(std::string(imgName),true,true);
	FilterImg->Reset();
	FilterImg->FillFromMat(dst);

	//## Clear allocated data
	if(img_norm) {
		delete img_norm;
		img_norm= 0;
	}

	return FilterImg;

}//close GetGuidedFilterImage()

Image* Image::GetSmoothedImage(int size_x,int size_y,double sigma_x,double sigma_y)
{

	//## Get OpenCV mat
	cv::Mat mat= this->GetOpenCVMat("64");
	
	//## Smooth matrix
	cv::Size smooth_size(size_x,size_y);
	cv::Mat smoothed_mat;
	cv::GaussianBlur(mat,smoothed_mat, smooth_size, sigma_x, sigma_y, cv::BORDER_DEFAULT);

	//## Fill smoothed image
	TString imgName= Form("%s_Smoothed",m_name.c_str());
	Image* SmoothedImg= this->GetCloned(std::string(imgName),true,true);
	SmoothedImg->Reset();
	SmoothedImg->FillFromMat(smoothed_mat);
	
	return SmoothedImg;

}//close GetSmoothedImage()

std::vector<Image*> Image::GetWaveletDecomposition(int nScales)
{
	std::vector<Image*> img_decomposition;
	img_decomposition= WTFilter::GetDecomposition(this,nScales);
	
	return img_decomposition;

}//close GetWaveletDecomposition()


Image* Image::GetSaliencyMap(int reso,double regFactor,int minRegionSize,double knnFactor,bool useRobust,double expFalloffPar,double distanceRegPar){

	//## Compute single-reso saliency map
	Image* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeSaliencyMap(this,reso,regFactor,minRegionSize,knnFactor,useRobust,expFalloffPar,distanceRegPar);
	if(!saliencyMap){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Saliency map estimation failed!");
		#endif
		return nullptr;
	}

	return saliencyMap;

}//close GetSaliencyMap()


Image* Image::GetMultiResoSaliencyMap(int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar,double salientMultiplicityThrFactor,bool addBkgMap,bool addNoiseMap,ImgBkgData* bkgData,double saliencyThrFactor,double imgThrFactor,bool useOptimalThr){

	//## Compute multi-reso saliency map
	Image* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeMultiResoSaliencyMap(this,resoMin,resoMax,resoStep,beta,minRegionSize,knnFactor,useRobustPars,expFalloffPar,distanceRegPar, salientMultiplicityThrFactor,addBkgMap,addNoiseMap,bkgData,saliencyThrFactor,imgThrFactor,useOptimalThr);
	if(!saliencyMap){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Multi-resolution saliency map estimation failed!");
		#endif
		return nullptr;
	}

	return saliencyMap;

}//close GetMultiResoSaliencyMap()


int Image::FindPeaks(std::vector<ImgPeak>& peakPoints,std::vector<int> kernelSizes, int peakShiftTolerance,bool skipBorders,int multiplicityThr)
{
	return MorphFilter::FindPeaks(peakPoints,this,kernelSizes,peakShiftTolerance,skipBorders,multiplicityThr);

}//close FindPeaks()

TGraph* Image::ComputePeakGraph(std::vector<int> kernelSizes,int peakShiftTolerance,bool skipBorders,int multiplicityThr)
{
	//Find peaks in image
	//std::vector<TVector2> peakPoints;
	std::vector<ImgPeak> peakPoints;
	if(this->FindPeaks(peakPoints,kernelSizes,peakShiftTolerance,skipBorders,multiplicityThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to find peaks in image!");
		#endif
		return nullptr;
	}

	//Fill peak graph
	TGraph* peakGraph= new TGraph(peakPoints.size());
	for(size_t i=0;i<peakPoints.size();i++){
		//peakGraph->SetPoint(i,peakPoints[i].X(),peakPoints[i].Y());
		peakGraph->SetPoint(i,peakPoints[i].x,peakPoints[i].y);
	}

	return peakGraph;

}//close ComputePeakGraph()

Image* Image::GetMorphDilatedImage(int kernSize,int niters,bool skipZeroPixels)
{
	return MorphFilter::ComputeMorphFilter(this,eMORPH_DILATION,kernSize,eMORPH_RECT,niters,skipZeroPixels);

}//close GetMorphDilatedImage()

Image* Image::GetMorphErodedImage(int kernSize,int niters,bool skipZeroPixels)
{
	return MorphFilter::ComputeMorphFilter(this,eMORPH_EROSION,kernSize,eMORPH_RECT,niters,skipZeroPixels);

}//close GetMorphErodedImage()

Image* Image::GetMorphClosingImage(int kernSize,int niters,bool skipZeroPixels)
{
	return MorphFilter::ComputeMorphFilter(this,eMORPH_CLOSING,kernSize,eMORPH_RECT,niters,skipZeroPixels);

}//close GetMorphClosingImage()

Image* Image::GetMorphOpeningImage(int kernSize,int niters,bool skipZeroPixels)
{
	return MorphFilter::ComputeMorphFilter(this,eMORPH_OPENING,kernSize,eMORPH_RECT,niters,skipZeroPixels);

}//close GetMorphOpeningImage()

Image* Image::GetMorphTopHatImage(int kernSize,int niters,bool skipZeroPixels)
{
	return MorphFilter::ComputeMorphFilter(this,eMORPH_TOPHAT,kernSize,eMORPH_RECT,niters,skipZeroPixels);

}//close GetMorphTopHatImage()

Image* Image::GetMorphGradientImage(int kernSize,int niters,bool skipZeroPixels)
{
	return MorphFilter::ComputeMorphFilter(this,eMORPH_GRADIENT,kernSize,eMORPH_RECT,niters,skipZeroPixels);

}//close GetMorphGradientImage()

Image* Image::GetMorphRecoImage(double baseline,int kernSize,double tol)
{
	return MorphFilter::ComputeMorphRecoFilter(this,baseline,kernSize,tol);

}//close GetMorphRecoImage()

Image* Image::GetHDomeImage(double baseline,int kernSize)
{
	return MorphFilter::ComputeHDomeFilter(this,baseline,kernSize);

}//close GetHDomeImage()

int Image::Add(Image* img,double c,bool computeStats)
{

	//Check input image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given input image!");
		#endif
		return -1;
	}
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	long int PixDataSize= img->GetPixelDataSize();
	if(m_Nx!=Nx || m_Ny!=Ny){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Image to be added has different size ("<<Nx<<","<<Ny<<") wrt to this image ("<<m_Nx<<","<<m_Ny<<")!");
		#endif
		return -1;
	}
	if(PixDataSize!=this->GetPixelDataSize()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Image to be added has different pixel vector size ("<<PixDataSize<<") wrt to this image ("<<m_pixels.size()<<")!");
		#endif
		return -1;
	}

	//Reset stats
	ResetImgStats(true,true);
	
	//Loop to sum vectors
	bool useNegativePixInStats= true;
	
	#ifdef OPENMP_ENABLED		
		Caesar::StatMoments<double> moments_t;	
		std::vector<Caesar::StatMoments<double>> parallel_moments;
	
		//#pragma omp declare reduction (merge : std::vector<Caesar::StatMoments<double>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
		#pragma omp parallel private(moments_t) //reduction(merge: parallel_moments)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Starting multithread image add (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");
			#endif

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			#pragma omp for 
			for(size_t i=0;i<m_pixels.size();i++){
				double w1= m_pixels[i];			
				double w2= img->GetPixelValue(i);	
				double w= w1 + c*w2;
				m_pixels[i]= 0;
				if(FillPixelMT(moments_t,i,w,useNegativePixInStats)<0) continue;
			}
			
			//parallel_moments.push_back(moments_t);
			parallel_moments[thread_id]= moments_t;
		}//close parallel section
		
		//Update moments from parallel estimates
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(m_StatMoments,parallel_moments)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
			#endif
			return -1;
		}

	#else 
		for(size_t i=0;i<m_pixels.size();i++){
			double w1= m_pixels[i];			
			double w2= img->GetPixelValue(i);	
			double w= w1 + c*w2;
			m_pixels[i]= 0;
			if(FillPixel(i,w,useNegativePixInStats)<0) continue;
		}
	#endif
	

	if(computeStats){
		bool computeRobustStats= true;
		//bool skipNegativePixels= false;
		bool useRange= false;
		bool forceRecomputing= false;
		//if(ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing)<0){
		if(ComputeStats(computeRobustStats,forceRecomputing,useRange)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute stats after adding the two images!");
			#endif
			return -1;
		}	
	}//close if computeStats

	return 0;

}//close Add()

int Image::Scale(double c)
{
	//Multiply pixel values by a factor c
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(size_t i=0;i<m_pixels.size();i++){
		double w_old= m_pixels[i];
		double w= w_old*c;
		m_pixels[i]= w;
	}

	//Force recomputation of stats if present, otherwise recompute only moments
	bool useRange= false;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	std::vector<float> maskedValues= {0};
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	if(this->HasStats()) status= this->ComputeStats(computeRobustStats,forceRecomputing,useRange);
	else status= this->ComputeMoments(useRange,minThr,maxThr,maskedValues);
		
	return status;

}//close Scale()

//================================
//==    THRESHOLDING METHODS
//================================
TH1D* Image::GetPixelHisto(int nbins,bool normalize){

	//Check if image has stats computed
	if(!HasStats()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No stats computed, returning nullptr!");
		#endif
		return nullptr;
	}

	double Smin= m_Stats->min;
	double Smax= m_Stats->max;
	double Srange= Smax-Smin;
	double tol= 0.0;
	double Smin_tol= Smin-tol*fabs(Srange);
	double Smax_tol= Smax+tol*fabs(Srange);

	TString histoName= Form("%s_histo",m_name.c_str());
	TH1D* histo= new TH1D(histoName,histoName,nbins,Smin_tol,Smax_tol);
		
	for(size_t i=0;i<m_pixels.size();i++){
		double w= m_pixels[i];
		if(w==0) continue;
		histo->Fill(w);
	}

	if(normalize) histo->Scale(1./histo->Integral());
	return histo;

}//close GetPixelHisto()


double Image::FindOtsuThreshold(TH1D* hist){

	//Check input histo
	if(!hist){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input pixel histo, returning inf!");
		#endif
		return TMath::Infinity();
	}

	//Compute Otsu threshold using pixel histo 
	int Nx= hist->GetNbinsX();
	double sum = 0;
	double Wxs[Nx];
	for (int t=0;t<Nx;t++) {
		double w= hist->GetBinContent(t+1);
		double x= hist->GetBinCenter(t+1);
		Wxs[t]= w*x;
		sum+= w*x;
	}

	double sumB = 0;
	double wB = 0;
	double wF = 0;

	double varMax = 0;
	double threshold = 0;
	
	for (int t=0;t<Nx;t++) {
		double binX= hist->GetBinCenter(t+1);
		double binContent= hist->GetBinContent(t+1);	
		double Wx= Wxs[t];
  	wB += binContent;               // Weight Background
   	if (wB == 0) continue;

   	wF = 1. - wB;                 // Weight Foreground
   	if (wF == 0) break;

		sumB+= Wx;
		double mB = sumB/wB;            // Mean Background
   	double mF = (sum - sumB)/wF;    // Mean Foreground

   	// Calculate Between Class Variance
   	double varBetween = wB*wF * (mB - mF) * (mB - mF);
		
   	// Check if new maximum found
   	if (varBetween > varMax) {
   		varMax = varBetween;
      threshold = binX;
   	}
	}//end loop bins

	if(hist) {
		hist->Delete();
		hist= 0;
	}

	return threshold;

}//close FindOtsuThreshold()


double Image::FindOtsuThreshold(int nbins){
	
	//## Get histo and normalize
	TH1D* hist= this->GetPixelHisto(nbins,true);
	if(!hist) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute pixel histo, return thr=0!");
		#endif
		return 0;
	}

	//## Compute Otsu threshold
	return FindOtsuThreshold(hist);

}//close FindOtsuThreshold()

double Image::FindMedianThreshold(double thrFactor){

	//Check if image has stats and get median
	if(!this->HasStats()){	
		#ifdef LOGGING_ENABLED
			INFO_LOG("Image has no stats computed, computing them...");
		#endif
		if(ComputeStats(true,false,false)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute image stats!");
			#endif
			return TMath::Infinity();
		}
	}
	double median= m_Stats->median;
	double medianThr= thrFactor*median;

	return medianThr;

}//close FindMedianThreshold()


double Image::FindOptimalGlobalThreshold(double thrFactor,int nbins,bool smooth){

	//Compute thresholds
	double medianThr= FindMedianThreshold(thrFactor);
	double otsuThr= FindOtsuThreshold(nbins);
	double valleyThr= FindValleyThreshold(nbins,smooth); 	
	if( TMath::IsNaN(medianThr) || fabs(medianThr)==TMath::Infinity() ||
			TMath::IsNaN(otsuThr) || fabs(otsuThr)==TMath::Infinity() ||
			TMath::IsNaN(valleyThr) || fabs(valleyThr)==TMath::Infinity() 
	)
	{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute one/more thresholds (inf values), returning inf!");
		#endif
		return TMath::Infinity();
	}

	//Compute optimal threshold
	double optimalThr= std::max(std::min(otsuThr,valleyThr),medianThr);

	return optimalThr;

}//close FindOptimalGlobalThreshold()


double Image::FindValleyThreshold(int nbins,bool smooth){

	//## Init dilate kernel size and histos
	const int nKernels= 3;
	int kernelSizes[]= {3,5,7};//{5,7,9}
	int maxStep= floor(kernelSizes[nKernels-1]/2.);

	//## Get pixel histo (invert to find peaks corresponding to valley in original histo)
	TH1D* histo= this->GetPixelHisto(nbins);
	if(!histo){		
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get pixel histo (hint: check if image stats were computed)!");
		#endif
		return 0;
	}

	if(smooth) histo->Smooth(1);
	histo->Scale(-1);
	double sMin= histo->GetMinimum();
	const double SMALL_NUMBER= 1.e-6;	
	for(int i=0;i<histo->GetNbinsX();i++){
		double w= histo->GetBinContent(i+1);
		double wnew= w-sMin+SMALL_NUMBER;
		histo->SetBinContent(i+1,wnew);
	}//end loop bins

	//## Get dilated histos
	TH1D* dilatedHisto= 0;
	std::vector<TH1D*> dilatedHistoList;

	for(int k=0;k<nKernels;k++){//loop over kernels
		TString histoName= Form("hdilate_%d",k+1);
		dilatedHisto= (TH1D*)histo->Clone(histoName);
		dilatedHisto->Reset();
		dilatedHistoList.push_back(dilatedHisto);
		
		int step= floor(kernelSizes[k]/2.);		

		for(int i=0;i<histo->GetNbinsX();i++){//loop over bins	
			double wmax= -1.e+99;
			for(int j=-step;j<step;j++){//Loop over kernel range
				int s= i+j;
				int binId= s+1;
				if(histo->IsBinUnderflow(binId) || histo->IsBinOverflow(binId)) continue;
				double w= histo->GetBinContent(binId);
				if(w>wmax) wmax= w;
			}//end loop kernel

			dilatedHistoList[k]->SetBinContent(i+1,wmax);
		}//end loop bins
	}//end loop kernels

	//## Find valleys
	TGraph* peaks= new TGraph;
	int npeaks= 0;
	double valleyThr= 0;
	for(int i=0;i<histo->GetNbinsX();i++){
		int binId= i+1;
		if(binId<maxStep+1) continue;
		double x= histo->GetXaxis()->GetBinCenter(i+1);
		double w= histo->GetBinContent(i+1);
		bool isPeak= true;
		for(int k=0;k<nKernels;k++) {	
			double wdilate= dilatedHistoList[k]->GetBinContent(binId); 
			if(wdilate!=w){
				isPeak= false;	
				break;
			}
		}//end loop kernels

		if(isPeak){
			peaks->SetPoint(npeaks,x,w);
			if(npeaks==0){
				valleyThr= x;
			}
			npeaks++;			
		}
	}//end loop bins

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<npeaks<<" valleys detected!");
	#endif

	//## Clear stuff
	if(histo){
		histo->Delete();
		histo= 0;
	}
	if(peaks) peaks->Delete();
	for(int k=0;k<nKernels;k++) {
		if(dilatedHistoList[k]) {
			dilatedHistoList[k]->Delete();
			dilatedHistoList[k]= 0;
		}
	}//end loop kernels
	dilatedHistoList.clear();

	return valleyThr;

}//close FindValleyThreshold()


double Image::FindCumulativeSumThr(double threshold,bool useRange,double minThr,double maxThr)
{
	//Check if image has pixels
	if(m_pixels.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Image has no pixels stored, returning 0!");
		#endif
		return 0;
	}

	//Copy pixel list
	//NB: This is required because of sorting
	std::vector<float> pixels;
	for(size_t i=0;i<m_pixels.size();i++){
		float w= m_pixels[i];
		//if( w==0 || (skipNegativePixels && w<0) ) continue;
		if( w==0 || (useRange && (w<=minThr || w>=maxThr)) ) continue;
		pixels.push_back(w);
	}

	//Find pixel value at which cumulative threshold is below desired threshold
	bool sorted= false;
	double thr= CodeUtils::FindCumulativeSumFractionThr(pixels,threshold,sorted);

	return thr;

}//close FindCumulativeSumThr()


Image* Image::GetBinarizedImage(double threshold,double fgValue,bool isLowerThreshold){

	TString imgName= Form("%s_Binarized",m_name.c_str());	
	Image* BinarizedImg= this->GetCloned(std::string(imgName),true,true);
	BinarizedImg->Reset();

	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(size_t i=0;i<m_pixels.size();i++){
		double w= m_pixels[i];
		if(w==0) continue;
		if(w>=threshold && !isLowerThreshold) {
			BinarizedImg->SetPixelValue(i,fgValue);
		}
		else if(w<threshold && isLowerThreshold) {
			BinarizedImg->SetPixelValue(i,fgValue);
		}
	}//end loop pixels
	
	//Force recomputation of stats if present, otherwise recompute only moments
	bool useRange= false;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	std::vector<float> maskedValues= {0};//CHECK
	if(BinarizedImg->ComputeMoments(useRange,minThr,maxThr,maskedValues)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to re-compute moments of binarized image!");
		#endif
		return nullptr;
	}	
	
	return BinarizedImg;

}//close GetBinarizedImage()


int Image::ApplyThreshold(double thr_min,double thr_max,double maskedValue)
{
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(size_t i=0;i<m_pixels.size();i++){
		double w= m_pixels[i];
		if(w==0) continue;
		if(w<thr_min || w>thr_max) {
			this->SetPixelValue(i,maskedValue);
		}
	}//end loop pixels
	
	//Force recomputation of stats if present, otherwise recompute only moments
	bool useRange= false;
	double minThr= -std::numeric_limits<double>::infinity();
	double maxThr= std::numeric_limits<double>::infinity();
	std::vector<float> maskedValues= {0};
	bool computeRobustStats= true;
	bool forceRecomputing= true;
	if(this->HasStats()){
		if(ComputeStats(computeRobustStats,forceRecomputing,useRange)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to re-compute stats of thresholded image!");
			#endif
			return -1;
		}
	}
	else{
		if(this->ComputeMoments(useRange,minThr,maxThr,maskedValues)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to re-compute moments of thresholded image!");
			#endif
			return -1;
		}	
	}

	return 0;

}//close ApplyThreshold()


//================================================================
//===    CONVERSION METHODS
//================================================================
TH2D* Image::GetHisto2D(std::string histoname)
{
	//Check if image is not empty
	if(m_pixels.empty() || m_Nx<=0 || m_Ny<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Image is empty or has no data/size stored, returning nullptr!");
		#endif
		return nullptr;
	}

	//Create 2D histo
	float xmin= m_Xmin - 0.5;
	float xmax= (m_Xmin+m_Nx-1) + 0.5;
	float ymin= m_Ymin - 0.5;
	float ymax= (m_Ymin+m_Ny-1) + 0.5;

	std::string hname= m_name;
	if(histoname!="") hname= histoname;
	TH2D* histo= new TH2D(hname.c_str(),hname.c_str(),m_Nx,xmin,xmax,m_Ny,ymin,ymax);
	histo->Sumw2();	
	
	//Fill histo
	for(long int j=0;j<m_Ny;j++){
		for(long int i=0;i<m_Nx;i++){
			long int gBin= this->GetBin(i,j);	
			double x= this->GetX(i); 
			double y= this->GetY(j); 
			double w= m_pixels[gBin];
			//histo->SetBinContent(i+1,j+1,w);
			histo->Fill(x,y,w);
		}
	}
	
	return histo;

}//close GetHisto2D()



TH2D* Image::GetWCSHisto2D(std::string histoname,WCS* wcs,bool useImageCoord)
{
	//Check if image is not empty
	if(m_pixels.empty() || m_Nx<=0 || m_Ny<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Image is empty or has no data/size stored, returning nullptr!");
		#endif
		return nullptr;
	}

	//Check WCS
	if(!wcs){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No WCS given!");
		#endif
		return nullptr;
	}

	//Compute bin x & y min values and convert to WCS
	std::vector<double> x_min;
	std::vector<double> y_min;
	for(long int i=0;i<m_Nx;i++){
		double x= i;
		double y= 0;

		//If using phys coords add xmin/ymin
		if(!useImageCoord){
			x+= m_Xmin;
			y+= m_Ymin;
		}	
		
		//Take pix low edge
		x-= 0.5;
		y-= 0.5;

		//Convert to WCS coords
		double x_wcs, y_wcs;
		if(AstroUtils::PixelToWCSCoords(x_wcs,y_wcs,wcs,x,y)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to convert x axis coords to WCS");
			#endif
			return nullptr;
		} 	
		x_min.push_back(x_wcs);
		
	}//end loop x axis
	
	for(long int i=0;i<m_Ny;i++){
		double x= 0;
		double y= i;

		//If using phys coords add xmin/ymin
		if(!useImageCoord){
			x+= m_Xmin;
			y+= m_Ymin;
		}	
		
		//Take pix low edge
		x-= 0.5;
		y-= 0.5;

		//Convert to WCS coords
		double x_wcs, y_wcs;
		if(AstroUtils::PixelToWCSCoords(x_wcs,y_wcs,wcs,x,y)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to convert y axis coords to WCS");
			#endif
			return nullptr;
		} 	
		y_min.push_back(y_wcs);
		
	}//end loop y axis

	//Sort axis ascending
	std::sort(x_min.begin(),x_min.end());
	std::sort(y_min.begin(),y_min.end());
	
	//Create histo 2D with variable bin width
	std::string hname= m_name;
	if(histoname!="") hname= histoname;
	int nBinsX= (int)(x_min.size()-1);
	int nBinsY= (int)(y_min.size()-1);
	TH2D* histo= new TH2D(hname.c_str(),hname.c_str(),nBinsX,x_min.data(),nBinsY,y_min.data());
	histo->Sumw2();	
	
	//Fill histo
	if(useImageCoord){
		for(long int j=0;j<m_Ny;j++){
			double y= j;
			for(long int i=0;i<m_Nx;i++){
				double x= i;
				long int gBin= this->GetBin(i,j);
				double w= m_pixels[gBin];
				double x_wcs, y_wcs;
				if(AstroUtils::PixelToWCSCoords(x_wcs,y_wcs,wcs,x,y)<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to convert bin ("<<x<<","<<y<<") coords to WCS!");
					#endif
					return nullptr;
				}
				histo->Fill(x,y,w);
			}//end loop x
		}//end loop y
	}//close if
	else{
		for(long int j=0;j<m_Ny;j++){
			double y= this->GetY(j);
			for(long int i=0;i<m_Nx;i++){
				double x= this->GetX(i); 
				long int gBin= this->GetBin(i,j);
				double w= m_pixels[gBin];
				double x_wcs, y_wcs;
				if(AstroUtils::PixelToWCSCoords(x_wcs,y_wcs,wcs,x,y)<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to convert bin ("<<x<<","<<y<<") coords to WCS!");
					#endif
					return nullptr;
				}
				histo->Fill(x,y,w);
			}//end loop x
		}//end loop y
	}//close else

	return histo;

}//close GetWCSHisto2D()

TMatrixD* Image::GetMatrix(){

	//Allocate matrix
	TMatrixD* M= new TMatrixD(m_Ny,m_Nx);
	M->Zero();

	//Fill matrix
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for collapse(2)
	#endif
	for(long int i=0;i<m_Nx;i++){//rows
		for(long int j=0;j<m_Ny;j++){//columns
			double w= this->GetPixelValue(i,j);
			long int rowId= j;
			long int colId= i;
			(*M)(rowId,colId)= w;
		}//end loop cols
	}//end loop rows

	return M;

}//close Image::GetMatrix()


std::string Image::GetPixelNumpyArrayStr()
{
	std::stringstream ss;
	ss<<"[";

	long int nRows = m_Ny;
  long int nCols = m_Nx;

	for(long int i=0;i<nRows;++i) {		
		long int rowId= i;
		long int iy= m_Ny-1-rowId;
  	
		ss<<"[";
    for (long int j=0;j<nCols;++j){
			int colId= j;
			int ix= colId;
			double w= this->GetPixelValue(ix,iy);
    	if(j==m_Ny-1) ss<<w;
			else ss<<w<<",";
    }//end loop cols
		if(i==m_Nx-1) ss<<"]";
		else ss<<"],";
  }//end loop rows

	/*
	for(long int i=0;i<m_Nx;i++){//rows
		ss<<"[";
		for(long int j=0;j<m_Ny;j++){//columns
			double w= this->GetPixelValue(i,j);
			long int rowId= j;
			long int colId= i;
			if(j==m_Ny-1) ss<<w;
			else ss<<w<<",";
		}//end loop cols
		if(i==m_Nx-1) ss<<"]";
		else ss<<"],";
	}//end loop rows
	*/
	ss<<"]";

	return ss.str();

}//close GetPixelNumpyArrayStr()

cv::Mat Image::GetOpenCVMat(std::string encoding){

	long int Nx= this->GetNx();
	long int Ny= this->GetNy();

	//## Fill OpenCV mat
	cv::Mat mat;
	if(encoding=="64") mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	else if(encoding=="32") mat= cv::Mat::zeros(Ny,Nx,CV_32FC1);
	else if(encoding=="32I") mat= cv::Mat::zeros(Ny,Nx,CV_32SC1);
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid encoding selected, using default 64bit encoding");
		#endif
		mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	}

	//The fast way
	long int nRows = mat.rows;
  long int nCols = mat.cols;
	
	if(encoding=="64"){
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(long int i=0;i<nRows;++i) {
			long int rowId= i;
			long int iy= Ny-1-rowId;
  		double* p = mat.ptr<double>(i);
    	for (long int j=0;j<nCols;++j){
				int colId= j;
				int ix= colId;
				double w= this->GetPixelValue(ix,iy);
    		p[j] = w;
    	}//end loop cols
  	}//end loop rows
	}//close if
	else if(encoding=="32"){
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(long int i=0;i<nRows;++i) {
			long int rowId= i;
			long int iy= Ny-1-rowId;
  		float* p = mat.ptr<float>(i);
    	for (long int j=0;j<nCols;++j){
				int colId= j;
				int ix= colId;
				float w= this->GetPixelValue(ix,iy);
    		p[j] = w;
    	}//end loop cols
  	}//end loop rows
	}//close else if
	else if(encoding=="32I"){
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(long int i=0;i<nRows;++i) {
			long int rowId= i;
			long int iy= Ny-1-rowId;
  		int* p = mat.ptr<int>(i);
    	for (long int j=0;j<nCols;++j){
				int colId= j;
				int ix= colId;
				int w= static_cast<int>(this->GetPixelValue(ix,iy));
    		p[j] = w;
    	}//end loop cols
  	}//end loop rows
	}//close else if
	else{
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(long int i=0;i<nRows;++i) {
			long int rowId= i;
			long int iy= Ny-1-rowId;
  		double* p = mat.ptr<double>(i);
    	for (long int j=0;j<nCols;++j){
				int colId= j;
				int ix= colId;
				double w= this->GetPixelValue(ix,iy);
    		p[j] = w;
    	}//end loop cols
  	}//end loop rows
	}//close else

	return mat;

}//close ImgToMat()

//================================================================
//===    DRAW  METHODS
//================================================================
int Image::Draw(int palette,bool drawFull,bool useCurrentCanvas,std::string units)
{
	//Set palette
	Caesar::GraphicsUtils::SetPalette(palette);

	//Delete existing object with name=htemp
	TObject* obj= gROOT->FindObject("htemp");
	if(obj){
		obj->Delete();
	}
	
	//Get temp histogram
	TH2D* htemp= GetHisto2D("htemp");
	if(!htemp){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get histo from this image!");
		#endif
		return -1;
	}
	
	//Set canvas
	int canvas_width= 720;
	int canvas_height= 700;
	TString canvasName= Form("%s_Plot",m_name.c_str());	
	TCanvas* canvas= 0;
	if(useCurrentCanvas && gPad) {
		canvas= gPad->GetCanvas();
		canvas->SetName(canvasName);
		canvas->SetTitle(canvasName);
	}
	else{
		canvas= new TCanvas(canvasName,canvasName,canvas_width,canvas_height);
	}

	if(!canvas){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to retrieve or set canvas!");
		#endif
		return -1;
	}

	//Draw full image (without borders)
	canvas->cd();
	htemp->SetStats(0);
	
	if(drawFull){
		canvas->ToggleEventStatus();
  	canvas->SetRightMargin(0.0);
  	canvas->SetLeftMargin(0.0);
  	canvas->SetTopMargin(0.0);
  	canvas->SetBottomMargin(0.0);
		htemp->Draw("COLA");
	}
	else{
		gStyle->SetPadTopMargin(0.1);
  	gStyle->SetPadBottomMargin(0.1);
  	gStyle->SetPadLeftMargin(0.15);
  	gStyle->SetPadRightMargin(0.15);

		gPad->SetTopMargin(0.1);
		gPad->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.15);
  	gPad->SetRightMargin(0.15);

		//Set palette axis title	
		htemp->GetZaxis()->SetTitle(units.c_str());
		htemp->GetZaxis()->SetTitleSize(0.05);
		htemp->GetZaxis()->SetTitleOffset(0.9);
		htemp->Draw();
		gPad->Update();
		htemp->Draw("COLZ");
		gPad->Update();
  }

	return 0;

}//close Draw()


int Image::Draw(std::vector<Source*>const& sources,int palette,bool drawFull,bool useCurrentCanvas,std::string units)
{
	
	//Draw image first
	Draw(palette,drawFull,useCurrentCanvas,units);

	//Retrieve canvas and draw sources	
	TCanvas* canvas= gPad->GetCanvas();
	if(!canvas){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to get access to current canvas!");
		#endif
		return -1;
	}

	//Draw sources in current canvas
	for(unsigned int k=0;k<sources.size();k++){	
		int type= sources[k]->Type;
		int lineColor= kBlack;
		if(type==eCompact)
			lineColor= kBlack;	
		else if(type==ePointLike)
			lineColor= kRed;
		else if(type==eExtended)
			lineColor= kGreen+1;	
		sources[k]->Draw(false,false,true,lineColor);
	}//end loop sources

	return 0;

}//close Draw()


int Image::Plot(std::vector<Source*>const& sources,bool useCurrentCanvas,bool drawFull,int paletteStyle,bool drawColorPalette,bool putWCSAxis,int coordSystem,std::string units){

	//Set palette
	Caesar::GraphicsUtils::SetPalette(paletteStyle);

	//Get temp histogram
	TH2D* htemp= GetHisto2D("htemp");
	if(!htemp){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get histo from this image!");
		#endif
		return -1;
	}
	
	//Set canvas
	TString canvasName= Form("%s_Plot",this->GetName().c_str());	
	TCanvas* canvas= 0;
	if(useCurrentCanvas && gPad) {
		canvas= gPad->GetCanvas();
		canvas->SetName(canvasName);
		canvas->SetTitle(canvasName);
	}
	else{
		canvas= new TCanvas(canvasName,canvasName,720,700);
	}

	if(!canvas){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to retrieve or set canvas!");
		#endif
		return -1;
	}

	//Draw full image (without borders)
	canvas->cd();
	htemp->SetStats(0);
	
	if(drawFull){
		canvas->ToggleEventStatus();
  	canvas->SetRightMargin(0.0);
  	canvas->SetLeftMargin(0.0);
  	canvas->SetTopMargin(0.0);
  	canvas->SetBottomMargin(0.0);
		htemp->Draw("COLA");
	}
	else{
		gStyle->SetPadTopMargin(0.1);
  	gStyle->SetPadBottomMargin(0.1);
  	gStyle->SetPadLeftMargin(0.15);
  	//gStyle->SetPadRightMargin(0.19);
		gStyle->SetPadRightMargin(0.15);

		
		gPad->SetTopMargin(0.1);
		gPad->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.15);
		//gPad->SetRightMargin(0.19);
  	gPad->SetRightMargin(0.15);
  		
		if(drawColorPalette) {
			//Set palette axis title	
			htemp->GetZaxis()->SetTitle(units.c_str());
			htemp->GetZaxis()->SetTitleSize(0.05);
			htemp->GetZaxis()->SetTitleOffset(0.9);
			htemp->Draw();
			gPad->Update();
			if(putWCSAxis) htemp->Draw("COLAZ");
			else htemp->Draw("COLZ");
			gPad->Update();
		}//close if
		else {
			if(putWCSAxis) htemp->Draw("COLA");
			else htemp->Draw("COL");
		}
	}


	//## Draw WCS axis?
	if(putWCSAxis){
		gPad->Update();
		TGaxis* xaxis_wcs= new TGaxis;
		TGaxis* yaxis_wcs= new TGaxis;
		int status= GraphicsUtils::SetWCSAxis(this,*xaxis_wcs,*yaxis_wcs,coordSystem);
		if(status>=0){
			TExec* ex = new TExec("ex","GraphicsUtils::PadUpdater()");
   		htemp->GetListOfFunctions()->Add(ex);

			xaxis_wcs->Draw("same");
			yaxis_wcs->Draw("same");
		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to set gAxis!");
			#endif
		}
	}//close if

	//## Draw sources
	for(unsigned int k=0;k<sources.size();k++){	
		int type= sources[k]->Type;
		int lineColor= kBlack;
		if(type==eCompact)
			lineColor= kBlack;	
		else if(type==ePointLike)
			lineColor= kRed;
		else if(type==eExtended)
			lineColor= kGreen+1;	
		sources[k]->Draw(false,false,true,lineColor);
	}//end loop sources
	
	return 0;

}//close Plot()


int Image::AddSimPointSources(long int nGenSources,float Smin,float Smax,float Sslope,int marginX,int marginY,std::string genFileName)
{
	//Check args
	if(nGenSources<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("nGenSources must be >0!");
		#endif
		return -1;
	}
	if(Smax<Smin){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Smax must be >= Smin!");
		#endif
		return -1;
	}
	if(Smax<0 || Smin<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Smax & Smin must be >0!");
		#endif
		return -1;
	}
	
	if (marginX<0 || marginX>=m_Nx/2.){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid marginX specified (<0 or larger than image half size)!");
		#endif
		return -1;
	}
	if (marginY<0 || marginY>=m_Ny/2.){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid marginY specified (<0 or larger than image half size)!");
		#endif
		return -1;
	}

	//Check metadata
	if(!m_HasMetaData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No metadata available in this image, cannot retrieve beam information needed for source generation!");
		#endif
		return -1;
	}

	//Check WCS
	WCS* wcs= m_MetaData->GetWCS(eJ2000);
	if(!wcs){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get WCS from metadata!");
		#endif
		return -1;
	}

	//Set generation pars
	bool randomizeFlux= (Smin!=Smax);
	double bmaj= m_MetaData->Bmaj*3600;//in arcsec
	double bmin= m_MetaData->Bmin*3600;//in arcsec
	double pa= m_MetaData->Bpa;//defined from north
	double pixSizeX= m_MetaData->dX*3600;//in arcsec 
	double pixSizeY= m_MetaData->dY*3600;//in arcsec 
	double pixSize= std::min(fabs(pixSizeX),fabs(pixSizeY));
	double bmaj_pix= bmaj/pixSize;//in pix
	double bmin_pix= bmin/pixSize;//in pix
	double sigmaX= bmaj_pix/GausSigma2FWHM;//in pix
	double sigmaY= bmin_pix/GausSigma2FWHM;//in pix
	double theta= pa + 90.;//# NB: BPA is the positional angle of the major axis measuring from North (up) counter clockwise, while theta is measured wrt to x axis
	double theta_rad= theta*TMath::DegToRad();//rad
	double lgS_min= log10(Smin);
	double lgS_max= log10(Smax);
	
	//Generate source positions
	
	long int index= 0;
	delete gRandom;
	gRandom= new TRandom3(0);
	std::vector<long int> takenBinIds;

	struct GenSourceInfo
	{
		std::string sname;
		double S;
		double x;
		double y;
		long int ix;
		long int iy;
		double ra;
		double dec;
		GenSourceInfo(std::string name,double _S,double _x,double _y,long int _ix,long int _iy,double _ra,double _dec)
			: sname(name),S(_S),x(_x),y(_y),ix(_ix),iy(_iy),ra(_ra),dec(_dec)
		{}
		
	};
	std::vector<GenSourceInfo> genSourceList;
	
	while(index < nGenSources)
	{				
		// Generate random coordinates
		double x0= m_Xmin + gRandom->Uniform(marginX,m_Nx-marginX-1);
		double y0= m_Ymin + gRandom->Uniform(marginY,m_Ny-marginY-1);
		long int gBin= this->FindBin(x0,y0);

		// Find if pixel was already taken by a previous source
		auto it= std::find(takenBinIds.begin(),takenBinIds.end(),gBin);
		if(!takenBinIds.empty() && it!=takenBinIds.end()){
			continue;
		}

		//Skip if flux in generation bin is =0 or NAN
		double w= this->GetPixelValue(gBin);
		if(w==0 || TMath::IsNaN(w) || fabs(w)==TMath::Infinity()){
			continue;
		}

		//Get bin WCS coords
		long int ix= this->GetBinX(gBin);
		long int iy= this->GetBinY(gBin);
		double ra= 0;
		double dec= 0;
		AstroUtils::PixelToWCSCoords(ra,dec,wcs,x0,y0);

		// Generate source Speak expo in log
		double S= Smin;
		if(randomizeFlux){
			double lgS_rand= gRandom->Exp(1./Sslope);
			double lgS= lgS_min + lgS_rand;
			if(lgS>lgS_max) continue;
			S= pow(10,lgS);
		}
			

		// Add source to list
		index++;
		std::string sname= Form("S%ld",index);
		genSourceList.push_back(GenSourceInfo(sname,S,x0,y0,ix,iy,ra,dec));

	}//end while loop

	//Add sources to map
	long int sourceGenBoxHalfSize= static_cast<long int>(std::round(10*max(sigmaX,sigmaY)));
	#ifdef LOGGING_ENABLED
		INFO_LOG("Generating sources: sigma("<<sigmaX<<","<<sigmaY<<"), theta="<<theta<<", genBoxSize="<<2*sourceGenBoxHalfSize<<" ...");
	#endif

	for(size_t k=0;k<genSourceList.size();k++)
	{
		double xs= genSourceList[k].x;
		double ys= genSourceList[k].y;
		long int ix= genSourceList[k].ix;
		long int iy= genSourceList[k].iy;
		double xs_bin= this->GetX(ix);
		double ys_bin= this->GetY(iy);
		double x0= xs-xs_bin;
		double y0= ys-ys_bin;

		double Speak= genSourceList[k].S;

		for(long int i=ix-sourceGenBoxHalfSize;i<=ix+sourceGenBoxHalfSize;i++)
		{
			double x_bin= this->GetX(i);
			double x= x_bin - xs_bin;

			for(long int j=iy-sourceGenBoxHalfSize;j<=iy+sourceGenBoxHalfSize;j++){
				long int id= this->GetBin(i,j);
				double y_bin= this->GetY(j);
				double y= y_bin - ys_bin;
				if(id<0) continue;
				double S_old= this->GetPixelValue(id);
				double w= MathUtils::EvalGaus2D(x,y,Speak,x0,y0,sigmaX,sigmaY,theta_rad);
				if(TMath::IsNaN(w) || fabs(w)==TMath::Infinity()) {
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Given source gen flux value (w="<<w<<") for bin "<<id<<" is NaN or inf, skipping!");
					#endif
					continue;
				}
				double S_new= S_old + w;
				this->SetPixelValue(id,S_new);

			}//end loop y
		}//end loop x

	}//end loop sources


	//Write true source info to ascii file
	FILE* fout= fopen(genFileName.c_str(),"w");
	for(size_t i=0;i<genSourceList.size();i++)
	{
		std::string sname= genSourceList[i].sname;
		double x= genSourceList[i].x;
		double y= genSourceList[i].y;
		double ra= genSourceList[i].ra;
		double dec= genSourceList[i].dec;
		double S= genSourceList[i].S;
		fprintf(fout,"%s %f %f %f %f %f\n",sname.c_str(),x,y,ra,dec,S);
	}
	fclose(fout);

	return 0;

}//close AddSimPointSources()

int Image::AddSimPointSourceDensity(double sourceDensity,float Smin,float Smax,float Sslope,int marginX,int marginY,std::string genFileName)
{
	//Check margin
	if (marginX<0 || marginX>=m_Nx/2.){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid marginX specified (<0 or larger than image half size)!");
		#endif
		return -1;
	}
	if (marginY<0 || marginY>=m_Ny/2.){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid marginY specified (<0 or larger than image half size)!");
		#endif
		return -1;
	}

	//Check metadata
	if(!m_HasMetaData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No metadata available in this image, cannot retrieve beam information needed for source generation!");
		#endif
		return -1;
	}
	double pixSizeX= fabs(m_MetaData->dX);//in deg
	double pixSizeY= fabs(m_MetaData->dY);//in deg 
	
	//Compute number of sources to be generated given map area in pixels
	long int nx_gen= m_Nx - 2*marginX;
	double nx_gen_deg= nx_gen*pixSizeX;
	long int ny_gen= m_Ny - 2*marginY;
	double ny_gen_deg= ny_gen*pixSizeY;
	  
	double genArea= (nx_gen_deg*ny_gen_deg);
	long int nSources= static_cast<long int>(std::round(sourceDensity*genArea));
	#ifdef LOGGING_ENABLED
		INFO_LOG("Generating #"<<nSources<<" (density="<<sourceDensity<<") ...");
	#endif

	return AddSimPointSources(nSources,Smin,Smax,Sslope,marginX,marginY,genFileName);

}//close AddSimPointSourceDensity()
	



}//close namespace
