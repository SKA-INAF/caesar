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

#include <Logger.h>

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


//OpenCV headers
#include <opencv/cv.h>
#include <opencv/highgui.h>
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
		ERROR_LOG(ss.str());
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
		ERROR_LOG(ss.str());
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
		ERROR_LOG(ss.str());
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
		ERROR_LOG(ss.str());
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
		ERROR_LOG(ss.str());
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
		ERROR_LOG(ss.str());
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

Image::Image(const Image &img)
{
	// Copy constructor
	DEBUG_LOG("Copy constuctor called...");
  ((Image&)img).Copy(*this);
}//close copy constructor



void Image::Copy(TObject& obj) const
{
	DEBUG_LOG("Copying parent TNamed...");
	TNamed::Copy((Image&)obj);

	DEBUG_LOG("Copying main vars...");
	((Image&)obj).m_name= m_name;
	((Image&)obj).m_Nx= m_Nx;
	((Image&)obj).m_Ny= m_Ny;
	((Image&)obj).m_Xmin= m_Xmin;
	((Image&)obj).m_Ymin= m_Ymin;
	
	DEBUG_LOG("Copying pixel collection...");
	((Image&)obj).m_pixels= m_pixels;
		
	DEBUG_LOG("Copying stats vars...");
	((Image&)obj).m_HasMetaData= m_HasMetaData;
	((Image&)obj).m_HasStats= m_HasStats;
	((Image&)obj).m_StatMoments= m_StatMoments;
	/*
	((Image&)obj).m_Npix= m_Npix;
	((Image&)obj).m_M1= m_M1;
	((Image&)obj).m_M2= m_M2;
	((Image&)obj).m_M3= m_M3;
	((Image&)obj).m_M4= m_M4;
	((Image&)obj).m_PixelMin= m_PixelMin;
	((Image&)obj).m_PixelMax= m_PixelMax;
	*/

	DEBUG_LOG("Copying meta data...");
	if(m_MetaData){
		((Image&)obj).m_MetaData= new ImgMetaData;
		*((Image&)obj).m_MetaData = *m_MetaData;
	}
	
	DEBUG_LOG("Copying stats data...");
	((Image&)obj).m_StatMoments= m_StatMoments;
	if(m_Stats){
		((Image&)obj).m_Stats= new ImgStats;
		*((Image&)obj).m_Stats = *m_Stats;
	}
	
}//close Copy()


Image::~Image(){

	if(m_Stats){
		DEBUG_LOG("Deleting stats...");	
		delete m_Stats;
		m_Stats= 0;
	}

	if(m_MetaData){
		DEBUG_LOG("Deleting meta-data...");
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
		ERROR_LOG("Invalid (empty or negative) sizes given, nothing will be done!");
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
		ERROR_LOG("C++ exception occurred while allocating memory (N="<<N<<") for pixel vector!");
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
	DEBUG_LOG("Read FITS image "<<filename<<" in "<<dt<<" ms");	

	//Check error status
	if(status<0){
		ERROR_LOG("Failed to fill image from FITS file!");
		return -1;
	}

	return 0;

}//close ReadFITS()


int Image::WriteFITS(std::string outfilename){
		
	//## Write fits to file
	if(FITSWriter::WriteFITS(this,outfilename)<0){
		ERROR_LOG("Failed to write image to FITS file!");
		return -1;
	}
	
	return 0;

}//close WriteFITS()


int Image::ReadFile(std::string filename,bool invert){

	//## Detect file extension
	std::string extension= filename.substr(filename.find_last_of(".") + 1);
	if(extension!= "png" && extension!="jpg" && extension!="bmp" && extension!="gif" ) {
		ERROR_LOG("Unknown file extension detected: ext="<<extension<<" (valid ones are png/jpg/bmp/gif)!");
		return -1;
	}
	
	//## Load image from file and set a matrix
	cv::Mat mat = cv::imread(filename.c_str(), CV_LOAD_IMAGE_COLOR);

	//## Convert to gray scale
	cv::Mat mat_gray;
  cvtColor( mat, mat_gray, CV_RGB2GRAY );
	
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

Image* Image::GetTile(long int ix_min,long int ix_max,long int iy_min,long int iy_max,std::string imgname){

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max)<0){
		ERROR_LOG("Failed to extract pixel rectangular data, returning nullptr!");
		return nullptr;
	}

	long int TileSizeX= ix_max - ix_min + 1;
	long int TileSizeY= iy_max - iy_min + 1;
	float xlow= ix_min;
	float ylow= iy_min;
	std::string name= imgname;
	if(name=="") name= "tileimg";

	DEBUG_LOG("ix min/max="<<ix_min<<"/"<<ix_max<<", iy min/max="<<iy_min<<"/"<<iy_max<<" TileSizeX="<<TileSizeX<<", TileSizeY="<<TileSizeY<<", pixels size="<<tile_pixels.size());
	
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
		WARN_LOG("Invalid pixel bin ("<<gbin<<") given to be filled (hint: check if image size is not initialized or given bin exceed size)!");
		return -1;
	}
	if(TMath::IsNaN(w) || fabs(w)==TMath::Infinity()) {
		DEBUG_LOG("Given value (w="<<w<<") for bin "<<gbin<<" is NaN or inf, skipping!");
		return 0;
	}

	//Compute binx & biny
	long int binx = GetBinX(gbin);
 	long int biny = GetBinY(gbin);

	//Check if pixel value has been already filled
	double w_old= m_pixels[gbin];
	if( w_old!=0 ){
		WARN_LOG("Pixel "<<gbin<<" ("<<binx<<","<<biny<<") has been already filled (w="<<w_old<<"), skipping...");
		return -1;
	}
	
	return 0;

}//close CheckFillPixel()


int Image::FillPixel(long int gbin,double w,bool useNegativePixInStats){

	//Check bin & value
	if(CheckFillPixel(gbin,w)<0) {
		WARN_LOG("Invalid bin given (gbin="<<gbin<<") or or pixel already filled!");
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


int Image::Fill(double x,double y,double w,bool useNegativePixInStats){

	//Find global bin
	long int gbin= FindBin(x,y);
	if(gbin<0 || gbin>=(long int)(m_pixels.size())) {
		WARN_LOG("Underflow/overflow bin ("<<gbin<<") corresponding to (x,y)=("<<x<<","<<y<<") given to be filled (hint: check if image size is not initialized)!");
		return -1;
	}
	
	//Fill pixel
	return FillPixel(gbin,w,useNegativePixInStats);

}//close FillPixel()


#ifdef OPENMP_ENABLED
int Image::FillPixelMT(Caesar::StatMoments<double>& moments,long int gbin,double w,bool useNegativePixInStats){

	//Check bin & value
	if(CheckFillPixel(gbin,w)<0) {
		WARN_LOG("Invalid bin given (gbin="<<gbin<<") or nan/inf value given (w="<<w<<") or pixel already filled!");
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
int Image::FillPixelMT(Caesar::StatMoments<double>& moments,long int ix,long int iy,double w,bool useNegativePixInStats){

	//Compute global bin
	long int gbin= GetBin(ix,iy);
	
	//Fill pixel
	return FillPixelMT(moments,gbin,w,useNegativePixInStats);

}//close FillPixelMT()
#endif

#ifdef OPENMP_ENABLED
int Image::FillMT(Caesar::StatMoments<double>& moments,double x,double y,double w,bool useNegativePixInStats){

	//Find global bin
	long int gbin= FindBin(x,y);
	if(gbin<0) {
		WARN_LOG("Underflow/overflow bin requested to be filled!");
	}

	//Fill pixel
	return FillPixelMT(moments,gbin,w,useNegativePixInStats);

}//close FillPixelMT()
#endif


int Image::FillFromMat(cv::Mat& mat,bool useNegativePixInStats){

	//Get image size
	long int nRows = mat.rows;
  long int nCols = mat.cols;
	long int Ny= this->GetNy();
	long int Nx= this->GetNx();
	if(nRows<=0 || nCols<=0){
		ERROR_LOG("Matrix has zero rows and/or cols!");
		return -1;
	}
	if(nCols!=Nx || nRows!=Ny){
		ERROR_LOG("Invalid nrows/ncols (hint: should be equal to image size "<<Nx<<" x "<<Ny<<")!");
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
			ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
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
		ERROR_LOG("Matrix has zero rows and/or cols!");
		return -1;
	}
	if(nCols!=Nx || nRows!=Ny){
		ERROR_LOG("Invalid nrows/ncols (hint: should be equal to image size "<<Nx<<" x "<<Ny<<")!");
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
			ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
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

#ifdef VTK_ENABLED
int Image::FillFromVtkImage(vtkSmartPointer<vtkImageData> imageData,bool useNegativePixInStats)
{
	//Check image data
	if(!imageData){
		ERROR_LOG("Null ptr to input VTK image data given!");
		return -1;
	}

	//Check dims
	int* dims = imageData->GetDimensions();
	if(dims[2]!=1){
		ERROR_LOG("Invalid z dimension (must be =1) in VTK image data!");
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
			ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
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

int Image::ComputeMoments(bool skipNegativePixels){

	//## Recompute stat moments
	//## NB: If OMP is enabled this is done in parallel and moments are aggregated to return the correct cumulative estimate 
	bool maskNanInfValues= true;
	int status= Caesar::StatsUtils::ComputeStatsMoments(m_StatMoments,m_pixels,skipNegativePixels,maskNanInfValues);
	if(status<0){
		ERROR_LOG("Failed to compute stat moments!");
	}

	return status;

}//close ComputeMoments()



void Image::ComputeStatsParams(bool computeRobustStats,bool skipNegativePixels,bool useParallelVersion){

	//## Reset previous stats params (Reset only Stats not moments!)
	ResetImgStats(false);

	//-- DEBUG
	//DEBUG_LOG("Npix="<<m_Npix<<", M1="<<m_M1<<", m_M2="<<m_M2<<", m_M3="<<m_M3<<", m_M4="<<m_M4<<", pixel sizes="<<m_pixels.size());
	DEBUG_LOG("Npix="<<m_StatMoments.N<<", M1="<<m_StatMoments.M1<<", M2="<<m_StatMoments.M2<<", M3="<<m_StatMoments.M3<<", M4="<<m_StatMoments.M4<<", pixel sizes="<<m_pixels.size());
	//----
	
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
	INFO_LOG("Image standard stats computed in "<<dt_stats<<" ms");
  
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
		if( w==0 || (skipNegativePixels && w<0) ) continue;
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
	INFO_LOG("Image median pars computed in "<<dt_median<<" ms");
  
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
	INFO_LOG("Image biweight pars computed in "<<dt_biweights<<" ms");
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
	INFO_LOG("Image clipped stats pars computed in "<<dt_clipped<<" ms");

	auto end_robuststats = chrono::steady_clock::now();
	double dt_robuststats= chrono::duration <double, milli> (end_robuststats-start_robuststats).count();
	INFO_LOG("Image robust stats computed in "<<dt_robuststats<<" ms");	
	
}//close ComputeStatsParams()



int Image::ComputeStats(bool computeRobustStats,bool skipNegativePixels,bool forceRecomputing,bool useParallelVersion){

	
	//## Start timer
	auto start = chrono::steady_clock::now();

	
	//## Check if image has already stats computed
	if(!HasStats()){
		m_Stats= new Caesar::ImgStats;
	}
	else{		
		WARN_LOG("Image has stats already computed...");
	}

	//## If recomputing is not requested (i.e. some pixels has been reset by the user, just set the stats params!
	if(!forceRecomputing){
		ComputeStatsParams(computeRobustStats,skipNegativePixels,useParallelVersion);
		m_HasStats= true;
		
		auto stop = chrono::steady_clock::now();
		double elapsed_time= chrono::duration <double, milli> (stop-start).count();
		INFO_LOG("Image stats computed in "<<elapsed_time<<" ms");
		return 0;
	}

	
	//## Recompute the moments and stats params
	DEBUG_LOG("Recomputing image stats...");

	//--> Reset stats
	ResetImgStats(true);

	//--> Recompute moments
	ComputeMoments(skipNegativePixels);

	//--> Recompute stats params
	ComputeStatsParams(computeRobustStats,skipNegativePixels,useParallelVersion);

	m_HasStats= true;

	//Stop timer
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	DEBUG_LOG("Image stats recomputed in "<<dt<<" ms");
	
	return 0;

}//close ComputeStats()


int Image::GetTilePixels(std::vector<float>& pixels,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{
	//Check range given
	if(ix_min<0 || ix_max<0 || ix_max>=m_Nx || ix_min>=ix_max){
		ERROR_LOG("Invalid x range given!");
		return -1;
	}
	if(iy_min<0 || iy_max<0 || iy_max>=m_Ny || iy_min>=iy_max){
		ERROR_LOG("Invalid y range given!");
		return -1;
	}
	
	//Extract pixel sub vector
	pixels.clear();

	if(skipNegativePixels){

		for(long int j=iy_min;j<=iy_max;j++){
			//Get row start/end iterators
			long int gBin_min= GetBin(ix_min,j);		
			long int gBin_max= GetBin(ix_max,j);

			std::copy_if (
				m_pixels.begin() + gBin_min, 
				m_pixels.begin() + gBin_max + 1, 
				std::back_inserter(pixels), 
				[](float w){
					return (w>=0);
				}
			);
		}//end loop rows
	}//close if

	else{

		for(long int j=iy_min;j<=iy_max;j++){
			//Get row start/end iterators
			long int gBin_min= GetBin(ix_min,j);		
			long int gBin_max= GetBin(ix_max,j);		

			//Extract row and append to vector
			pixels.insert(pixels.end(), m_pixels.begin() + gBin_min, m_pixels.begin() + gBin_max + 1);

		}//end loop row
	}//close else
	
	return 0;

}//close GetTilePixels()


int Image::GetTileMeanStats(float& mean,float& stddev,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{
	//Init
	mean= 0;
	stddev= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
		return -1;
	}
	
	//Compute mean/rms
	Caesar::StatsUtils::ComputeMeanAndRMS(mean,stddev,tile_pixels);
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileMeanStats()
		
int Image::GetTileMedianStats(float& median,float& mad_rms,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{

	//Init
	median= 0;
	mad_rms= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
		return -1;
	}
	
	//Compute median/mad
	bool useParallelVersion= false;
	median= Caesar::StatsUtils::GetMedianFast(tile_pixels,useParallelVersion);
	mad_rms= Caesar::StatsUtils::GetMADFast(tile_pixels,median,useParallelVersion)*1.4826;
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileMedianStats()


int Image::GetTileClippedStats(ClippedStats<float>& clipped_stats,long int& npix,double clipSigma,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{
	//Init
	npix= 0;
	
	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
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
		

int Image::GetTileBiWeightStats(float& bwLocation,float& bwScale,long int& npix,double C,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{
	//Init
	bwLocation= 0;
	bwScale= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
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
ImgBkgData* Image::ComputeBkg(int estimator,bool computeLocalBkg,int boxSizeX,int boxSizeY, double gridStepSizeX,double gridStepSizeY,bool use2ndPass,bool skipOutliers,double seedThr,double mergeThr,int minPixels){
	
	//## Compute bkg data
	DEBUG_LOG("Using grid bkg method...");
	ImgBkgData* bkgData= BkgFinder::FindBkg(this,estimator,computeLocalBkg, boxSizeX, boxSizeY, gridStepSizeX,gridStepSizeY,use2ndPass,skipOutliers,seedThr,mergeThr,minPixels);	
	if(!bkgData){
		ERROR_LOG("Computation of local background failed for this image!");
		return 0;
	}
	return bkgData;
	
}//close ComputeBkg()


Image* Image::GetSignificanceMap(ImgBkgData* bkgData,bool useLocalBkg){

	//Check image
	if(!bkgData){
		ERROR_LOG("Null ptr to bkg data!");
		return 0;
	}

	//Integrity check for local bkg
	long int Nx= this->GetNx();
	long int Ny= this->GetNy();
	if(useLocalBkg){
		if(!bkgData->HasLocalBkg()){
			ERROR_LOG("Local bkg option requested but no local bkg data are available!");
			return 0;
		}
		if( Nx!=(bkgData->BkgMap)->GetNx() || Ny!=(bkgData->BkgMap)->GetNy() ||
				Nx!=(bkgData->NoiseMap)->GetNx() || Ny!=(bkgData->NoiseMap)->GetNy()
		){
			ERROR_LOG("Bkg/Noise maps have different size!");		
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
							DEBUG_LOG("Empty pixel or invalid bkg  ("<<i<<","<<j<<") skip it!");
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
			ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
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
						DEBUG_LOG("Empty pixel or invalid bkg  ("<<i<<","<<j<<") skip it!");
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
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	if(this->FindCompactSource(sources,this,0,thr,thr,minPixels,findNegativeExcess,mergeBelowSeed,findNestedSources)<0){
		ERROR_LOG("Compact source finder failed!");
		return -1;
	}
	return 0;

}//close FindCompactSource()


int Image::FindCompactSource(std::vector<Source*>& sources,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,bool findNestedSources,double nestedBlobThreshold,double minNestedMotherDist,double maxMatchingPixFraction,Image* curvMap)
{

	//Find sources
	int status= BlobFinder::FindBlobs(this,sources,floodImg,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed,curvMap);
	if(status<0){
		ERROR_LOG("Blob finder failed!");
		for(unsigned int k=0;k<sources.size();k++){
			if(sources[k]){
				delete sources[k];
				sources[k]= 0;
			}	
		}//end loop sources
		sources.clear();
		return -1;
	}

	//Find nested sources?
	if(findNestedSources && sources.size()>0){
		int status= FindNestedSource(sources,bkgData,minPixels,nestedBlobThreshold,minNestedMotherDist,maxMatchingPixFraction);
		if(status<0){
			WARN_LOG("Nested source search failed!");
		}
	}//close if

	return 0;

}//close FindCompactSource()


int Image::FindNestedSource(std::vector<Source*>& sources,ImgBkgData* bkgData,int minPixels,double nestedBlobThreshold,double minNestedMotherDist,double maxMatchingPixFraction){

	//Check if given mother source list is empty
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0){
		WARN_LOG("Empty source list given!");
		return 0;
	}

	//Find image mask of found sources
	Image* sourceMask= this->GetSourceMask(sources,false);
	if(!sourceMask){
		ERROR_LOG("Null ptr to computed source mask!");
		return -1;
	}

	//Find curvature map
	Image* curvMap= this->GetLaplacianImage(true);
	if(!curvMap){
		ERROR_LOG("Failed to compute curvature map!");
		delete sourceMask;
		sourceMask= 0;
		return -1;
	}

	//Compute curvature map stats
	if(curvMap->ComputeStats(true,false,false)<0){
		ERROR_LOG("Failed to compute curvature map stats!");
		delete sourceMask;
		sourceMask= 0;
		delete curvMap;
		curvMap= 0;
		return -1;
	}

	//Thresholding the curvature map
	double curvMapRMS= curvMap->GetPixelStats()->medianRMS;
	double curvMapThr= curvMapRMS*nestedBlobThreshold;
	Image* blobMask= curvMap->GetBinarizedImage(curvMapThr);
	if(!blobMask){
		ERROR_LOG("Failed to compute curvature blob mask!");
		delete sourceMask;
		sourceMask= 0;
		delete curvMap;
		curvMap= 0;
		return -1;
	}

	//Find blob+source mask
	Image* sourcePlusBlobMask= sourceMask->GetMask(blobMask,true);
	if(!sourcePlusBlobMask){
		ERROR_LOG("Failed to compute (source+blob) mask!");
		delete sourceMask;
		sourceMask= 0;
		delete curvMap;
		curvMap= 0;
		delete blobMask;
		blobMask= 0;
		return -1;
	}

	//Find nested blobs 
	std::vector<Source*> NestedSources;
	double fgValue= 1;
	int status= BlobFinder::FindBlobs(this,NestedSources,sourcePlusBlobMask,bkgData,fgValue,fgValue,minPixels,false,false);
	if(status<0){
		ERROR_LOG("Nested blob finder failed!");
		delete sourceMask;
		sourceMask= 0;
		delete curvMap;
		curvMap= 0;
		delete blobMask;
		blobMask= 0;
		delete sourcePlusBlobMask;
		sourcePlusBlobMask= 0;
		for(size_t k=0;k<NestedSources.size();k++){
			if(NestedSources[k]){
				delete NestedSources[k];
				NestedSources[k]= 0;
			}	
		}//end loop sources
		NestedSources.clear();
		return -1;
	}//close if

	//Add nested sources to mother source
	int nNestedSources= static_cast<int>(NestedSources.size());
	
	
	if(nNestedSources>=0){
		INFO_LOG("#"<<nNestedSources<<" blobs found in curvature map!");

		//## Init mother-nested association list
		std::vector<int> nestedSourcesToBeRemoved;
		std::vector< std::vector<int> > MotherNestedAssociationList;
		for(int i=0;i<nSources;i++) MotherNestedAssociationList.push_back( std::vector<int>() );

		//## Find matching between mother and nested sources
		for(int j=0;j<nNestedSources;j++){
			bool isMotherFound= false;
			DEBUG_LOG("Finding matching for nested source no. "<<j);
			
			for(int i=0;i<nSources;i++){
				int sourceId= sources[i]->Id;
				bool isInside= NestedSources[j]->IsInsideSource(sources[i]);
				
				if(isInside){
					DEBUG_LOG("Nested source no. "<<j<<" added to source id="<<sourceId<<" ...");
					MotherNestedAssociationList[i].push_back(j);
					//NestedSources[j]->ComputeStats();
					//NestedSources[j]->ComputeMorphologyParams();
					//sources[i]->AddNestedSource(NestedSources[j]);
					isMotherFound= true;
					break;
				}
				
			}//end loop mother sources

			//If nested is not associated to any mother source, mark for removal
			if(!isMotherFound){
				WARN_LOG("Cannot find mother source for nested source no. "<<j<<", will remove it from the list of nested sources...");
				nestedSourcesToBeRemoved.push_back(j);
			}			

		}//end loop nested sources	

		//## Select nested
		int nSelNestedSources= 0;
		for(size_t i=0;i<MotherNestedAssociationList.size();i++){
			long int NPix= sources[i]->GetNPixels(); 
			int nComponents= static_cast<int>(MotherNestedAssociationList[i].size());
			if(nComponents<=0) continue;

			//If only one component is present select it if:
			//  1) mother and nested distance is > thr (e.g. 
			//  2) mother and nested pix superposition is <thr (e.g. 50%)
			for(int j=0;j<nComponents;j++){
				int nestedIndex= MotherNestedAssociationList[i][j];

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
					INFO_LOG("areOffset? "<<areOffset<<" (dist_x="<<centroidDistX<<", dist_y="<<centroidDistY<<"), isNestedSmaller?"<<isNestedSmaller<<" (matchingPixFraction="<<matchingPixFraction<<", maxMatchingPixFraction="<<maxMatchingPixFraction<<")");					

					if( areOffset || isNestedSmaller){//Add nested to mother source
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

		INFO_LOG("#"<<nSelNestedSources<<" nested sources found and added to mother sources...");

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

	//Clear
	if(sourceMask) {
		delete sourceMask;
		sourceMask= 0;
	}
	if(curvMap) {
		delete curvMap;
		curvMap= 0;
	}
	if(blobMask) {
		delete blobMask;
		blobMask= 0;
	}
	if(sourcePlusBlobMask) {
		delete sourcePlusBlobMask;
		sourcePlusBlobMask= 0;
	}

	return 0;

}//close FindNestedSources()

int Image::FindExtendedSource_CV(std::vector<Source*>& sources,Image* initSegmImg,ImgBkgData* bkgData,int minPixels,bool findNegativeExcess,double dt,double h,double lambda1,double lambda2,double mu,double nu,double p,int niters){

	//## Compute segmented image
	Image* segmentedImg= ChanVeseSegmenter::FindSegmentation(this,initSegmImg,false,dt,h,lambda1,lambda2,mu,nu,p,niters);
	if(!segmentedImg){
		ERROR_LOG("Failed to compute ChanVese image segmentation!");
		return -1;
	}
	
	//## Finding blobs in masked image
	double fgValue= 1;	
	int status= this->FindCompactSource(sources,segmentedImg,bkgData,fgValue,fgValue,minPixels,false,false,false);
	if(status<0){
		ERROR_LOG("Finding sources in Chan-Vese segmented mask failed!");
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
			ERROR_LOG("Failed to compute image stats, clearing all detected sources and returning empty list!");
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
			sources[k]->SetType(Source::eExtended);

			//Check if source is a "negative excess"
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);
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
		ERROR_LOG("Null ptr to given image mask!");
		return 0;
	}
		
	//## Check mask bins
	long int Nx= mask->GetNx();
	long int Ny= mask->GetNy();
	if(Nx!=this->GetNx() || Ny!=this->GetNy()){
		ERROR_LOG("Mask binning is different than current image!");
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


int Image::MaskSources(std::vector<Source*>const& sources,float maskValue)
{
	//## Check source list
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0) {
		WARN_LOG("Source list is empty, nothing to be done...");
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
	bool skipNegativePixels= false;
	bool forceRecomputing= true;
	this->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);

	return 0;

}//close MaskSources()

Image* Image::GetSourceMask(std::vector<Source*>const& sources,bool isBinary,bool invert){

	//## Clone map
	long int Nx= this->GetNx();
	long int Ny= this->GetNy();
	bool copyMetaData= true;
	bool resetStats= true;
	TString imgName= Form("%s_SourceMask",m_name.c_str());	
	Image* maskedImage= this->GetCloned(std::string(imgName),copyMetaData,resetStats);
	
	//## Check source list
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0) {
		WARN_LOG("Source list is empty, returning same image!");
		return maskedImage;	
	}

	if(invert){
		if(isBinary){
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					maskedImage->SetPixelValue(id,0);
				}//end loop pixels
			}//end loop sources	

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
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					maskedImage->SetPixelValue(id,0);
				}//end loop pixels
			}//end loop sources		
		}//close else

		//Force re-computation of stats
		maskedImage->ComputeStats(true,false,true);

	}//close if invert
	else{
		//Reset map and loop over sources
		maskedImage->Reset();

		if(isBinary){	
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif		
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					//maskedImage->FillPixel(id,1);
					maskedImage->SetPixelValue(id,1);
				}//end loop pixels
			}//end loop sources		
		}//close if
		else{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					double w= this->GetPixelValue(id);
					//maskedImage->FillPixel(id,w);
					maskedImage->SetPixelValue(id,w);
				}//end loop pixels
			}//end loop sources		
		}//close else

		//Force re-computation of stats
		maskedImage->ComputeStats(true,false,true);

	}//close else

	return maskedImage;

}//close GetSourceMask()

Image* Image::GetSourceResidual(std::vector<Source*>const& sources,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr){

	//Clone input image	
	TString imgName= Form("%s_Residual",m_name.c_str());	
	Image* residualImg= this->GetCloned(std::string(imgName),true,true);
	if(!residualImg){
		ERROR_LOG("Failed to clone input image, returning nullptr!");
		return nullptr;
	}

	//Dilate source pixels
	int status= MorphFilter::DilateAroundSources(residualImg,sources,KernSize,dilateModel,dilateSourceType,skipToNested,bkgData,useLocalBkg,randomize,zThr);

	if(status<0){
		ERROR_LOG("Failed to dilate sources!");
		if(residualImg) {
			delete residualImg;
			residualImg= 0;
		}
		return 0;		
	}

	//Re-Compute stats
	residualImg->ComputeStats(true,false,true);

	return residualImg;

}//close GetSourceResidual()

Image* Image::GetNormalizedImage(std::string normScale,int normmin,int normmax,bool skipEmptyBins){

	//Check if image has data
	if(m_pixels.empty()){
		ERROR_LOG("Image has no data stored!");
		return nullptr;
	}

	//Get image content min/max
	double wmin= (this->m_StatMoments).minVal;
	double wmax= (this->m_StatMoments).maxVal;
	if(TMath::IsNaN(wmin) || fabs(wmin)==TMath::Infinity() || TMath::IsNaN(wmax) || fabs(wmax)==TMath::Infinity()){
		ERROR_LOG("Min ("<<wmin<<") or max ("<<wmax<<") image values are nan or inf, this should not occur (hint: check moment calculation!)");
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
		WARN_LOG("Invalid norm scale option selected ("<<normScale<<") no transform applied to original image!");
		return norm_img;
	}

	//Force recomputation of stats if present, otherwise recompute only moments
	bool skipNegativePixels= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	if(this->HasStats()) status= norm_img->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
	else status= norm_img->ComputeMoments(skipNegativePixels);
	if(status<0){
		WARN_LOG("Failed to re-compute moments/stats for normalized image!");
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
		ERROR_LOG("Failed to get normalized image!");
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

std::vector<Image*> Image::GetWaveletDecomposition(int nScales){

	DEBUG_LOG("Computing wavelet decomposition up to scale J="<<nScales<<" ...");
	std::vector<Image*> img_decomposition;
	img_decomposition= WTFilter::GetDecomposition(this,nScales);
	
	return img_decomposition;

}//close GetWaveletDecomposition()


Image* Image::GetSaliencyMap(int reso,double regFactor,int minRegionSize,double knnFactor,bool useRobust,double expFalloffPar,double distanceRegPar){

	//## Compute single-reso saliency map
	Image* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeSaliencyMap(this,reso,regFactor,minRegionSize,knnFactor,useRobust,expFalloffPar,distanceRegPar);
	if(!saliencyMap){
		ERROR_LOG("Saliency map estimation failed!");
		return nullptr;
	}

	return saliencyMap;

}//close GetSaliencyMap()


Image* Image::GetMultiResoSaliencyMap(int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar,double salientMultiplicityThrFactor,bool addBkgMap,bool addNoiseMap,ImgBkgData* bkgData,double saliencyThrFactor,double imgThrFactor){

	//## Compute multi-reso saliency map
	Image* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeMultiResoSaliencyMap(this,resoMin,resoMax,resoStep,beta,minRegionSize,knnFactor,useRobustPars,expFalloffPar,distanceRegPar, salientMultiplicityThrFactor,addBkgMap,addNoiseMap,bkgData,saliencyThrFactor,imgThrFactor);
	if(!saliencyMap){
		ERROR_LOG("Multi-resolution saliency map estimation failed!");
		return nullptr;
	}

	return saliencyMap;

}//close GetMultiResoSaliencyMap()


int Image::FindPeaks(std::vector<TVector2>& peakPoints,int peakShiftTolerance,bool skipBorders)
{
	return MorphFilter::FindPeaks(peakPoints,this,peakShiftTolerance,skipBorders);

}//close FindPeaks()

TGraph* Image::ComputePeakGraph(int peakShiftTolerance,bool skipBorders)
{
	//Find peaks in image
	std::vector<TVector2> peakPoints;
	if(this->FindPeaks(peakPoints,skipBorders)<0){
		ERROR_LOG("Failed to find peaks in image!");
		return nullptr;
	}

	//Fill peak graph
	TGraph* peakGraph= new TGraph(peakPoints.size());
	for(size_t i=0;i<peakPoints.size();i++){
		peakGraph->SetPoint(i,peakPoints[i].X(),peakPoints[i].Y());
	}

	return peakGraph;

}//close ComputePeakGraph()


int Image::Add(Image* img,double c,bool computeStats)
{

	//Check input image
	if(!img){
		ERROR_LOG("Null ptr to given input image!");
		return -1;
	}
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	long int PixDataSize= img->GetPixelDataSize();
	if(m_Nx!=Nx || m_Ny!=Ny){
		ERROR_LOG("Image to be added has different size ("<<Nx<<","<<Ny<<") wrt to this image ("<<m_Nx<<","<<m_Ny<<")!");
		return -1;
	}
	if(PixDataSize!=this->GetPixelDataSize()){
		ERROR_LOG("Image to be added has different pixel vector size ("<<PixDataSize<<") wrt to this image ("<<m_pixels.size()<<")!");
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
			DEBUG_LOG("Starting multithread image add (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");

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
			ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
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
		bool skipNegativePixels= false;
		bool forceRecomputing= false;
		if(ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing)<0){
			WARN_LOG("Failed to compute stats after adding the two images!");
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
	bool skipNegativePixels= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	if(this->HasStats()) status= this->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
	else status= this->ComputeMoments(skipNegativePixels);
		
	return status;

}//close Scale()

//================================
//==    THRESHOLDING METHODS
//================================
TH1D* Image::GetPixelHisto(int nbins,bool normalize){

	//Check if image has stats computed
	if(!HasStats()){
		WARN_LOG("No stats computed, returning nullptr!");
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
		ERROR_LOG("Null ptr to input pixel histo, returning inf!");
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
		ERROR_LOG("Failed to compute pixel histo, return thr=0!");
		return 0;
	}

	//## Compute Otsu threshold
	return FindOtsuThreshold(hist);

}//close FindOtsuThreshold()

double Image::FindMedianThreshold(double thrFactor){

	//Check if image has stats and get median
	if(!this->HasStats()){
		INFO_LOG("Image has no stats computed, computing them...");
		if(ComputeStats(true,false,false)<0){
			ERROR_LOG("Failed to compute image stats!");
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
		ERROR_LOG("Failed to compute one/more thresholds (inf values), returning inf!");
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

	DEBUG_LOG("#"<<npeaks<<" valleys detected!");
	
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
	bool skipNegativePixels= false;
	if(BinarizedImg->ComputeMoments(skipNegativePixels)<0){
		ERROR_LOG("Failed to re-compute moments of binarized image!");
		return nullptr;
	}	
	
	return BinarizedImg;

}//close GetBinarizedImage()

//================================================================
//===    CONVERSION METHODS
//================================================================
TH2D* Image::GetHisto2D(std::string histoname){

	//Check if image is not empty
	if(m_pixels.empty() || m_Nx<=0 || m_Ny<=0){
		WARN_LOG("Image is empty or has no data/size stored, returning nullptr!");
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

}//close Img::GetMatrix()

cv::Mat Image::GetOpenCVMat(std::string encoding){

	long int Nx= this->GetNx();
	long int Ny= this->GetNy();

	//## Fill OpenCV mat
	cv::Mat mat;
	if(encoding=="64") mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	else if(encoding=="32") mat= cv::Mat::zeros(Ny,Nx,CV_32FC1);
	else{
		WARN_LOG("Invalid encoding selected, using default 64bit encoding");
		mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	}

	//The fast way
	long int nRows = mat.rows;
  long int nCols = mat.cols;
	
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
    }
  }

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
		ERROR_LOG("Failed to get histo from this image!");
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
		ERROR_LOG("Failed to retrieve or set canvas!");
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
		WARN_LOG("Failed to get access to current canvas!");
		return -1;
	}

	//Draw sources in current canvas
	for(unsigned int k=0;k<sources.size();k++){	
		int type= sources[k]->Type;
		int lineColor= kBlack;
		if(type==Source::eCompact)
			lineColor= kBlack;	
		else if(type==Source::ePointLike)
			lineColor= kRed;
		else if(type==Source::eExtended)
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
		ERROR_LOG("Failed to get histo from this image!");
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
		ERROR_LOG("Failed to retrieve or set canvas!");
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
			WARN_LOG("Failed to set gAxis!");
		}
	}//close if

	//## Draw sources
	for(unsigned int k=0;k<sources.size();k++){	
		int type= sources[k]->Type;
		int lineColor= kBlack;
		if(type==Source::eCompact)
			lineColor= kBlack;	
		else if(type==Source::ePointLike)
			lineColor= kRed;
		else if(type==Source::eExtended)
			lineColor= kGreen+1;	
		sources[k]->Draw(false,false,true,lineColor);
	}//end loop sources
	
	return 0;

}//close Plot()

/*
void Image::SetDrawRange(double zmin,double zmax){

	if(zmin>=zmax) return;

	for(long int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double w= this->GetBinContent(i+1,j+1);
			if(w<zmin) this->SetBinContent(i+1,j+1,zmin);
			if(w>zmax) this->SetBinContent(i+1,j+1,zmax);
		}//end loop bins y
	}//end loop bins x

}//close Img::SetDrawRange()
*/


}//close namespace
