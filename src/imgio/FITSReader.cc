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
* @file FITSReader.cc
* @class FITSReader
* @brief FITSReader
*
* FITS Image Reader class
* @author S. Riggi
* @date 20/01/2015
*/

#include <FITSReader.h>
#include <Image.h>
#include <SysUtils.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>


#include <TObject.h>
#include <TMath.h>
#include <TVectorD.h>


#ifdef OPENMP_ENABLED
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif

//CFITSIO headers
#include <fitsio.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <chrono>
using namespace std;

ClassImp(Caesar::FITSHeader)
ClassImp(Caesar::FITSFileInfo)
ClassImp(Caesar::FITSReader)

namespace Caesar {


FITSReader::FITSReader() {

}//close costructor


FITSReader::~FITSReader() {

}//close destructor

void FITSReader::HandleError(int& status,fitsfile* fp)
{
	char errdescr[FLEN_STATUS+1];
	fits_get_errstatus(status, errdescr);
	#ifdef LOGGING_ENABLED
 		ERROR_LOG("Error while processing FITS file (reason="<<errdescr<<")");
	#endif
  status = 0;
  if (fp) fits_close_file(fp, &status);

}//close HandleError()


int FITSReader::Read(Caesar::Image& img,Caesar::FITSFileInfo& fits_info,std::string filename,int hdu_id,int ix_min,int ix_max,int iy_min,int iy_max,bool checkFile)
{
	//## Check file
	if(checkFile && !SysUtils::CheckFile(filename,fits_info.info,true,".fits")){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read file "<<filename<<"!");
		#endif
		return -1;
	}

	//## Check if tile reading option was given
	//## If so append coord filters to input file (see CFITSIO for details)
	bool readFull= (ix_min==-1 || ix_max==-1 || iy_min==-1 || iy_max==-1);

	std::string filename_wfilter= filename;
	if(!readFull){
		std::stringstream ss;
		ss<<filename<<"["<<ix_min+1<<":"<<ix_max+1<<","<<iy_min+1<<":"<<iy_max+1<<"]";
		filename_wfilter= ss.str();
	}

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading image "<<filename_wfilter<<"...");
	#endif

	//## Init vars
	fitsfile* fp= 0;
  int status = 0;

	//## Open file with filter
  fits_open_file(&fp, filename_wfilter.c_str(), READONLY, &status);
  if (status) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open FITS file "<<filename_wfilter<<"!");
		#endif
		HandleError(status,fp);
		return -1;
	}

	//## Read HDU number
  int hdu_number= 0;
  fits_get_hdu_num(fp, &hdu_number);
	
	//## Move to the desired HDU in input file
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading hdu no. "<<hdu_id<<"/"<<hdu_number<<" in FITS file "<<filename<<"...");
	#endif
	int hdutype;
  if ( fits_movabs_hdu(fp, hdu_id, &hdutype, &status) ){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to access HDU "<<filename<<" (hint: check if HDU exists)!");
		#endif
		HandleError(status,fp);
		return -1;
	}

  //## Read HDU type and check if it is an IMAGE_HDU
	//## NB: Table not supported
  int hdu_type;
  fits_get_hdu_type(fp, &hdu_type, &status);
  if (status) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get HDU type!");
		#endif
   	HandleError(status,fp);
		return -1;
	}

	if(hdu_type!=IMAGE_HDU){
		#ifdef LOGGING_ENABLED
  	 	WARN_LOG("HDU not supported (hduType="<<hdu_type<<"), only IMAGE_HDU supported.");
		#endif
   	if (fp) fits_close_file(fp, &status);
		return -1;
	}

  //## Read HDU header records
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading HDU header ...");
	#endif
  if(ReadHeader(fits_info,fp)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read header of FITS file "<<filename<<"!");
		#endif
		if (fp) fits_close_file(fp, &status);
		return -1;
	}

	
  //## Read HDU's image data
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading HDU image data ...");
	#endif

	#ifdef OPENMP_ENABLED
		//Multi-thread reading
		if(ReadImageMT(img,fits_info,fp,filename,ix_min,ix_max,iy_min,iy_max)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read image data of FITS file "<<filename<<"!");
			#endif
			if (fp) fits_close_file(fp, &status);
			return -1;
		}
	#else
		//Single-thread reading
		if(ReadImage(img,fits_info,fp,filename,ix_min,ix_max,iy_min,iy_max)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read image data of FITS file "<<filename<<"!");
			#endif
			if (fp) fits_close_file(fp, &status);
			return -1;
		}
	#endif
	
	
	//## Close FITS file
	if (fp) {
		fits_close_file(fp, &status);
		fp= 0;
	}

	return 0;

}//close Read()



int FITSReader::ReadImage(Image& img,Caesar::FITSFileInfo& fits_info,fitsfile* fp,std::string filename,int ix_min,int ix_max,int iy_min,int iy_max)
{
	auto start = chrono::steady_clock::now();

	//## Check if ptr is null
	if(!fp) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to FITS file given!");
		#endif
		return -1;
	}

	//Read the image size from header
	long int Nx= (fits_info.header).Nx; 
	long int Ny= (fits_info.header).Ny;

	//Check if tile reading mode is given
	bool readTile= (ix_min!=-1 && ix_max!=-1 && iy_min!=-1 && iy_max!=-1);
	long int ImgSizeX= Nx;
	long int ImgSizeY= Ny;
	long int xlow= 0;
	long int ylow= 0;
	std::string filename_wfilter= filename;
	if(readTile){
		ImgSizeX= ix_max-ix_min+1;
		ImgSizeY= iy_max-iy_min+1;
		xlow= ix_min;
		ylow= iy_min;
		std::stringstream ss;
		ss<<filename<<"["<<ix_min+1<<":"<<ix_max+1<<","<<iy_min+1<<":"<<iy_max+1<<"]";
		filename_wfilter= ss.str();
	}

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading image data (Nx x Ny)=("<<ImgSizeX<<","<<ImgSizeY<<")");
	#endif

	//Set image size (Allocate memory for image pixels)	
	if(img.SetSize(ImgSizeX,ImgSizeY,xlow,ylow)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to set image size to ("<<ImgSizeX<<","<<ImgSizeY<<") (hint: check if enough memory is available)!");
		#endif
		return -1;
	}

	//Set image meta-data	
	Caesar::ImgMetaData* metadata= new Caesar::ImgMetaData;
	metadata->SetFITSCards(fits_info);
	img.SetMetaData(metadata);

  
	//For multithreading we cannot share the same file pointer across different thread (see CFITSIO multithreading ref in doc)
	//Open file multiple times?
	fitsfile* fp_safe= 0;
	int errflag= 0;	
	int readdata_errflag= 0;

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Starting serial-version image data reading ...");
	#endif

	//For serial single thread use the already opened file ptr
	fp_safe= fp;

	if(ReadAndFillImageDataFast(img,ImgSizeX,ImgSizeY,fp_safe,readdata_errflag)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and fill the image data!");
		#endif
		errflag++;
	}
	
	//Check for errors
	int close_file_status= 0;
	if(errflag>0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("One/more errors occurred while reading and filling the image!");
		#endif
		if (fp) fits_close_file(fp, &close_file_status);		
		return -1;
	}	

	//Stop timer
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image read in "<<dt<<" ms");
	#endif

	return 0;

}//close ReadImage()


#ifdef OPENMP_ENABLED
int FITSReader::ReadImageMT(Image& img,Caesar::FITSFileInfo& fits_info,fitsfile* fp,std::string filename,int ix_min,int ix_max,int iy_min,int iy_max)
{
	//## Check if ptr is null
	if(!fp) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to FITS file given!");
		#endif
		return -1;
	}

	//Read the image size from header
	long int Nx= (fits_info.header).Nx; 
	long int Ny= (fits_info.header).Ny;

	//Check if tile reading mode is given
	bool readTile= (ix_min!=-1 && ix_max!=-1 && iy_min!=-1 && iy_max!=-1);
	long int ImgSizeX= Nx;
	long int ImgSizeY= Ny;
	long int xlow= 0;
	long int ylow= 0;
	std::string filename_wfilter= filename;
	if(readTile){
		ImgSizeX= ix_max-ix_min+1;
		ImgSizeY= iy_max-iy_min+1;
		xlow= ix_min;
		ylow= iy_min;
		std::stringstream ss;
		ss<<filename<<"["<<ix_min+1<<":"<<ix_max+1<<","<<iy_min+1<<":"<<iy_max+1<<"]";
		filename_wfilter= ss.str();
	}

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading image data (Nx x Ny)=("<<ImgSizeX<<","<<ImgSizeY<<")");
	#endif

	//Set image size (Allocate memory for image pixels)	
	if(img.SetSize(ImgSizeX,ImgSizeY,xlow,ylow)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to set image size to ("<<ImgSizeX<<","<<ImgSizeY<<") (hint: check if enough memory is available)!");
		#endif
		return -1;
	}

	//Set image meta-data	
	Caesar::ImgMetaData* metadata= new Caesar::ImgMetaData;
	metadata->SetFITSCards(fits_info);
	img.SetMetaData(metadata);

  //Close existing file pointer
	/*
	int close_status= 0;
	if(fp) {
		fits_close_file(fp, &close_status);
		fp= 0;
	}
	*/
	
	//For multithreading we cannot share the same file pointer across different thread (see CFITSIO multithreading ref in doc)
	//Open file multiple times?
	fitsfile* fp_safe= 0;
	int errflag= 0;	
	int readdata_errflag= 0;
	
	Caesar::StatMoments<double> moments_t;	
	//std::vector<Caesar::StatMoments<double>> parallel_moments;
	
	//#pragma omp declare reduction (merge : std::vector<Caesar::StatMoments<double>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
	#pragma omp parallel private(fp_safe,moments_t) reduction(+: errflag)	// reduction(merge: parallel_moments)
	{
		int thread_id= omp_get_thread_num();
		int nthreads= SysUtils::GetOMPThreads();
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Starting multithread image data reading (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");
		#endif

		//Open file (each thread should open its own)
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Opening FITS file ptr (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");
		#endif

		int status= 0;
		#pragma omp critical 
		{
			fits_open_file(&fp_safe, filename_wfilter.c_str(), READONLY, &status);
		}			
  	if (status) {
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to open FITS file "<<filename_wfilter<<" in thread "<<omp_get_thread_num()<<"!");
			#endif
			errflag++;
		}

		/*
		#pragma omp single
   	{
     	parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   	}
		*/

		//Read image data in parallel (each thread with its own file pointer)
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Reading image data (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");
		#endif

		if(ReadAndFillImageDataFastMT(img,ImgSizeX,ImgSizeY,fp_safe,readdata_errflag)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read and fill the image data!");
			#endif
			errflag++;
		}

		//Close open file (executed by each thread)
		if(fp_safe) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Closing FITS file (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");
			#endif
			fits_close_file(fp_safe, &status);
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Closed FITS file (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");
			#endif
		}

	}//close parallel section		

	//Check for errors
	int close_file_status= 0;
	if(errflag>0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("One/more errors occurred while reading and filling the image!");
		#endif
		if (fp) {
			fits_close_file(fp, &close_file_status);		
			fp= 0;
		}
		return -1;
	}	

	//Compute moments	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing image moments after pixel read...");
	#endif
	bool skipNegativePixels= false;
	if(img.ComputeMoments()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute image moments!");
		#endif
		return -1;
	}

	return 0;

}//close ReadImageMT()
#endif


int FITSReader::ReadAndFillImageDataFast(Image& img,long int Nx,long int Ny,fitsfile* fp,int& err_flag)
{
	auto start = chrono::steady_clock::now();

	//Check pointer
	if(!fp) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to  FITS file given!");
		#endif
		return -1;
	}

	//Get number of image dimensions	
	int status= 0;
	int naxis= 0;
	fits_get_img_dim(fp, &naxis, &status);
	if(status){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get FITS image size!");
		#endif
		return -1;
	}

	//Initialize fpixel 
	long int fpixel[naxis];
	for(int i=0;i<naxis;i++) fpixel[i]= 1;
	
	//Fill image data
	for(long int j=0;j<Ny;j++){
		
		//Init buffer
		//long int fpixel[2] = {1,j+1};//start from column 0th, row j-th 
		fpixel[0]= 1;
		fpixel[1]= j+1;
		long int bufsize= Nx;
		long int gBin= 0 + j*Nx;
		
		//Read all columns	
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Reading row "<<j<<" (nelements="<<bufsize<<")...");
		#endif

		float nullval= 0; //don't check for null values in the image
		int anynull;
		int read_status= 0;
		fits_read_pix(fp,TFLOAT,fpixel,bufsize,(void*)&nullval,(void*)( (img.m_pixels).data()+gBin ),&anynull, &read_status);
		if(read_status){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read FITS image data row "<<j<<" (bufsize="<<bufsize<<"), skip...");
			#endif
			HandleError(read_status,fp);
			err_flag++;
			continue;
		}

	}//end loop row

	if(err_flag>0) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read fits and fill image pixels!");
		#endif
		return -1;
	}

	//Compute moments	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing image moments after pixel read...");
	#endif
	bool skipNegativePixels= false;
	if(img.ComputeMoments()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute image moments!");
		#endif
		return -1;
	}

	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image data read in "<<dt<<" ms");
	#endif

	return 0;

}//close ReadAndFillImageDataFast()


int FITSReader::ReadAndFillImageData(Image& img,long int Nx,long int Ny,fitsfile* fp,int& err_flag)
{
	//Check pointer
	if(!fp) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to  FITS file given!");
		#endif
		return -1;
	}

	//Get number of image dimensions	
	int status= 0;
	int naxis= 0;
	fits_get_img_dim(fp, &naxis, &status);
	if(status){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get FITS image size!");
		#endif
		return -1;
	}

	//Initialize fpixel 
	long int fpixel[naxis];
	for(int i=0;i<naxis;i++) fpixel[i]= 1;
	
	//Fill image data
	for(long int j=0;j<Ny;j++){
		
		//Init buffer
		//long int fpixel[2] = {1,j+1};//start from column 0th, row j-th 
		fpixel[0]= 1;
		fpixel[1]= j+1;
		long int bufsize= Nx;
		float* buffer= new float[bufsize];
		
		//Read all columns
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Reading row "<<j<<" (nelements="<<bufsize<<")...");
		#endif

		float nullval= 0; //don't check for null values in the image
		int anynull;
		int read_status= 0;
		fits_read_pix(fp,TFLOAT,fpixel,bufsize,(void*)&nullval,(void*)buffer,&anynull, &read_status);
		if(read_status){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read FITS image data row "<<j<<" (bufsize="<<bufsize<<"), skip...");
			#endif
			HandleError(read_status,fp);
			err_flag++;
			continue;
		}
		
		//Loop over buffer and fill image
		for(long int i=0;i<bufsize;i++){
			double w= buffer[i];
			if(img.FillPixel(i,j,w)<0) continue;
		}//end loop columns (bufsize)
		
		//Clear allocated buffer
		if(buffer) delete [] buffer;

	}//end loop row

	if(err_flag>0) return -1;
	return 0;

}//close ReadAndFillImageData()

#ifdef OPENMP_ENABLED
int FITSReader::ReadAndFillImageDataMT(Image& img,long int Nx,long int Ny,fitsfile* fp,int& err_flag,std::vector<Caesar::StatMoments<double>>& parallel_moments)
{
	//Check pointer
	if(!fp) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to FITS file given!");
		#endif	
		return -1;
	}
	
	//Get number of image dimensions	
	int status= 0;
	int naxis= 0;
	fits_get_img_dim(fp, &naxis, &status);
	if(status){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get FITS image size!");
		#endif
		return -1;
	}

	//Initialize fpixel 
	long int fpixel[naxis];
	for(int i=0;i<naxis;i++) fpixel[i]= 1;

	//Define parallel stat moments list	
	Caesar::StatMoments<double> moments_t;	
	int nthreads= SysUtils::GetOMPThreads();
	int thread_id= omp_get_thread_num();
	
	#pragma omp for private(moments_t) // reduction(+: err_flag)	
	for(long int j=0;j<Ny;j++){
		
		//Init buffer
		//long int fpixel[2] = {1,j+1};//start from column 0th, row j-th 
		fpixel[0]= 1;
		fpixel[1]= j+1;

		long int bufsize= Nx;
		float* buffer= new float[bufsize];
			
		//Read all columns
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Reading row "<<j<<" (nelements="<<bufsize<<")...");
		#endif
			
		float nullval= 0; //don't check for null values in the image
		int anynull;
		int read_status= 0;
		fits_read_pix(fp,TFLOAT,fpixel,bufsize,(void*)&nullval,(void*)buffer,&anynull, &read_status);
		if(read_status){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read FITS image data row "<<j<<", skip...");
			#endif
	
			#pragma omp critical
			{
				err_flag++;
			}
			continue;
		}
		
		//Loop over buffer and fill image
		for(long int i=0;i<bufsize;i++){
			double w= buffer[i];
			if(img.FillPixelMT(moments_t,i,j,w)<0) continue;
		}//end loop columns (bufsize)

		//Fill moment list 
		parallel_moments[thread_id]= moments_t;

		//Clear allocated buffer
		if(buffer) delete [] buffer;

	}//end loop row

	if(err_flag>0) return -1;
	return 0;

}//close ReadAndFillImageDataMT()
#endif


#ifdef OPENMP_ENABLED
int FITSReader::ReadAndFillImageDataFastMT(Image& img,long int Nx,long int Ny,fitsfile* fp,int& err_flag)
{
	//Check pointer
	if(!fp) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to FITS file given!");
		#endif
		return -1;
	}

	//Get number of image dimensions	
	int status= 0;
	int naxis= 0;
	fits_get_img_dim(fp, &naxis, &status);
	if(status){
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Failed to get FITS image size!");
		#endif
		return -1;
	}

	//Initialize fpixel 
	long int fpixel[naxis];
	for(int i=0;i<naxis;i++) fpixel[i]= 1;
	
	//Define parallel stat moments list	
	int nthreads= SysUtils::GetOMPThreads();
	int thread_id= omp_get_thread_num();
	
	//#pragma omp for private(moments_t) 
	#pragma omp for
	for(long int j=0;j<Ny;j++){
		
		//Init buffer
		//long int fpixel[2] = {1,j+1};//start from column 0th, row j-th 
		fpixel[0]= 1;
		fpixel[1]= j+1;
		long int bufsize= Nx;
		long int gBin= 0 + j*Nx;
		
		//Read all columns
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Reading row "<<j<<" (nelements="<<bufsize<<")...");
		#endif

		float nullval= 0; //don't check for null values in the image
		int anynull;
		int read_status= 0;
		fits_read_pix(fp,TFLOAT,fpixel,bufsize,(void*)&nullval,(void*)( (img.m_pixels).data()+gBin ),&anynull, &read_status);
		if(read_status){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read FITS image data row "<<j<<", skip...");
			#endif

			#pragma omp critical
			{
				err_flag++;
			}
			continue;
		}

	}//end loop row

	if(err_flag>0) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read fits and fill image pixels!");
		#endif
		return -1;
	}	


	return 0;

}//close ReadAndFillImageDataFastMT()
#endif


int FITSReader::ReadHeader(FITSFileInfo& fits_info,fitsfile* fp)
{
	//Check input FITS file ptr
	if(!fp) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given FITS file!");
		#endif
		return -1;
	}

	//Read FITS number of keywords
	int status= 0;
	int nkeys, morekeys;
  
  fits_get_hdrspace(fp, &nkeys, &morekeys, &status);
  if (status) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get the number of existing keywords in FITS file!");
		#endif
		HandleError(status,fp);
		return -1;
	}

	//Fill keywords in map
	std::map<std::string,std::string> records; 
	char keyname[FLEN_KEYWORD+1];
  char keyvalue[FLEN_VALUE+1];
  char comment[FLEN_COMMENT+1];

  for (int i=1;i<=nkeys;i++) {
  	fits_read_keyn(fp, i, keyname, keyvalue, comment, &status);
    if (status) {
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed read FITS header keyword no. "<<i<<"!");
			#endif
			HandleError(status,fp);
			return -1;
		}
		records.insert( std::make_pair<std::string,std::string>(keyname,keyvalue) );
 	}//end loop keywords

	if(records.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("FITS header keywords list is empty!");
		#endif
		return -1;
	}

 
	//##### GET STANDARD & MANDATORY KEYWORDS (if not existing set an invalid header...)
	// Get number of axis
	int Nchannels= 0;
	if ( records.find("NAXIS") == records.end() ) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid header detected (no NAXIS keyword)!");
		#endif
		return -1;
	}
	Nchannels= std::stoi(records["NAXIS"]);
		
	if(Nchannels<2){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of channels ("<<Nchannels<<"), at least 2 channels must be present!");
		#endif
		return -1;	
	}


	//Get dimensions
	long int Nx= 0;
	long int Ny= 0;
	if( records.find("NAXIS1") == records.end() || records.find("NAXIS2") == records.end() ) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid header detected (no NAXIS1/NAXIS2 keywords found!");
		#endif
		return -1;
	}
	Nx= std::stol(records["NAXIS1"]);//image size X
	Ny= std::stol(records["NAXIS2"]);//image size Y
	
	//Check for degenerate axis
	if(Nchannels>2){
		for(int i=3;i<Nchannels;i++){
			std::stringstream ss;
			ss<<"NAXIS"<<i;
			if( records.find(ss.str().c_str()) == records.end() ){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Cannot find keyword "<<ss.str()<<"!");
				#endif
				return -1;
			}
			long int NpixDegAxis= std::stol(records[ss.str().c_str()]);//pixels in current degenerate axis
			if(NpixDegAxis>1 || NpixDegAxis<=0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Axis "<<ss.str()<<" is invalid or not degenerate (N="<<NpixDegAxis<<"), cube or hypercube are not supported!");
				#endif
				return -1;
			}
		}//end loop axis
	}//close if degenerate axis check
	

	// Get pixel content unit
	std::string BUnit= "";
	if(records.find("BUNIT") == records.end()) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("BUNIT keyword not found in header!");
		#endif
	}
	else{
		BUnit= records["BUNIT"];	
		if(BUnit[0]=='\'') BUnit.erase(0,1);
		if(BUnit[BUnit.length()-1]=='\'') BUnit.erase(BUnit.length()-1,1);		
	}

	//Get Coordinate type and projection
	std::string CoordTypeX= "";
	std::string CoordTypeY= "";
	if(records.find("CTYPE1") == records.end() || records.find("CTYPE2") == records.end()) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("CTYPE keyword not found in header!");
		#endif
	}
	else{
		CoordTypeX= records["CTYPE1"];
		CoordTypeY= records["CTYPE2"];

		if(CoordTypeX[0]=='\'') CoordTypeX.erase(0,1);
		if(CoordTypeX[CoordTypeX.length()-1]=='\'') CoordTypeX.erase(CoordTypeX.length()-1,1);
		if(CoordTypeY[0]=='\'') CoordTypeY.erase(0,1);
		if(CoordTypeY[CoordTypeY.length()-1]=='\'') CoordTypeY.erase(CoordTypeY.length()-1,1);
	}

	//Get reference pixel id
	double CenterPixIdX= 0;
	double CenterPixIdY= 0;
	if(records.find("CRPIX1") == records.end() || records.find("CRPIX2") == records.end()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("CRPIX keyword(s) not found in header!");
		#endif
	}
	else{
		CenterPixIdX= std::stod(records["CRPIX1"]);
		CenterPixIdY= std::stod(records["CRPIX2"]);
	}
	
	//Get reference pixel coordinates (RA/DEC)
	double CenterPixX= 0;
	double CenterPixY= 0;
	if(records.find("CRVAL1") == records.end() || records.find("CRVAL2") == records.end()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("CRVAL keyword(s) not found in header!");
		#endif
	}
	else{
		CenterPixX= std::stod(records["CRVAL1"]);
		CenterPixY= std::stod(records["CRVAL2"]);
	}
	
	double PixStepX= 0;
	double PixStepY= 0;
	if(records.find("CDELT1") == records.end() || records.find("CDELT2") == records.end()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("CDELT keyword(s) not found in header!");
		#endif
	}
	else{
		PixStepX= std::stod(records["CDELT1"]);
		PixStepY= std::stod(records["CDELT2"]);
	}

	//Beam information
	double Bmaj= 0;
	double Bmin= 0;
	double Bpa= 0;
	if(records.find("BMAJ") == records.end() || records.find("BMIN") == records.end()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("BMAJ/BMIN keywords not found in header, beam info not available...");
		#endif
	}
	else{
		Bmaj= std::stod(records["BMAJ"]); 
		Bmin= std::stod(records["BMIN"]); 
	}

	if(records.find("BPA") == records.end()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("BPA keyword not found in header, setting it to zero!");
		#endif
		Bpa= 0;
	}
	else{
		Bpa= std::stod(records["BPA"]);
	}
	
	
	//########  GET NON STANDARD/NON-MANDATORY KEYWORDS
	// Get rotation matrix
	double RotX= 0;
	double RotY= 0;
	if(records.find("CROTA1") == records.end() || records.find("CROTA2") == records.end()){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("CROTA keyword(s) not found in header, setting rotation to 0!");
		#endif
	}
	else{
		RotX= std::stod(records["CROTA1"]);
		RotY= std::stod(records["CROTA2"]);
	}

	// Get Observation Location
	double ObsRA= -999;
	double ObsDEC= -999;
	if(records.find("OBSRA") == records.end() || records.find("OBSDEC") == records.end()){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("OBSRA/OBSDEC keywords not found in header, setting them to fake values!");
		#endif
	}
	else{
		ObsRA= std::stod(records["OBSRA"]);
		ObsDEC= std::stod(records["OBSDEC"]);
	}

	// Get epoch information
	double Epoch= 2000;
	if(records.find("EPOCH") == records.end()){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("EPOCH keyword not found in header, setting it to 2000!");
		#endif
	}
	else{
		Epoch= std::stod(records["EPOCH"]);
	}

	//## Get frequency info
	std::string FreqUnit= "";
	double Freq= 0;
	double dFreq= 0;
	double FreqRef= 0;
	if( records.find("CUNIT3") == records.end() || 
			records.find("CRVAL3") == records.end() || 
			records.find("CDELT3") == records.end() || 
			records.find("CRPIX3") == records.end()
	)
	{
		#ifdef LOGGING_ENABLED
			WARN_LOG("CRVAL3/CDELT3/CRPIX3 keywords not found in header (frequency info not available)...");
		#endif
	}
	else{	
		FreqUnit= records["CUNIT3"];	
		if(FreqUnit[0]=='\'') FreqUnit.erase(0,1);
		if(FreqUnit[FreqUnit.length()-1]=='\'') FreqUnit.erase(FreqUnit.length()-1,1);	
		Freq= std::stod(records["CRVAL3"]); 
		dFreq= std::stod(records["CDELT3"]);
		FreqRef= std::stod(records["CRPIX3"]); 
	}

	//## Fill header info in FITS INFO struct
	(fits_info.header).Nx= Nx;
	(fits_info.header).Ny= Ny;
	(fits_info.header).ObsRA= ObsRA;
	(fits_info.header).ObsDEC= ObsDEC;
	(fits_info.header).BUnit= BUnit;
	(fits_info.header).CoordTypeX= CoordTypeX;
	(fits_info.header).CoordTypeY= CoordTypeY;
	(fits_info.header).Cx= CenterPixIdX;
	(fits_info.header).Cy= CenterPixIdY;
	(fits_info.header).Xc= CenterPixX;
	(fits_info.header).Yc= CenterPixY;
	(fits_info.header).dX= PixStepX;
	(fits_info.header).dY= PixStepY;
	(fits_info.header).Bmaj= Bmaj;
	(fits_info.header).Bmin= Bmin;
	(fits_info.header).Bpa= Bpa;
	(fits_info.header).RotX= RotX;
	(fits_info.header).RotY= RotY;
	(fits_info.header).Epoch= Epoch;
	(fits_info.header).FreqUnit= FreqUnit;
	(fits_info.header).Freq= Freq;
	(fits_info.header).dFreq= dFreq;
	(fits_info.header).FreqRef= FreqRef;
	
	return 0;

}//close ReadHeader()





}//close namespace

